# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 08:34:11 2016

Classes for GEOS-Chem Dataset and Sonde data set

Me das tu numero?
@author: jesse
"""

import numpy as np
from datetime import datetime#, timedelta
from tau_to_date import tau_to_date as ttd

##################################################################
###### GEOS-Chem, created to hold stuff         ##################
##################################################################
class GChem:
    '''
    Class to hold and access GC data
    '''
    def __init__(self, data):
        '''
        Initiallise using dictionary read in by fio.py
        '''
        
        # one dimensional things
        self.lats = data['latitude']            #  -90 to 90 degrees
        self.latbounds = data['latbounds']    
        self.lons = data['longitude']           # -180 to 180 degrees
        self.lonbounds = data['lonbounds']
        self.taus = data['time']                # Hours since 19850101 000000
        self.dates = np.array(ttd(self.taus))   # Array of datetimes
        
        # latlon indexes of davis, macca, melb
        self.siteinds=[(0,0),(0,0),(0,0)]
        
        # More Dimensional Things: [ time, lev, lat, lon]
        self.O3ppb = data['O3ppb']              # ppb
        self.airdensity = data['airdensity']    # molecules/cm3
        self.tppressure = data['tppressure']    # hPa
        self.tpaltitude = data['tpaltitude']    # m
        self.tplevel = data['tplevel']          # dimensionless
        self.boxheight = data['boxheight']      # m
        self.psurf = data['psurf']              # hpa at bottom of each vertical level
        
        n_t, n_z, n_y, n_x=len(data['time']), 72, len(data['latitude']), len(data['longitude'])
        self.dims=np.array([n_t,n_z,n_y,n_x])
        
        # geometric pressure mid points, and pressure edges: hPa
        self.pedges=np.append(data['psurf'],np.zeros([n_t,1,n_y,n_x]),axis=1)
        self.pmids = np.sqrt(self.pedges[:, 0:n_z,:,:] * self.pedges[:, 1:(n_z+1),:,:]) # geometric mid point
        
        # Density Column = VMR * AIRDEN [ O3 molecules/cm3 ]
        self.O3density = data['O3ppb'] * 1e-9 * data['airdensity'] # 1e-9*ppb = vmr
                
        # Altitude mid points from boxheights: m
        self.altitudes = np.cumsum(data['boxheight'],axis=1)-data['boxheight']/2.0
        
        # Tropospheric O3 Column!!!!
        TVC = np.ndarray([n_t,n_y,n_x]) + np.NaN
        Atrop=np.ndarray([n_t,n_y,n_x]) + np.NaN
        altm=np.cumsum(data['boxheight'], axis=1) # m
        # For each vertical column of data
        for i in range(n_t):
            tpi=data['tplevel'][i,0,:,:]
            
            # tpinds = tropospheric part of column
            levs=np.arange(1,72.5)
            levs=np.transpose(np.tile(levs[:], (n_x,n_y,1)))
            #inds=np.where(levs<tpi)
            
            # ozone and boxheight for the entire troposphere
            O3i=self.O3density[i,...] # molecules/cm3
            bhi=data['boxheight'][i,...]*100 # metres -> centimeters
            
            # working out how to do this dimensionally was too hard, just loop over..
            for x in range(n_x):
                for y in range(n_y):
                    inds=np.where(levs[:,y,x]<=tpi[y,x])[0]
                    #TVC[i,y,x] = np.sum(O3i[inds,y,x]*bhi[inds,y,x],axis=0) # molecules/cm2
                    TVC[i,y,x] = np.inner(O3i[inds,y,x], np.transpose(bhi[inds,y,x]))
                    Atrop[i,y,x]=altm[i,inds[-1],y,x] # m
                    ## Add fractional part?
            # np.sum takes ~ 0.32 seconds, while np.inner takes ~ 0.27 seconds
        # sanity check:
        nearzero=np.squeeze(data['tpaltitude'])-Atrop
        assert np.mean(nearzero) < 500, "TP Altitude changing by more than 500m"
        
        self.O3tropVC = TVC                 # molecs/cm2
        self.tpaltitudeforTropVC = Atrop    # m
        
    def count(self):
        return len(self.taus)
    def len(self):
        return len(self.taus)
    def mapData(self, data, label=None,lon0=-179,lat0=-80,lon1=179,lat1=80):
        from mpl_toolkits.basemap import Basemap
        
        lons,lats=np.meshgrid(self.lonbounds, self.latbounds)
        m=Basemap(lon0,lat0,lon1,lat1)
        m.pcolormesh(np.transpose(lons),np.transpose(lats),np.transpose(data), latlon=True)
        m.drawcoastlines()
        cb= m.colorbar()
        if label is not None:
            cb.set_label(label)
        return m,cb
    
    def averagedTVC_monthly(self, Region):
        '''
        Average tropospheric vertical column for Region
        Region: [S ,W ,N ,E]        
        Returns: data, date, std
        Return units: molecules/cm2
        Return dimension: X months of data
        '''
        allmonths= np.array([ d.month for d in self.dates ])
        allyears= np.array([d.year for d in self.dates])
        data=[]
        std=[]
        date=[]
        south,west,north,east=Region
        assert north > south, "North needs to be greater than south"
        assert east > west, "East needs to be greater than west"
        
        loninds=np.where( (east > self.lons) * (self.lons > west) )[0]
        latinds=np.where( (north > self.lats) * (self.lats > south) )[0]
        TVC=self.O3tropVC[:, latinds, :]
        TVC=TVC[:, :, loninds]
        for year in set(allyears):
            for month in range(12):
                inds= (allmonths == month+1) * (allyears==year)
                data.append(np.mean(TVC[inds,:,:]))
                std.append(np.std(TVC[inds,:,:]))
                date.append(datetime(year,month+1,1))
        return data, date, std
    
    def averagedTVC(self, Region, homeoskedastic=True):
        '''
        Average tropospheric vertical column for Region
        Region: [S ,W ,N ,E]        
        Returns: data, std
            data
            std. deviation
        Return units: molecules/cm2
        Return dimension: 12 months
        '''
        data=np.zeros(12)
        std=np.zeros(12)
        allmonths= np.array([ d.month for d in self.dates ])
        south,west,north,east=Region
        assert north > south, "North needs to be greater than south"
        assert east > west, "East needs to be greater than west"
        
        loninds=np.where( (east > self.lons) * (self.lons > west) )[0]
        latinds=np.where( (north > self.lats) * (self.lats > south) )[0]
        TVC=self.O3tropVC[:, latinds, :]
        TVC=TVC[:, :, loninds]
        for i in range(12):
            minds=np.where(allmonths == i+1)[0]
            data[i]=np.mean(TVC[minds])
            std[i]=np.std(TVC[minds])
        if homeoskedastic:
            std=np.std(TVC)
        return data, std
    
    def southernOceanTVC(self, north=-35,south=-75):
        '''
        Get the Tropospheric Vertical Column averaged into months
        By default averaged over 35S to 75S ( range includes our 3 data sites )
        Returns:
            molecules/cm2 [12 months] 
        '''
        data=np.zeros(12)
        allmonths= np.array([ d.month for d in self.dates ])
        
        latinds=np.where( (north > self.lats) * (self.lats > south) )[0]
        SOTVC=self.O3tropVC[:,latinds,:]
        for i in range(12):
            minds=np.where(allmonths == i+1)[0]
            data[i]=np.mean(SOTVC[minds])
        return data
    
    def saveToFile(self):
        '''
        Save the shit to a netcdf file
        '''
        assert False, "Not Implemented"
    
    def readFromFile(self, filename):
        '''
        Read from a specific file, written to using the SaveToFile() function
        '''
        assert False, "Not Implemented"
        return 0
    
    def subset(self, subset=[-180,-90,180,90]):
        '''
        Return copy of class subsetted to some space...
        '''
        assert False, "Not Implemented"
        return 0

##################################################################
###### GEOS-Chem area                           ##################
##################################################################
class GCArea:
    def __init__(self):
        fname="/home/jesse/Desktop/Repos/stations/data/GC/GC_area.nc"
        from netCDF4 import Dataset as dset 
        ncd=dset(fname,mode='r')
        self.lats=ncd['latitude'][...] # [91]
        self.lons=ncd['longitude'][...] # [144]
        self.latedges=ncd['latedges'][...]
        self.lonedges=ncd['lonedges'][...]
        self.area=ncd['area'][...] # square metres [144, 91]
        
    def band_area(self,lat0,lat1):
        inds= (self.lats >= lat0) * (self.lats <= lat1)
        #        print ('BUGGGGG')
        #        print (np.sum(self.area[inds,:])) # 7.1e13
        #        print (np.sum(self.area[:,inds])) # 1e14
        #        print ('Changed?')
        return np.sum(self.area[:,inds])
    
    def region_area(self, Region):
        S,W,N,E=Region
        lats = (self.lats >= S) * (self.lats <= N) 
        lons = (self.lons >= W) * (self.lons <= E)
        Area=self.area[lons,:]
        Area=Area[:,lats]
        return np.sum(Area)
    
        
        