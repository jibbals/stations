# File input output for ozone stations work 
#
########################
######## Imports #######
########################

import numpy as np
from datetime import datetime
from datetime import timedelta

# directory and file listing
from os import walk as _walk
# library to read single line only
import linecache as _linecache
# CSV Library for reading/writing CSV files
from pandas.io.parsers import read_csv as _preadcsv
import netCDF4 as nc
import csv
from glob import glob 
# local module, converts GC tau to datetimes
from tau_to_date import tau_to_date as ttd

########################
######## Globals #######
########################

#directory paths
_Datadir='./data/'
_sondesdir='./data/sondes/'
_sitenames = ['Davis','Macquarie','Melbourne']
_Datafiles=[ _Datadir+d+'.nc' for d in _sitenames]
_sondesdirs=[_sondesdir+s for s in _sitenames]
_trac_avgs=_Datadir+"GC/trac_avg/*.nc"

_event_type={0:"misc",1:"front",2:"cutoff"}

# method to find all the csv paths:
def _list_files():
    sondecsvs=[[],[],[]]
    # for each location
    for dind, path in enumerate(_sondesdirs):
        # look in each year folder and pull out the csvs
        for yearpath, years, blahfiles in _walk(path):
            if 'lowres' in yearpath: # ignore low res data
                continue
            for fpath, blahdirs, fnames in _walk(yearpath):
                for fname in fnames:
                    if 'lowres' in fpath: continue
                    if (fname[0] != '.') and (fname[-4:] == '.csv') :
                        sondecsvs[dind].append(fpath+'/'+fname)
        uniqs=set(sondecsvs[dind])
        sondecsvs[dind] = list(uniqs)
        sondecsvs[dind].sort()
    return(sondecsvs)
# davis, macquarie, melbourne data locations
_sondecsvs=_list_files()
# now _sondecsvs[0] is all the davis csv files.. etc

def save_to_netcdf(outfilename, dimsdict, arraydict, unitsdict):
    '''
    Takes a dict of arrays and a dict of units...
    save to netcdf
    '''
    fid = nc.Dataset(outfilename, 'w', format='NETCDF4')
    fid.description = "NetCDF4 File written by Jesse"
    
    # add dimensions, these will also be one dimensional variables
    for key in dimsdict.keys():
        fid.createDimension(key, dimsdict[key].shape)
        var = fid.createVariable(key, dimsdict[key].dtype, (key,))
        var[:] = dimsdict[key]
        var.units=unitsdict[key]
    # add variables
    for key in arraydict.keys():
        #writedata and add units attribute
        # easy way to match variable to it's dimensionals?
        var = fid.createVariable(key, 'f8', arraydict[key].shape)
        var[:] = arraydict[key]
        var.units=unitsdict[key]
        #var.setncatts({'units': unitsdict[key]})
    
    fid.close()
    print("saved to "+outfilename)

def read_GC_station(station):
    '''
    Read a (GEOS-CHEM) station file, and the met field TROPP
    Return dictionary:
        AIRDEN: molecules/cm3
        PSURF: hpa
        Tau: Hrs since sliced bread
        O3: ppb, (GEOS_CHEM SAYS PPBV)
        BXHEIGHT: m
        Pressure: hPa
        # non station file stuff(added here)
        TropopausePressure: hPa #(interpolated from .a3 met fields)
        O3Density: molecules/cm3
        O3TropColumn: molecules/cm2
        TropopauseAltitude: m
        Altitudes: m
        PMids: hPa
        PEdges: hPa
        Date: datetime structure of converted Taus
        Station: string with stn name and lat/lon
    '''
    path=_Datafiles[station]
    fh=nc.Dataset(path, mode='r')
    stndict=dict()
    for key in ['Tau','Pressure','O3','PSURF','BXHEIGHT','AIRDEN']:
        stndict[key]=fh.variables[key][...]
    ## Station name and lat/lon
    stationstr="%s [%3.2fN, %3.2fE]"%(fh.station, float(fh.latitude), float(fh.longitude))
    troppfile= _Datadir+fh.station.lower()+'_2x25GEOS5_TROPP.csv'
    fh.close()
    
    # Density Column = VMR * AIRDEN [ O3 molecules/cm3 ]
    stndict['O3Density']= stndict['O3'] * 1e-9 * stndict['AIRDEN']
    
    # Get Date from Tau
    stndict['Date'] = np.array(ttd(stndict['Tau']))
    stndict['Station'] = stationstr
    
    # PSURF is pressure at bottom of boxes
    Pressure_bottoms=stndict['PSURF'][:,:]
    n_t=len(stndict['Tau'])
    n_y=72
    
    # geometric pressure mid points
    pedges=np.zeros( [n_t, n_y+1] )
    pedges[:, 0:n_y]=Pressure_bottoms
    pmids = np.sqrt(pedges[:, 0:n_y] * pedges[:, 1:(n_y+1)]) # geometric mid point
    
    # Altitude mid points from boxheights
    Altitude_mids=np.cumsum(stndict['BXHEIGHT'],axis=1)-stndict['BXHEIGHT']/2.0
    
    # Read station tropopause level:
    with open(troppfile) as inf:
        tropps=np.array(list(csv.reader(inf)))
    tropTaus=tropps[:,0].astype(float) # trops are every 3 hours, so we will only want half
    tropPres=tropps[:,1].astype(float)
    TPP = np.interp(stndict['Tau'], tropTaus, tropPres)
    stndict['TropopausePressure']=TPP
    
    # Tropospheric O3 Column!!!!
    TVC = []
    Atrop=[]
    TPL = []
    # For each vertical column of data
    for i in range(n_t):
        # tpinds = tropospheric part of column
        tpinds = np.where(pedges[i,:] > TPP[i])[0]
        tpinds = tpinds[0:-1] # drop last edge
        
        # Which layer has the tropopause
        TPLi=tpinds[-1]+1
        TPL.append(TPLi) # tropopause level
        
        # Fraction of TP level which is tropospheric
        pb, pt= pedges[i,TPLi], pedges[i,TPLi+1]
        frac= (pb - TPP[i])/(pb-pt)
        assert (frac>0) & (frac < 1), 'frac is wrong, check indices of TROPP'
        
        ## Find Trop VC of O3
        # sum of (molecules/cm3 * height(cm)) in the troposphere
        TVCi=np.sum(stndict['O3Density'][i,tpinds]*stndict['BXHEIGHT'][i,tpinds]*100)
        TVCi=TVCi+ frac*stndict['O3Density'][i,TPLi]*stndict['BXHEIGHT'][i,TPLi]*100
        TVC.append(TVCi)
        
        ## Altitude of trop
        #
        Atropi=np.sum(stndict['BXHEIGHT'][i, tpinds])
        Atropi = Atropi+frac*stndict['BXHEIGHT'][i, TPLi]
        Atrop.append(Atropi)
        
    stndict['O3TropColumn']=np.array(TVC)
    stndict['TropopauseAltitude']=np.array(Atrop)
    stndict['TropopauseLevel']=np.array(TPL)
    # Add pressure info
    stndict['PEdges']=pedges
    stndict['PMids']=pmids
    stndict['Altitudes']=Altitude_mids
    return(stndict)

def read_GC_global(subset=[-180,-90,180,90], savetofile=None, readfromfile=None):
    '''
    Read the GEOS-CHEM dataset, created by running pncgen on the GC trac_avg files produced by UCX_2x25 updated run
    Optionally subset to some region
    Optionally write to new file (will be smaller than all the GC files)
    Optionally read from new file (written with save to file)
    Return dictionary:
        psurf: hpa at level surface
        time: Hrs since sliced bread
        altitudes: m
        O3ppb: ppb, (GEOS_CHEM SAYS PPBV)
        boxheight: m
        tppressure: hPa
        tpaltitude: m
        tpaltitudeforTropVC: m used to calc tropospheric VC
        tplevel: level of tropopause
        O3density: molecules/cm3
        O3tropVC: molecules/cm2
        airdensity: molecules/m3
        pmids: hPa
        pedges: hPa
        date: datetime structure of converted Taus
    '''
    if readfromfile is not None:
        # Structure is saved in this file:
        assert False, "IMPLEMENT THIS"
    
    files=glob(_trac_avgs)
    files.sort()
    data=dict()
    
    # Read in the data then close the file
    with nc.MFDataset(files) as fh:
        
        # simple dimensions to read in
        for key in ['time','latitude','longitude']:
            data[key]=fh.variables[key][...]
        
        # also grab and rename these things
        data['latbounds']=np.append(fh.variables['latitude_bounds'][:][:,0], fh.variables['latitude_bounds'][:][-1,1])
        data['lonbounds']=np.append(fh.variables['longitude_bounds'][:][:,0], fh.variables['longitude_bounds'][:][-1,1])
        data['O3ppb']=fh.variables['IJ-AVG-$_O3'][...] # ppb
        data['tppressure']=fh.variables['TR-PAUSE_TP-PRESS'][...] # hPa
        data['tpaltitude']=fh.variables['TR-PAUSE_TP-HGHT'][...]*1e3 # km -> m
        data['tplevel']=fh.variables['TR-PAUSE_TP-LEVEL'][...] # dimensionless
        data['boxheight'] = fh.variables['BXHGHT-$_BXHEIGHT'][...] # m
        data['airdensity']=fh.variables['BXHGHT-$_N(AIR)'][...] # molecules / m3
        data['psurf']=fh.variables['PEDGE-$_PSURF'][...] # hpa at bottom of each vertical level
        
        
    
    # dimensions
    n_t, n_z, n_y, n_x=len(data['time']), 72, len(data['latitude']), len(data['longitude'])
    
    # Density Column = VMR * AIRDEN [ O3 molecules/cm3 ]
    data['O3density']= data['O3ppb'] * 1e-9 * data['airdensity'] * 1e-6 # ppb -> vmr and /m3 -> /cm3
    
    # Get Date from Tau
    data['date'] = np.array(ttd(data['time']))
    
    # geometric pressure mid points
    pedges=np.append(data['psurf'],np.zeros([n_t,1,n_y,n_x]),axis=1)
    pmids = np.sqrt(pedges[:, 0:n_z,:,:] * pedges[:, 1:(n_z+1),:,:]) # geometric mid point
    data['pmids']=pmids
    data['pedges']=pedges
    
    # Altitude mid points from boxheights
    data['altitudes']=np.cumsum(data['boxheight'],axis=1)-data['boxheight']/2.0
    
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
        O3i=data['O3density'][i,...] # molecules/cm3
        bhi=data['boxheight'][i,...]*100 # metres -> centimeters
        
        # Fraction of TP level which is tropospheric
        # frac= tpi - tpi.astype(int)
        
        #starttime=datetime.now()
        # working out how to do this dimensionally was too hard, just loop over..
        for x in range(n_x):
            for y in range(n_y):
                inds=np.where(levs[:,y,x]<=tpi[y,x])[0]
                #TVC[i,y,x] = np.sum(O3i[inds,y,x]*bhi[inds,y,x],axis=0) # molecules/cm2
                TVC[i,y,x] = np.inner(O3i[inds,y,x], np.transpose(bhi[inds,y,x]))
                Atrop[i,y,x]=altm[i,inds[-1],y,x] # m
                ## Add fractional part?
        #endtime=datetime.now()
        #print("Took %8.2f seconds to get Trop VC for one month"%(endtime-starttime).total_seconds())
        # np.sum takes ~ 0.32 seconds
        # np.inner takes ~ 0.27 seconds
    # sanity check:
    nearzero=np.squeeze(data['tpaltitude'])-Atrop
    assert np.mean(nearzero) < 1000, "TP Altitude changing by more than 1km"
    data['O3tropVC']=TVC # molecs/cm2
    data['tpaltitudeforTropVC']=Atrop # m
    
    if savetofile is not None:
        assert False, 'IMPLEMENT THIS'
    
    return(data)
    
# Sonde data class, created to hold stuff
#
class sondes:
    '''
    Class for holding sondes information at particular site
    '''
    def __init__(self):
        # class variables
        self.name = ''
        self.lat = 0
        self.lon = 0
        self.alt = 0 # site altitude (m)
        
        # these have one dim
        self.dates = []
        self.edates= [] # event dates
        self.edatedeltas= [] # sonde datetime - event datetime
        self.etype = [] # event is either from a cutoff low, a front, or something else
        self.einds = [] # event indices
        self.tpo3 = [] #ozone tropopause (km)
        self.tplr = [] # lapse rate tropopause (km)
        self.tp = [] # minimum tropopause (km)
        self.tpinds = [] # index of min tropopause
        self.tropvc = [] # tropospheric vertical column molecs/cm2
        
        # the following will have 2 dims: [dates,heights]
        self.press = 0 # heights(hPa) 
        self.gph = 0 # heights(m)
        self.o3pp = 0 # partial pressure(mPa)
        self.o3ppbv = 0 # ppbv[dates,heights]
        self.o3dens = 0 # molecules/cm3
        self.rh = 0 # rel humidity
        self.temp = 0 # temperature(Celcius)
        # extras
        self.seasonalstds = 0 # when making seasonal averages store the stds
        
    def count(self):
        return len(self.dates)
    def len(self):
        return len(self.dates)
    def __repr__(self): # if callled within prompt do this
        return self.__str__()
    def __str__(self): # what gets called when printing this class
        return("sondes(class):\n name, lat, lon, alt(m), "+\
                "dates, press(hPa), gph(km), o3pp(mPa), o3ppbv, rh, temp(C)")
    def utc_offset(self):
        # map the inputs to the function blocks
        offsets = {"Melbourne":11, "Macquarie":11,"Davis":7}
        return(offsets[self.name])
    def local_dates(self):
        # return stored dates plus utc offset
        offset = timedelta(hours=self.utc_offset())
        return([x+offset for x in self.dates])
    def get_index(self, date):
        # return index of closest matching profile
        closest=sorted(self.dates, key=lambda d : abs (d-date))[0]
        di = self.dates.index(closest)
        return di
    def get_profile(self, date):
        # return (gph, o3ppbv, temp, rh, dates, tp, tplr, tpo3) profile
        if np.size(self.tp) < 2:
            self._set_tps()
        # find closest matching date
        closest = sorted(self.dates, key = lambda d : abs( d - date))[0]
        di = self.dates.index(closest)
        profile=(self.gph[di,:], self.o3ppbv[di,:], self.temp[di,:],self.rh[di,:], \
            self.dates[di], \
            self.tp[di], self.tplr[di], self.tpo3[di])
        return profile
    def plot_profile(self, date, ytop=14, xtop=130, size=(8,16)):
        import matplotlib.pyplot as plt
        prof=self.get_profile(date)
        date2=prof[4]
        xl,xr=0,xtop
        yl,yr=0,ytop
        yaxi=prof[0]/1e3
        
        fig=plt.figure(figsize=size)
        ax1= fig.add_subplot(111)
        ax1.plot(prof[1],yaxi,'k',linewidth=2.0)
        ax2 = ax1.twiny()
        ax2.plot(prof[2],yaxi,'r')
        ax1.plot([xl,xr],[prof[5],prof[5]],'--k')
        ax1.set_ylim(yl,yr)
        ax1.set_xlim(xl,xr)
        ax2.set_xlim(-55,25)
        #plt.yscale('log')
        
        ax1.set_ylabel('GPH (km)')
        ax2.set_xlabel('Temp (C)')
        ax1.set_xlabel('Ozone (ppbv)')

        title=self.name+' '+date2.strftime('%Y-%m-%d')
        ax1.set_title(title,x=0.5,y=0.93)
        return(fig)
    def _set_tps(self):
        
        polar = (np.abs(self.lat) > 60)

        #for each profile
        ns= len(self.dates)
        for si in np.arange(0,ns):
            ppbv=np.array(self.o3ppbv[si,:])
            Z=np.array(self.gph[si,:]) / 1e3
            tpo3=-1.0

            ## FIRST
            ## set the ozone tropopause
            dmrdz = (ppbv[0:-1]-ppbv[1:])/(Z[0:-1]-Z[1:])
            
            
            # for each dvmr/dz gt 60 with ppbv > 80 check the .5 to 2 k higher values are gt 110
            # check ranges:
            z1=np.where(dmrdz>60)[0]
            z2=np.where(ppbv[0:-2] > 80)[0]
            testrange=np.intersect1d(z1,z2)
            upper = [ 2.0, 1.5 ][polar]
            if np.size(testrange) < 2 or np.size(Z) < 2 :
                self.tpo3.append(-1)
                self.tplr.append(-1)
                continue
            for ind in testrange:
                
                alt=Z[ind]
                #print "shape(Z):", np.shape(Z)
                #print "ind:", ind
                #print "testrange:", testrange
                #print "alt:", alt
                if np.isnan(alt) : continue
                z1=np.where(Z >= (alt+0.5))[0]
                z2=np.where(Z <= (alt+upper))[0]
                checks = np.intersect1d(z1,z2)
                #print "checks:", checks
                checks=list(checks)
                ##if all the indices in our check range are gt 110 ppbv we are finished
                if ((ppbv[checks] > 110 ).all()):
                    tpo3=alt
                    break
            self.tpo3.append(tpo3)
            
            ## SECOND 
            ## find temperature tropopause
            tplr=-1.0
            rate=-2.0
            minh=2.0
            temp=np.array(self.temp[si,:])
            temp=temp[Z > minh]
            Z=Z[Z>minh]
            lapse = (temp[0:-1]-temp[1:])/(Z[0:-1]-Z[1:])
            lapse = np.append(lapse,0)
            # lapse rate should be greater than -2 over two kilometers
            testrange=np.where(lapse > rate)[0]
            for ind in testrange[0:-1]:
                alt=Z[ind]
                z1=np.where(Z > alt)[0]
                z2=np.where(Z < (alt+2.0))[0]
                checks =np.intersect1d(z1,z2) 
                
                if np.mean(lapse[checks]) > rate :
                    tplr=alt
                    break
            self.tplr.append(tplr)
        ## FINALLY
        # tp is minimum of lapse rate and ozone tropopause
        self.tp = np.minimum(self.tplr,self.tpo3).tolist()
        for i in range(ns):
            Z=np.array(self.gph[i,:])/1e3
            if self.tp[i]  == -1:
                self.tpinds.append(np.NaN)
            else:
                self.tpinds.append(np.where(Z == self.tp[i])[0][0])

    def get_seasonal(self):
        '''
        Return a sondes class which has the average seasonal data for this class
        Not Yet Implemented
        '''
        
        # need to rebin everything using a standard height array
        #
        new_profile_points=1000
        
        ndat=sondes()
        ndat.name=self.name
        ndat.lat=self.lat
        ndat.lon=self.lon
        ndat.alt=self.alt
        ndat.dates=['J','F','M','A','M','J','J','A','S','O','N','D']
        ndat.o3pp = np.ndarray([12,new_profile_points])
        
        # create seasonal averages of other variables 
        #for each month
        for i in range(1,13): # 1..12
            #indexes of month i
            mi=self.month_indices(i)
              
            #get average partial pressure
            # for now just surface pp
            ndat.o3pp[i,0]= np.nanmean(self.o3pp[mi,0])
            
            # average of tps:
            #ndat.tps[i] = np.nanmean(tps[mi])
            
    
        return(ndat)
    
    def _set_density(self):
        '''
        PV=nRT:
        n/V [moles/cm3] = P/RT
        n/V * N_av [molecules/cm3] = P*N_av/RT
        n_O3 / V * N_av = P_O3 * N_av/RT
        
        P is pressure [hPa]
        P_O3 = vmr * P
        N_av = avocado's number [molecules/mole]
        R = gas constant [ cm3 hPa / ( K mole ) ]
        T = Temperature [K]
        
        Tropospheric VC= sum up to TP of (density * height)
        '''
        R=8.3144621e4 # [ cm3 hPa / K / mole ]
        N_av=6.0221409e+23 # Avegadro
        vmr=self.o3ppbv * 1e-9 # ppb to vmr
        P = np.array(self.press) # hPa
        T = np.array(self.temp) + 273.15 # celcius to Kelvin        
        self.o3dens = vmr * P / R / T * N_av
        
        # vertical column profile [molecs/cm2] = dens[molecs/cm3]*height[cm]
        Z=self.gph * 100
        edgesshape=list(Z.shape)
        edgesshape[1]+=1
        Zedges=np.zeros(edgesshape)
        Zedges[:,0] = Zedges[:,0] + self.alt*100
        Zedges[:,1:-1] = (Z[:, 0:-1] + Z[:, 1:]) / 2.0
        Zedges[:,-1] =Z[:,-1]
        # height of each vertical level
        Zheights= Zedges[:,1:]-Zedges[:, 0:-1] 
        vcp = self.o3dens * Zheights
        for i in range(len(self.dates)):
            tpind=self.tpinds[i]
            if tpind is np.NaN:
                self.tropvc.append(np.NaN)
            else:
                # doesn't include tropopause layer
                self.tropvc.append( np.nansum( vcp[i, 0:tpind] ) )
        
    def _set_events(self):
        '''
        return indices of dates matching event days in the CSVs (created from IDL)
        In this case the csv's are close to UTC+0 I think
        '''
        # read CSVs
        import csv
        csvs={"Melbourne":"events_melb.csv",
              "Macquarie":"events_mac.csv",
              "Davis":"events_dav.csv"}

        filename="data/"+csvs[self.name]
        with open(filename, 'rb') as csvfile:
            reader = csv.reader(csvfile)
            for row in reader:
                ir=[int(x) for x in row] # change strings to ints
                edate=datetime(ir[0],ir[1],ir[2],ir[3])
                # sort by difference between csv date and sonde date
                closest = sorted(self.dates, key = lambda d : abs( d - edate))[0]
                # add closest date index
                cind=self.dates.index(closest)
                self.einds.append(cind)
                self.edates.append(edate)
                self.edatedeltas.append(self.dates[cind].now() - edate.now())
    
    def degrade_vertically(self, gc):
        '''
        Degrade data to match GEOS Chem vertical levels
        input: gc
            station which get's returned from read_GC_station
        '''
        newinds=[]
        newppbv =[]
        newgph  =[]
        # For each sonde date, if there is a matching GEOS-Chem profile then we degrade the sonde
        for sind, date in enumerate(self.dates):
            gcdates=gc['Date']
            ymd = np.array([ d.strftime('%Y%m%d') for d in gcdates])
            ind=np.where( ymd=='%4d%02d%02d'%(date.year,date.month,date.day) )[0]
            
            # if no matching model data then skip date
            if len(ind) == 0:
                continue
            ind=ind[0]
            
            # Degrade sonde to GEOS vertical resolution
            # first what are the GEOS edges
            GCedges=np.zeros(73)
            GCedges[1:]=np.cumsum(gc['BXHEIGHT'][ind,:]) # in metres
            GCmids=gc['Altitudes'][ind,:] # in metres
            salts=self.gph[sind,:] # in metres
            prof = np.zeros(72)
            #degrade prof to average within GC box
            for i in range(72):
                prof[i]=np.mean(self.o3ppbv[sind, (GCedges[i] < salts) * ( salts < GCedges[i+1]) ])
            newinds.append(sind)
            newppbv.append(prof)
            newgph.append(GCmids)
        
        # update sonde where simple to do so
        self.dates = [self.dates[i] for i in newinds]
        self.edates= [] # event dates
        self.edatedeltas= [] # sonde datetime - event datetime
        self.einds = [] # event indices
        self.tpo3 = [self.tpo3[i] for i in newinds] #ozone tropopause (km)
        self.tplr = [self.tplr[i] for i in newinds]# lapse rate tropopause (km)
        self.tp = [self.tp[i] for i in newinds] # minimum tropopause (km)
        self.tpinds = [self.tpinds[i] for i in newinds] # index of min tropopause
        self.tropvc = [self.tropvc[i] for i in newinds] # tropospheric vertical column molecs/cm2
        
        # the following will have 2 dims: [dates,heights]
        self.press = 0 # heights(hPa) 
        self.gph = np.array(newgph) # heights(m)
        self.o3pp = 0 # partial pressure(mPa)
        self.o3ppbv = np.array(newppbv) # ppbv[dates,heights]
        self.o3dens = 0 # molecules/cm3
        self.rh = 0 # rel humidity
        self.temp = 0 # temperature(Celcius)
        
        

# end of header
_endofheader="#PROFILE"
_headerlines=51

########################
######## METHODS #######
########################
def _read_file(fpath):
    '''
    return data for single file
    RETURNS: (date, press, o3pp, temp, gph, rh)
    '''
    
    # date=datetime(y,m,d,h) at utc+0
    # first work out date:
    timeline=_linecache.getline(fpath, 26)
    date=0
    try:    
        (y,m,d) = map(int, timeline.split(',')[1].split('-'))
        h = int(timeline.split(',')[2].split(':')[0])
        date=datetime(y,m,d,h)    
    except  IndexError as ieerr:
        print("Index Error:"+ieerr.message)
        print(fpath)
        raise
    
    with open(fpath) as f:
        # 10 columns: [press,o3pp,temp,windspeed,winddir,levelcode,duration,gph,rh,sampletemp]
        usecols=[0,1,2,7,8]
        frame=_preadcsv(f, sep=',', header=_headerlines, usecols=usecols)
        cols=frame.columns.values
        press= frame[cols[0]].values
        o3pp = frame[cols[1]].values
        temp = frame[cols[2]].values
        gph  = frame[cols[3]].values
        rh   = frame[cols[4]].values
    
    return [date, press, o3pp, temp, gph, rh]
    
def ozonetp(ppbv,gph,polar=0):
    # Not Implemented
    return(0)
    
def read_site(site=0):
    '''
    read data from sondes csv files
    PARAMETERS: site=0
        site: 0 for davis, 1 for macca, 2 for maqcuarie
    '''
    assert(0<=site<3)
    # up to how many datapoints in a profile
    _profile_points = 1500
    
    # array to build up with data
    sondescount=len(_sondecsvs[site])
    sondesdata=np.ndarray([5,sondescount,_profile_points])+np.NaN
    
    # location from csv:
    locline=_linecache.getline(_sondecsvs[site][0], 22)
    lat,lon,alt = map(float,locline.split(','))
    
    dates=[]
    # loop through files appending data to lists
    for (fnum,fpath) in enumerate(_sondecsvs[site]):
        #(date, press, o3pp, temp, gph, rh)
        filedata=_read_file(fpath)
        dates.append(filedata[0])
        levels=len(filedata[1])
        for i in range(0,5):
            sondesdata[i,fnum,0:levels] = filedata[i+1]
        
    snd = sondes()
    snd.lat=lat
    snd.lon=lon
    snd.alt=alt
    snd.name=_sitenames[site]
    snd.dates=dates
    snd.press=sondesdata[0]
    snd.o3pp=sondesdata[1]
    snd.temp=sondesdata[2]
    snd.gph=sondesdata[3]
    snd.rh=sondesdata[4]
    snd.o3ppbv = snd.o3pp/snd.press*1.0e-5 *1.0e9
    snd._set_tps()
    snd._set_density()
    snd._set_events()
    return(snd)
