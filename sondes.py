# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 09:03:57 2016

@author: jesse
"""

import numpy as np
from datetime import datetime, timedelta
# read CSVs
import csv

_event_type={0:"misc",1:"front",2:"cutoff"}

##################################################################
###### Sonde data class, created to hold stuff  ##################
##################################################################
class sondes:
    '''
    Class for holding sondes information at particular site
    '''
    # Event Types:
    _event_types={0:'misc',1:'front',2:'cutoff'}
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
        self.edepth= [] # sonde depth in km, ported from IDL
        self.eflux = [] # EVENT FLUX FROM GEOS-CHEM IN MOLECULES/CM2
        self.einds = [] # event indices
        self.epeak = [] # event peak altitudes (km)
        self.etype = [] # event is either from a cutoff low, a front, or something else
        self.etropvc = [] # TROP VC CALCULATED IN IDL
        self.fireflagged=[] # Events could be due to biomass burning
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
    
    def month_indices(self, month):
        ''' Return indexes of data from a particular month '''
        inds= [ i for i,d in enumerate(self.dates) if d.month==month ]
        return inds
        #for i,d in enumerate(self.dates):
        #    if d.month == month:
        #        inds.append(i)
    
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
        
        csvs={"Melbourne":"events_melb.csv",
              "Macquarie":"events_mac.csv",
              "Davis":"events_dav.csv"}

        filename="data/"+csvs[self.name]
        with open(filename, 'rb') as csvfile:
            reader = csv.reader(csvfile)
            #YYYY	MM	DD	HH	tropozone	flux	peak	tp	fire	etype
            header=reader.next()
            for row in reader:
                # Y, M, D, H
                ir=[int(x) for x in row[0:4]] # change strings to ints
                edate=datetime(ir[0],ir[1],ir[2],ir[3])
                # sort by difference between csv date and sonde date
                closest = sorted(self.dates, key = lambda d : abs( d - edate))[0]
                # add closest date index
                cind=self.dates.index(closest)
                self.einds.append(cind)
                self.edates.append(edate)
                self.edatedeltas.append(self.dates[cind].now() - edate.now())
                
                # Read in the flux and trop column worked out in IDL getevents.pro
                # Tropospheric VC of ozone (molecules/cm2), flux (molecules/cm2)
                vcs = [float(x) for x in row[4:6]]
                self.etropvc.append(vcs[0])
                self.eflux.append(vcs[1])
                
                # read the event peak altitude (km)
                self.epeak.append(float(row[6]))
                self.edepth.append(float(row[7])-float(row[6]))
                
                # read weather the event may be due to biomass burning
                self.fireflagged.append(row[8] == "1")
                
                # finally grab the likely cause of the event (misc,front,cutoff)
                self.etype.append(int(row[9]))
    
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
        
        