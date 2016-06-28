# -*- coding: utf-8 -*-

########################
######## Imports #######
########################

import numpy as _np
from datetime import datetime
from datetime import timedelta

# directory and file listing
from os import walk as _walk
# library to read single line only
import linecache as _linecache
# CSV Library for reading/writing CSV files
from pandas.io.parsers import read_csv as _preadcsv

########################
######## Globals #######
########################
_sondesdir='./data/sondes/'
_sitenames = ['Davis','Macquarie','Melbourne']
_sondesdirs=[_sondesdir+s for s in _sitenames]

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
        sondecsvs[dind].sort()
    return(sondecsvs)
# davis, macquarie, melbourne data locations
_sondecsvs=_list_files()
# now _sondecsvs[0] is all the davis csv files.. etc

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
        self.einds = [] # event indices
        self.tpo3 = [] #ozone tropopause
        self.tplr = [] # lapse rate tropopause
        self.tp = [] # minimum tropopause
        
        # the following will have 2 dims: [dates,heights]
        self.press = 0 # heights(hPa) 
        self.gph = 0 # heights(m)
        self.o3pp = 0 # partial pressure(mPa)
        self.o3ppbv = 0 # ppbv[dates,heights]
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
        if _np.size(self.tp) < 2:
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
        
        polar = (_np.abs(self.lat) > 60)

        #for each profile
        ns= len(self.dates)
        for si in _np.arange(0,ns):
            ppbv=_np.array(self.o3ppbv[si,:])
            Z=_np.array(self.gph[si,:])/1e3
            tpo3=-1.0

            ## FIRST
            ## set the ozone tropopause
            dmrdz = (ppbv[0:-1]-ppbv[1:])/(Z[0:-1]-Z[1:])
            
            
            # for each dvmr/dz gt 60 with ppbv > 80 check the .5 to 2 k higher values are gt 110
            # check ranges:
            z1=_np.where(dmrdz>60)[0]
            z2=_np.where(ppbv[0:-2] > 80)[0]
            testrange=_np.intersect1d(z1,z2)
            upper = [ 2.0, 1.5 ][polar]
            if _np.size(testrange) < 2 or _np.size(Z) < 2 :
                self.tpo3.append(-1)
                self.tplr.append(-1)
                continue
            for ind in testrange:
                
                alt=Z[ind]
                #print "shape(Z):", _np.shape(Z)
                #print "ind:", ind
                #print "testrange:", testrange
                #print "alt:", alt
                if _np.isnan(alt) : continue
                z1=_np.where(Z >= (alt+0.5))[0]
                z2=_np.where(Z <= (alt+upper))[0]
                checks = _np.intersect1d(z1,z2)
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
            temp=_np.array(self.temp[si,:])
            temp=temp[Z > minh]
            Z=Z[Z>minh]
            lapse = (temp[0:-1]-temp[1:])/(Z[0:-1]-Z[1:])
            lapse = _np.append(lapse,0)
            # lapse rate should be greater than -2 over two kilometers
            testrange=_np.where(lapse > rate)[0]
            for ind in testrange[0:-1]:
                alt=Z[ind]
                z1=_np.where(Z > alt)[0]
                z2=_np.where(Z < (alt+2.0))[0]
                checks =_np.intersect1d(z1,z2) 
                
                if _np.mean(lapse[checks]) > rate :
                    tplr=alt
                    break
            self.tplr.append(tplr)
        ## FINALLY
        # tp is minimum of lapse rate and ozone tropopause
        self.tp = _np.minimum(self.tplr,self.tpo3).tolist()

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
        ndat.o3pp = _np.ndarray([12,new_profile_points])
        
        # create seasonal averages of other variables 
        #for each month
        for i in range(1,13): # 1..12
            #indexes of month i
            mi=self.month_indices(i)
              
            #get average partial pressure
            # for now just surface pp
            ndat.o3pp[i,0]= _np.nanmean(self.o3pp[mi,0])
            
            # average of tps:
            #ndat.tps[i] = _np.nanmean(tps[mi])
            
    
        return(ndat)
        
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
    sondesdata=_np.ndarray([5,sondescount,_profile_points])+_np.NaN
    
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
    snd._set_events()
    return(snd)
