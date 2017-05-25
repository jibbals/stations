# File input output for ozone stations work 
# uses python 2.7 'maths' environment
#
########################
######## Imports #######
########################

import numpy as np
from datetime import datetime

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

# data classes
from GChem import GChem, GCArea
from sondes import sondes

########################
######## Globals #######
########################

#directory paths
_Datadir='./data/'
_sondesdir='./data/sondes/'
_sitenames = ['Davis','Macquarie','Melbourne']
GC_stn_from_name={'Davis':0,'Macquarie':1,'Melbourne':2}
_Datafiles=[ _Datadir+d+'.nc' for d in _sitenames]
_sondesdirs=[_sondesdir+s for s in _sitenames]
_trac_avgs=_Datadir+"GC/trac_avg/*.nc"

# end of header
_endofheader="#PROFILE"
_headerlines=51

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

########################
######## METHODS #######
########################

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

def get_GC_area():
    ''' Read the monthly gridded AREA (m2) for GEOS-Chem. '''
    return GCArea()

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

def read_GC_global():
    '''
    Read the GEOS-CHEM dataset, created by running pncgen on the GC trac_avg files produced by UCX_2x25 updated run
    Optionally subset to some region
    Optionally write to new file (will be smaller than all the GC files)
    Optionally read from new file (written with save to file)
    Returns: GChem class
    '''
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
        data['airdensity']=fh.variables['BXHGHT-$_N(AIR)'][...]*1e-6 # molecules / m3 -> molecules/cm3
        data['psurf']=fh.variables['PEDGE-$_PSURF'][...] # hpa at bottom of each vertical level
    
    return(GChem(data))

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
    
def read_sonde(site=0):
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

def read_ANDREW():
    stts=np.zeros([4,12]) # 4 sites, 12 months
    sdic={'Davis':0,'Macquarie Island':1,'Melbourne':2,'Laverton':3}
    snames=['Davis','Macquarie Island','Melbourne','Laverton']
    with open('data/AndrewOzonesondeSTEProxy.csv') as fin:
        reader=csv.reader(fin)
        for line in reader:
            site,month,stt=line[0],int(line[1])-1,float(line[2])
            sind=sdic[site]
            stts[sind,month]=float(stt)
    return (snames, stts)
            
