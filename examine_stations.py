# use python to examine dataset
# Created 20160616 by three gunmen in a salon

## modules
#

# These stop python from displaying images, faster to save images this way
import matplotlib
matplotlib.use('Agg')

# plotting, reading ncdf, csv, maths
import matplotlib.pyplot as plt
import netCDF4 as nc
import csv
import numpy as np

# these two help plot over datetimes
#from matplotlib.dates import date2num as d2n
from matplotlib.dates import DateFormatter
from datetime import datetime

# local module, converts GC tau to datetimes
from tau_to_date import tau_to_date as ttd

# Local module for reading sonde dataset
import fio

# set all fonts to 15
matplotlib.rcParams.update({'font.size': 15})

_Datadir='./data/'
_Datafiles=[ _Datadir+d for d in ['Davis.nc','Macquarie.nc', 'Melbourne.nc']]

# directory and file listing
_sondesdir='./data/sondes/'
_sondesdirs=[_sondesdir+ site for site in ['Davis/','Macquarie/','Melbourne/']]
_sitenames = ['Davis','Macquarie','Melbourne']

def read_station(station):
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

def event_profiles(station=2, data=None):
    '''
    Get O3, altitude, pressure profile at specific date.
    Returns: O3, Altitude_mids, Pressure_mids, Pressure_edges
    
    '''
    ## First read the station data
    if data is None:
        data=read_station(station)
    stn_name=data['Station'].split(' ')[0]
    dates=data['Date']
    
    ## Read the list of event dates from csv, created by blah.pro
    sondes= fio.read_site(station)
    SO3s  = sondes.o3ppbv # [time, altitude]
    SAlts = sondes.gph/1000 # [time, altitude] GPH in km
    eventdates=sondes.edates # this info is read from a csv I create in IDL
    
    ## At each event date plot the profiles side by side.
    for date in eventdates:
        outfile='images/eventprofiles/%s/%s_%s.ps'%(stn_name,stn_name, date.strftime("%Y%m%d"))
        
        # find matching model profile
        ymd = np.array([ d.strftime('%Y%m%d') for d in dates])
        ind=np.where( ymd=='%4d%02d%02d'%(date.year,date.month,date.day) )[0]
        
        if len(ind) == 0:
            outfile='images/eventprofiles/%s/missing_%s_%s.ps'%(stn_name,stn_name, date.strftime("%Y%m%d"))
            outf=open(outfile,'w')
            print('Missing %s'%date.strftime('%Y%m%d'))
            outf.close()
            continue
        
        # pick out midday hour TODO:
        #
        ind=ind[0]
        
        O3 = data['O3'][ind,:]
        Altitude=data['Altitudes'][ind,:] *1e-3 # m to km
        
        # plot the modelled event
        f=plt.figure(1,figsize=(6,10))
        plt.plot(O3,Altitude,color='r', 
                 linewidth=3, label='Modelled Profile')
        plt.ylabel('Altitude(km)')
        plt.xlabel('O3 (ppbv)')
        plt.xlim([5,120])
        plt.ylim([0,14])
        
        ## event profile on same plot
        Sind  = sondes.get_index(date)
        plt.plot(SO3s[Sind,:], SAlts[Sind,:], color='k',
                 linewidth=3, label='Sonde Profile')
        plt.legend()
        plt.savefig(outfile)
        print("plotted %s"%date.strftime('%Y%m%d'))
        plt.close(f)    
    
def time_series(outfile='images/StationSeries.png'):
    '''
    Plot timeseries for each station
    TODO: Also show dobson units...
    '''
    f3, f3axes = plt.subplots(3, 1, sharex=True, figsize=(14,10))
    f3axes[2].set_xlabel('Date')
    f3axes[1].set_ylabel('Ozone molecules/cm2')
    
    # read the three files into a list
    files=[ read_station(f) for f in range(3) ]
    
    # for each station do this
    for fh, f3ax, i in zip(files, f3axes, range(len(files))):
        ## grab variables
        
        # Ozone data array [time] in molecules/cm2
        data=fh['O3TropColumn']
        
        # create string title and turn tau's into matplotlib date numbers
        station=fh['Station']
        dates = fh['Date']
        
        # plot time series for each station
        f3ax.plot_date(dates, data[:])
        f3ax.set_title(station)
        
        # add dobson units
        newax=f3ax.twinx()
        if i == 1: newax.set_ylabel('Dobson Units')
        newylim= [1/2.69e16 * oldlim for oldlim in f3ax.get_ylim()]
        newax.set_ylim(newylim)
    
    # set plot titles
    f3.suptitle('Tropospheric ozone column from GEOS-Chem',fontsize=21)
    
    # set xticks nicely
    f3axes[2].xaxis.set_major_formatter( DateFormatter('%Y-%m') )
    f3axes[2].autoscale_view()
    f3.autofmt_xdate()
    
    # save then close plots
    #plt.show()
    f3.savefig(outfile)
    print("Image saved to %s"%outfile)
    plt.close(f3)

def yearly_cycle(hour=None, outfile='images/Yearly_cycle.png'):
    '''
    Show yearly cycle of surface ozone for our three sites
    If hour is set, only consider values at that hour (0, 6, 12, or 18)
    '''
    # read site data
    sites = [ read_station(p) for p in range(3) ]
    
    # Set up plot
    f = plt.figure(figsize=(16,9))
    f.suptitle("Seasonal cycle for tropospheric ozone from GEOS-Chem")
    plt.xlabel('Month')
    X=range(12)
    plt.xticks(X,['J','F','M','A','M','J','J','A','S','O','N','D'])
    plt.xlim([-0.5, 11.5])
    plt.ylabel('ozone molecules/cm2')
    colors=['magenta','cyan','orange']
    #for each station do this
    for site, color in zip(sites, colors):
        
        # Grab surface Ozone for all times
        O3 = site['O3TropColumn'][:]
        
        means=np.zeros(12)
        stds =np.zeros(12)
        # bin data into 12 months
        for month in range(12):
            allmonths=np.array([ d.month for d in site['Date'] ])
            inds=np.where(allmonths == month+1)[0]
            if hour is not None:
                allhours =np.array([ d.hour for d in site['Date'] ])
                inds = np.where(allmonths == month+1 & allhours==hour)[0]
            means[month]=np.mean(O3[inds])
            stds[month] =np.std(O3[inds])
        
        # plot the mean cycle and shade the area
        plt.plot(X, means, color=color, linewidth=3, label=site['Station'])
        plt.fill_between(X, means+stds, means-stds, alpha=0.3, color=color)
        
    leg=plt.legend(loc='best')
    for l in leg.legendHandles:
        l.set_linewidth(10)
    
    # add dobson units
    oldlims=plt.axes().get_ylim()
    newax=plt.twinx()
    newax.set_ylabel('Dobson Units')
    newylim= [1/2.69e16 * oldlim for oldlim in oldlims]
    newax.set_ylim(newylim)
    plt.savefig(outfile)
    print("Image saved to %s"%outfile)
    plt.close(f)

def monthly_profiles(hour=None):
    '''
    Profile mean and std deviation for each month for each site
    If you only want to consider a particular hour then set the hour parameter
        hour can be one of [0, 6, 12, 18, None]
    '''
    
    # read site data
    sites = [ read_station(p) for p in range(3) ]
    
    #for each station do this
    for site in sites:
        
        # Set up plot
        f, axes = plt.subplots(4,3, sharex=True, sharey=True, figsize=(16,16))
        axes[3,1].set_xlabel('Ozone (ppb)')
        xlim=[0,125]
        axes[3,1].set_xlim(xlim)
        axes[1,0].set_ylabel('Altitude (km)')
        ylim=[0,14]
        axes[1,0].set_ylim(ylim)
        
        # Grab Ozone
        O3 = site['O3']
        # metres to kilometres
        TP = site['TropopauseAltitude'] / 1000.0
        Z  = site['Altitudes']/1000.0 
        Znew= np.linspace(0,14,100)
        # need to vertically bin the O3 profiles,
        # interpolated at 100 points up to 14km
        means=np.zeros([12,100])
        stds =np.zeros([12,100])
        TPm = np.zeros(12)
        TPs = np.zeros(12)
        
        # bin data into 12 months
        allmonths=np.array([ d.month for d in site['Date'] ])
        for month in range(12):            
            inds=np.where(allmonths == month+1)[0]
            if hour is not None:
                allhours =np.array([ d.hour for d in site['Date'] ])
                inds = np.where( (allmonths == month+1) * (allhours==hour) )[0]
            # each profile needs to be interpolated up to 14km
            profs=np.zeros([len(inds),100])
            for i in range(len(inds)):
                profs[i,:] = np.interp(Znew, Z[inds[i],:], O3[inds[i],:])
            means[month,:]=np.mean(profs,axis=0)
            stds[month,:] =np.std(profs,axis=0)
            TPm[month] = np.mean(TP[inds])
            TPs[month] = np.std(TP[inds])
        
        # plot the mean profiles and shade the area of 1 stdev
        titles=np.array([['Dec','Jan','Feb'],['Mar','Apr','May'],['Jun','Jul','Aug'],['Sep','Oct','Nov']])
        months=np.array([[11,0,1],[2,3,4],[5,6,7],[8,9,10]])
        # colour each season!
        colours=['red','magenta','blue','green']
        for i in range(4):
            for j in range(3):
                X=means[months[i,j],:]
                Xl=X-stds[months[i,j],:]
                Xr=X+stds[months[i,j],:]
                axes[i,j].plot(X, Znew , linewidth=3, color='k')
                axes[i,j].fill_betweenx(Znew, Xl, Xr, alpha=0.3, color=colours[i])
                axes[i,j].set_title(titles[i,j])
                Y=TPm[months[i,j]]
                Ye=TPs[months[i,j]]
                axes[i,j].plot(xlim, [Y, Y], 'k--', linewidth=2 )
                axes[i,j].fill_between(xlim, [Y+Ye,Y+Ye], [Y-Ye,Y-Ye], alpha=0.2, color='k')
        
        # set title, and layout, then save figure
        stn_name=site['Station'].split(' ')[0]
        if hour is not None: stn_name+='_H%02d'%hour
        f.suptitle("Monthly Averaged Profiles from GEOS-Chem over "+stn_name)
        outfile='images/%s_profiles.png'%stn_name
        plt.tight_layout()
        plt.subplots_adjust(top=0.94)
        plt.savefig(outfile)
        print("Image saved to %s"%outfile)
        plt.close(f)

if __name__ == "__main__":
    print ("Running")
    
    event_profiles()
    #time_series()
    #yearly_cycle()
    #monthly_profiles()
    #[monthly_profiles(hour=h) for h in [0,6,12,18] ]
    
    
