# use python to examine dataset
# Created 20160616 by three gunmen in a salon

## modules
#

# These stop python from displaying images, faster to save images this way
import matplotlib
matplotlib.use('Agg')


# plotting, reading ncdf, csv, maths
import matplotlib.pyplot as plt # plotting module
import matplotlib.cm as cm # colour maps
from matplotlib.dates import DateFormatter
import numpy as np
from datetime import datetime

# Local module for reading sonde dataset
import fio

# set all fonts to 15
matplotlib.rcParams.update({'font.size': 15})

def anomaly_corellation(outfile='images/corellations_anomalies.png'):
    '''
    plot correlation of monthly average anomalies between sondes and GC where both occur on the same day.
    '''
    from scipy import stats
    
    # read the three files into a list
    files=[ fio.read_GC_station(f) for f in range(3) ]
    o3sondes = [fio.read_site(s) for s in range(3) ]
    
    f3, f3axes = plt.subplots(3, 1, figsize=(12,16))
    f3axes[2].set_xlabel('Sonde tropospheric O3 relative anomaly')#(molecules/cm2)')
    f3axes[1].set_ylabel('GEOS-Chem tropospheric O3 relative anomaly')#(molecules/cm2)')
    #xlims=[[1e17,1e18], [1e17,1e18], [1e17,1.5e18]]
    #ylims=[[1e17,1.5e18], [1e17,1.5e18], [1e17, 2e18]]
    
    # set up colorbar
    cmap=plt.get_cmap('rainbow', 12)
    
    # for each station do this
    # gc, os, f3ax, i = files[0], o3sondes[0], f3axes[0], range(len(files))[0]
    for gc, os, f3ax, i in zip(files, o3sondes, f3axes, range(len(files))):
        #xlim,ylim=xlims[i],ylims[i]
        xlim,ylim=[-1,1],[-1.5,1.5]
        plt.sca(f3ax) # set current axis
        ## grab variables
        # Ozone data array [time] in molecules/cm2
        # Note I only want UTC0 data from the model
        allhours=np.array([ d.hour for d in gc['Date'] ])
        H0inds=np.where(allhours==0)[0]
        
        data=gc['O3TropColumn'][H0inds]
        sdata=np.array(os.tropvc)
        
        # create string title and turn tau's into matplotlib date numbers
        station=gc['Station']
        dates = gc['Date'][H0inds]
        allmonths=np.array([ d.month for d in dates ])
        sdates= np.array(os.dates)
        sallmonths= np.array([ d.month for d in sdates ])
        mean = np.zeros(12)
        smean= np.zeros(12)
        for i in range(12):
            minds=np.where(allmonths == i+1)[0]
            sminds=np.where(sallmonths == i+1)[0]
            mean[i]=np.nanmean(data[minds])
            smean[i]=np.nanmean(sdata[sminds])
        # only look at data with matching dates
        osh0dates= [datetime(d.year,d.month,d.day,0,0) for d in sdates]
        Yos=[]
        Ygc=[]
        mnth=[]
        for si, sdate in enumerate(osh0dates):
            if sdate in dates:
                ind=np.where(dates == sdate)[0]
                if len(ind) == 0 or np.isnan(sdata[si]):
                    continue
                else: 
                    ind=ind[0]
                #print(ind, dates[ind])
                #print(si, sdates[si])
                #Work out anomaly and append to list
                mind=sdate.month-1
                Yos.append((sdata[si]-smean[mind])/smean[mind])
                Ygc.append((data[ind]-mean[mind])/mean[mind])
                mnth.append(sdate.month)
        
        
        Yos,Ygc = np.array(Yos),np.array(Ygc)
        # plot correlation coefficient
        slope,intercept,r_value,p_value,std_err= stats.linregress(Yos,Ygc)
        colors=cmap(mnth)
        f3ax.scatter(Yos, Ygc, color=colors, cmap=cmap, label='Trop O3')
        f3ax.plot(Yos, intercept+slope*Yos, 'k-', label='Regression')
        #f3ax.plot(xlim, xlim, 'k--', label='1-1 line')
        f3ax.plot([-1,1], [-1, 1], 'k--', label='1-1 line')
        f3ax.set_title(station)
        f3ax.set_ylim(ylim)
        f3ax.set_xlim(xlim)
        txty=ylim[0]+0.1*(ylim[1]-ylim[0])
        txtx=xlim[0]+0.81*(xlim[1]-xlim[0])
        txty2=ylim[0]+0.17*(ylim[1]-ylim[0])
        txty3=ylim[0]+.24*(ylim[1]-ylim[0])
        plt.text(txtx,txty,"N=%d"%len(Yos))
        plt.text(txtx,txty2,"r=%5.3f"%r_value)
        plt.text(txtx,txty3,"slope=%5.3f"%slope)
        if i==1: plt.legend()
    # set plot titles
    f3.suptitle('Relative anomaly from monthly mean',fontsize=21)
    
    # add colourbar space to the right
    f3.subplots_adjust(right=0.85)
    cbar_ax = f3.add_axes([0.9, 0.25, 0.04, 0.5])
    
    # add colourbar, force the stupid thing to be the same as the one used in plotting...
    sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(vmin=1,vmax=12))
    sm._A=[]
    cb=f3.colorbar(sm,cax=cbar_ax)
    
    # set the 12 ticks nicely and roughly centred
    cb.set_ticks(np.linspace(1,12,12) + (6.5-np.linspace(1,12,12))/12.)
    cb.set_ticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    #cb.set_label('month')
    
    # save then close plots
    #plt.show()
    f3.savefig(outfile)
    print("Image saved to %s"%outfile)
    plt.close(f3)

def corellation(outfile='images/corellations.png'):
    '''
    plot correlation between sondes and GC where both occur on the same day.
    '''
    from scipy import stats
    
    # read the three files into a list
    files=[ fio.read_GC_station(f) for f in range(3) ]
    o3sondes = [fio.read_site(s) for s in range(3) ]
    
    f3, f3axes = plt.subplots(3, 1, figsize=(12,16))
    f3axes[2].set_xlabel('Sonde tropospheric O3 (molecules/cm2)')
    f3axes[1].set_ylabel('GEOS-Chem tropospheric O3 (molecules/cm2)')
    xlims=[[1e17,1e18], [1e17,1e18], [1e17,1.5e18]]
    ylims=[[1e17,1.5e18], [1e17,1.5e18], [1e17, 2e18]]
    
    # set up colorbar
    cmap=plt.get_cmap('rainbow', 12)
    
    # for each station do this
    # gc, os, f3ax, i = files[0], o3sondes[0], f3axes[0], range(len(files))[0]
    for gc, os, f3ax, i in zip(files, o3sondes, f3axes, range(len(files))):
        xlim,ylim=xlims[i],ylims[i]
        plt.sca(f3ax) # set current axis
        ## grab variables
        # Ozone data array [time] in molecules/cm2
        data=gc['O3TropColumn']
        sdata=np.array(os.tropvc)
        
        # create string title and turn tau's into matplotlib date numbers
        station=gc['Station']
        dates = gc['Date']
        
        sdates= np.array(os.dates)
        
        # only look at data with matching dates
        osh0dates= [datetime(d.year,d.month,d.day,0,0) for d in sdates]
        Yos=[]
        Ygc=[]
        mnth=[]
        for sdate, si in zip(osh0dates, range(len(osh0dates))):
            if sdate in dates:
                ind=np.where(dates == sdate)[0]
                if len(ind) == 0 or np.isnan(sdata[si]):
                    continue
                else: 
                    ind=ind[0]
                #print(ind, dates[ind])
                #print(si, sdates[si])
                Yos.append(sdata[si])
                Ygc.append(data[ind])
                mnth.append(sdate.month)
        
        # Instead of data, look at correlation of difference from the monthly mean
        #TODO:
        Yos,Ygc = np.array(Yos),np.array(Ygc)
        # plot correlation coefficient
        slope,intercept,r_value,p_value,std_err= stats.linregress(Yos,Ygc)
        #f3ax.plot(Yos, Ygc, 'o', c=mnth, label='Trop O3')
        
        #mcolours= cm.rainbow(np.array(mnth)/12.0)
        colors=cmap(mnth)
        f3ax.scatter(Yos, Ygc, color=colors, cmap=cmap, label='Trop O3')
        f3ax.plot(Yos, intercept+slope*Yos, 'k-', label='Regression')
        f3ax.plot(xlim, xlim, 'k--', label='1-1 line')
        f3ax.set_title(station)
        f3ax.set_ylim(ylim)
        f3ax.set_xlim(xlim)
        txty=ylim[0]+0.1*(ylim[1]-ylim[0])
        txtx=xlim[0]+0.81*(xlim[1]-xlim[0])
        txty2=ylim[0]+0.17*(ylim[1]-ylim[0])
        txty3=ylim[0]+.24*(ylim[1]-ylim[0])
        plt.text(txtx,txty,"N=%d"%len(Yos))
        plt.text(txtx,txty2,"r=%5.3f"%r_value)
        plt.text(txtx,txty3,"slope=%5.3f"%slope)
        if i==1: plt.legend()
    # set plot titles
    f3.suptitle('Corellation',fontsize=21)
    
    # add colourbar space to the right
    f3.subplots_adjust(right=0.85)
    cbar_ax = f3.add_axes([0.9, 0.25, 0.04, 0.5])
    
    # add colourbar, force the stupid thing to be the same as the one used in plotting...
    sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(vmin=1,vmax=12))
    sm._A=[]
    cb=f3.colorbar(sm,cax=cbar_ax)
    
    # set the 12 ticks nicely and roughly centred
    cb.set_ticks(np.linspace(1,12,12) + (6.5-np.linspace(1,12,12))/12.)
    cb.set_ticklabels(['J','F','M','A','M','J','J','A','S','O','N','D'])
    #cb.set_label('month')
    
    # save then close plots
    #plt.show()
    f3.savefig(outfile)
    print("Image saved to %s"%outfile)
    plt.close(f3)



def SO_extrapolation():
    '''
    Rough estimate of Southern Oceanic STT flux
    '''
    # Read sondes data
    sonde_files=[fio.read_site(s) for s in range(2)]
    # monthly likelihood = occurrences/sondecount each month, for davis and macquarie
    likelihood=np.zeros([12,2])
    # monthly flux percentage
    fluxpct=np.zeros([12,2])
    # monthly modelled O3 in SO
    SOTropO3 = np.zeros(12)
    
    # for davis and macquarie:
    for os,i in zip(sonde_files,range(2)):
        # for each month
        for month in range(12):
            # indices of all events and all sondes within this month
            allmonths=np.array([ d.month for d in os.dates ])
            alleventsmonths=np.array([ d.month for d in os.edates ])
            inds=np.where(allmonths == month+1)[0]
            einds=np.where(alleventsmonths == month+1)[0]
            likelihood[month,i] = len(einds)/ float(len(inds))
            # TODO: add flux amount to events_site.csv, will need to do this in IDL
            fluxpct[month,i] = 0.03
    
    # model SO tropospheric O3 Column:
    for month in range(12):
        # TODO: get this data from GEOS-Chem output... try to do the bpch to coards or else more IDL
        SOTropO3[month] = 0.5e18 # molecules/cm2
    
    # plot estimated flux on left hand axis
    # plot likelihoods and flux pct on the right hand axis
    f=np.mean(fluxpct, axis=1)
    l=np.mean(likelihood, axis=1)
    flux = SOTropO3 * f * l
    
    plt.clf()
    ax=plt.subplot(111)
    X=range(12)
    plt.plot(X,flux, color='black', linewidth=3, label="Flux")
    
    # plot percentages
    newax=plt.twinx()
    newax.plot(X,f*100, color='teal', label='STT contribution')
    newax.plot(X,l*100, color='magenta', label='STT likelihood')
    
    # axes and legends
    newax.legend(loc=1)
    newax.set_ylabel('percent')
    newax.set_ylim([0,50])
    ax.legend(loc=2)
    ax.set_ylim([1e14, 6e15])
    ax.set_ylabel('Flux (molecules/cm2)')
    plt.xlim([-0.5, 11.5])
    ax.set_xlabel('Month')
    plt.xticks(X,['J','F','M','A','M','J','J','A','S','O','N','D'])
    plt.title('Tropospheric ozone due to STT to over the Southern Ocean')
    print("Created image at image/SO_extrapolation.ps")
    plt.savefig('images/SO_extrapolation.ps')
    
def event_profiles(station=2, data=None):
    '''
    Plot all the modelled profiles alongside sonde profiles for event dates
    inputs: station=2, data=None
        station: 0,1,or 2 for davis,mac,melb
        data: if the model data is already read, can be passed in as argument
    '''
    ## First read the station data
    if data is None:
        data=fio.read_GC_station(station)
    stn_name=data['Station'].split(' ')[0]
    dates=data['Date']
    
    ## Read the list of event dates from csv, created by blah.pro
    sondes= fio.read_site(station)
    SO3s  = sondes.o3ppbv # [time, altitude]
    SAlts = sondes.gph/1000 # [time, altitude] GPH in km
    eventdates=sondes.edates # this info is read from a csv I create in IDL
    print("saving profiles for station %s"%stn_name)
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
    Plot timeseries for each station, also shows sonde measurements
    '''
    f3, f3axes = plt.subplots(3, 1, sharex=True, figsize=(14,10))
    f3axes[2].set_xlabel('Date')
    f3axes[1].set_ylabel('Ozone molecules/cm2')
    
    # read the three files into a list
    files=[ fio.read_GC_station(f) for f in range(3) ]
    o3sondes = [fio.read_site(s) for s in range(3) ]
    
    # for each station do this
    # gc, os, f3ax, i = files[0], o3sondes[0], f3axes[0], range(len(files))[0]
    for gc, os, f3ax, i in zip(files, o3sondes, f3axes, range(len(files))):
        ## grab variables
        
        # Ozone data array [time] in molecules/cm2
        data=gc['O3TropColumn']
        sdata=np.array(os.tropvc)
        
        # create string title and turn tau's into matplotlib date numbers
        station=gc['Station']
        dates = gc['Date']
        sdates= np.array(os.dates)
        
        # plot time series for each station
        f3ax.plot_date(dates, data[:], '-b', label='GEOS-Chem')
        f3ax.plot_date(sdates, sdata[:], 'k*', label='sondes')
        f3ax.set_title(station)
        if i == 0: 
            f3ax.legend()
            f3ax.set_xlim([datetime(2003,9,1),datetime(2014,3,1)])
        
        # add dobson units
        newax=f3ax.twinx()
        if i == 1: newax.set_ylabel('Dobson Units')
        newylim= [1/2.69e16 * oldlim for oldlim in f3ax.get_ylim()]
        newax.set_ylim(newylim)
    
    
    # set plot titles
    f3.suptitle('Tropospheric ozone column',fontsize=21)
    
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
    Show yearly cycle of tropospheric ozone for our three sites
    If hour is set, only consider values at that hour (0, 6, 12, or 18)
    '''
    # read site data
    sites = [ fio.read_GC_station(p) for p in range(3) ]
    
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

def monthly_GC_profiles(hour=None):
    '''
    Profile mean and std deviation for each month for each site
    If you only want to consider a particular hour then set the hour parameter
        hour can be one of [0, 6, 12, 18, None]
    '''
    
    # read site data
    sites = [ fio.read_GC_station(p) for p in range(3) ]
    
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
        outfile='images/eventprofiles/%s_GC_monthprofiles.png'%stn_name
        plt.tight_layout()
        plt.subplots_adjust(top=0.94)
        plt.savefig(outfile)
        print("Image saved to %s"%outfile)
        plt.close(f)

def monthly_sonde_profiles():
    '''
    Profile mean and std deviation for each month for each site
    '''
    
    # read site data
    o3sondes = [ fio.read_site(p) for p in range(3) ]
    
    #for each station do this
    for site in o3sondes:
        
        # Set up plot
        f, axes = plt.subplots(4,3, sharex=True, sharey=True, figsize=(16,16))
        axes[3,1].set_xlabel('Ozone (ppb)')
        xlim=[0,125]
        axes[3,1].set_xlim(xlim)
        axes[1,0].set_ylabel('Altitude (km)')
        ylim=[0,14]
        axes[1,0].set_ylim(ylim)
        
        # Grab Ozone
        O3 = np.array(site.o3ppbv)
        # metres to kilometres
        TP = np.array(site.tp) # TP is already in km
        Z  = np.array(site.gph) / 1000.0 
        Znew= np.linspace(0,14,100)
        # need to vertically bin the O3 profiles,
        # interpolated at 100 points up to 14km
        means=np.zeros([12,100])
        stds =np.zeros([12,100])
        TPm = np.zeros(12)
        TPs = np.zeros(12)
        counts=np.zeros(12)
        
        # bin data into 12 months
        allmonths=np.array([ d.month for d in site.dates ])
        for month in range(12):
            inds=np.where(allmonths == month+1)[0]
            n = len(inds)
            # each profile needs to be interpolated up to 14km
            profs=np.zeros([n,100])
            for i in range(n):
                profs[i,:] = np.interp(Znew, Z[inds[i],:], O3[inds[i],:])
            means[month,:]=np.nanmean(profs,axis=0)
            stds[month,:] =np.nanstd(profs,axis=0)
            TPm[month] = np.nanmean(TP[inds])
            TPs[month] = np.nanstd(TP[inds])
            counts[month]=n
        
        # plot the mean profiles and shade the area of 1 stdev
        titles=np.array([['Dec','Jan','Feb'],['Mar','Apr','May'],['Jun','Jul','Aug'],['Sep','Oct','Nov']])
        months=np.array([[11,0,1],[2,3,4],[5,6,7],[8,9,10]])
        # colour each season!
        colours=['red','magenta','blue','green']
        for i in range(4):
            for j in range(3):
                plt.sca(axes[i,j]) # set current axis
                mind=months[i,j] # which month are we plotting
                X=means[mind,:]
                Xl=X-stds[mind,:]
                Xr=X+stds[mind,:]
                plt.plot(X, Znew , linewidth=3, color='k')
                plt.fill_betweenx(Znew, Xl, Xr, alpha=0.3, color=colours[i])
                plt.title(titles[i,j])
                Y=TPm[mind]
                Ye=TPs[mind]
                plt.plot(xlim, [Y, Y], 'k--', linewidth=2 )
                plt.fill_between(xlim, [Y+Ye,Y+Ye], [Y-Ye,Y-Ye], alpha=0.2, color='k')
                # add count text to upper corner
                plt.text(.75*xlim[1], .5, "N=%d"%counts[mind])
        
        # set title, and layout, then save figure
        stn_name=site.name
        f.suptitle("Monthly averaged profiles from ozonesondes over "+stn_name)
        outfile='images/%s_sonde_monthprofiles.png'%stn_name
        plt.tight_layout()
        plt.subplots_adjust(top=0.94)
        plt.savefig(outfile)
        print("Image saved to %s"%outfile)
        plt.close(f)

if __name__ == "__main__":
    print ("Running")
    
    #[event_profiles(s) for s in [0,1,2]]
    #time_series()
    #anomaly_corellation()
    corellation()
    #yearly_cycle()
    #monthly_GC_profiles()
    #monthly_sonde_profiles()
    #SO_extrapolation()
    #[monthly_GC_profiles(hour=h) for h in [0,6,12,18] ]
    