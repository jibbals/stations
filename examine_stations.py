# use python to examine dataset
# Created 20160616 by three gunmen in a salon

###########################################################################
#####################    Imports and Globals               ################
###########################################################################

#IMPORTS:
#

import matplotlib
# Stop python from displaying images, faster to save images this way
matplotlib.use('Agg')

# font size global, also legend markers should only have 1 points in them
matplotlib.rcParams.update({'font.size': 15,'legend.numpoints':1,'legend.scatterpoints':1})

# plotting, dates, maths
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.dates import DateFormatter, date2num
import numpy as np
from datetime import datetime
from scipy.constants import N_A # avegadro's number

# Local module for reading sonde dataset
import fio
from JesseRegression import RMA

# GLOBALS:
#
__DEBUG__=False

# seasonal colormap
seasonal_cmap=matplotlib.colors.ListedColormap(['fuchsia','chocolate','cyan','darkgreen'])
# colours for model vs sondes
model_colour='red'
data_colour='k'
col={'GEOS':model_colour,'Sonde':data_colour}
###########################################################################
#####################    Functions                         ################
###########################################################################
def brief_summary():
    ''' Print a brief summary '''
    
    sondes=[fio.read_sonde(s) for s in range(3)]
    # print a brief summary of sondes
    for s in sondes:
        print(s.name)
        print("%d releases, %d events, "%(len(s.dates),len(s.edates)))
        print("first date: %s, end date: %s"%(s.dates[0],s.dates[-1]))
        print("misc, front, cutoff, fire")
        arr=np.array(s.etype)
        fires=np.array(s.fireflagged)
        arr[fires]=-1
        print("%4d,%5d,%4d,%4d"%(np.sum(arr == 0),np.sum(arr==1),np.sum(arr==2),np.sum(arr==-1)))
        
        month_inds=[]
        print("Monthly sonde releases")
        for i in range(12):
            month_inds.append(str(len(s.month_indices(i+1))))
        print ','.join(month_inds)
        
        # check points below 10km:
        #
        hinds = s.gph < 1e4 # truw when height < 10000m
        under10=np.nansum(hinds,axis=1) # sum how many measurements below 10km
        print(np.mean(under10),np.median(under10))
        #hist=plt.hist(under10,bins=12)
        
def summary_plots():
    '''
    Summary of seasonality, STT altitude, STT depth for each site.
    '''
    # read sonde event data
    sondes=[fio.read_sonde(s) for s in range(3)]
    
    # set up plot limits, labels etc..
    # three seperate summary bar plots.
    titles=['Event '+s for s in ['seasonality','altitude','depth']]
    pltnames=[ 'images/summary_'+s+'.png' for s in ['season','altitude','depth'] ]
    evcolours=['cyan','darkblue','darkcyan','red'] # colours based on event type or burning
    
    linewidth=0.0
    krange=[1,2,0] # Order for barplot stacks
    label=['frontal','cutoff','misc','fire']
    
    # Three kinds of plot will require 3 functions
    def plot_season(sonde,xlab=False,ylab=False, legend=False):
        X = np.arange(12)    # the x locations for the barchart
        left=X-0.5 # left side for barplot bins
        bins=np.linspace(-0.5,11.5,13)
        width=0.9 # bar widths can be array
        
        edates=sonde.edates
        etypes=np.array(sonde.etype)
        fireinds=np.where(np.array(sonde.fireflagged))[0]
        
        # loop over the three types we want to plot.
        ph=[] # plot handles
        obs=np.zeros(12)
        for m in range(12):
            obs[m]=np.sum([ s.month == m+1 for s in sonde.dates ])
        prev=np.zeros(12)
        rprev=np.zeros(12)
        for k in krange:
            inds=list(set(np.where(etypes == k)[0]) - set(fireinds))
            mons=[edates[i].month-1 for i in inds]
            histo=np.histogram(mons,bins)[0]
            # relative monthly occurrence
            rhisto=histo/obs*100.0
            #ph.append(plt.bar(left, histo, width, color=evcolours[k],bottom=prev,linewidth=linewidth))
            ph.append(plt.bar(left, rhisto, width, color=evcolours[k],bottom=rprev,linewidth=linewidth))
            prev=histo+prev
            rprev=rhisto+rprev
        # add fires 
        firem=[edates[i].month-1 for i in fireinds]
        histo=np.histogram(firem,bins)[0]
        rhisto=histo/obs*100.0
        #ph.append(plt.bar(left,histo,width,color=evcolours[3],bottom=prev,linewidth=linewidth))
        ph.append(plt.bar(left,rhisto,width,color=evcolours[3],bottom=rprev,linewidth=linewidth))
        
        plt.xlim([-0.5, 11.5])
        if xlab: plt.xlabel('month')
        plt.xticks(X,['J','F','M','A','M','J','J','A','S','O','N','D'])
        if legend:
            plt.legend( [ p[0] for p in ph ], label, loc='upper center')
        if ylab: plt.ylabel('Event frequency [%]')
    
    def plot_altitude(sonde,depth=False,xlab=False,ylab=False,legend=False):
        xlims=[4,14]
        bins=np.linspace(4,14,21)
        left=bins[0:20]
        width=0.5
        prev=np.zeros(20)
        epeaks=sonde.epeak
        if depth:
            xlims=[0,9]
            bins=np.linspace(0,9,19)
            left=bins[0:18]
            prev=np.zeros(18)
            # due to difference in how I interpolate the sondes in IDL this doesn't work yet
            epeaks=np.array([sonde.tp[s] for s in sonde.einds])-np.array(sonde.epeak)
            # just use idl provided data
            epeaks=sonde.edepth
        
        etypes=np.array(sonde.etype)
        fireinds=np.where(np.array(sonde.fireflagged))[0]
        
        # loop over the three types we want to plot.
        ph=[] # plot handles
        for k in krange:
            inds=list(set(np.where(etypes == k)[0]) - set(fireinds))
            alts=[epeaks[i] for i in inds]
            histo=np.histogram(alts,bins)[0]
            ph.append(plt.bar(left, histo, width, color=evcolours[k],bottom=prev,linewidth=linewidth))
            prev=histo+prev
        # add fires 
        firem=[epeaks[i] for i in fireinds]
        histo=np.histogram(firem,bins)[0]
        ph.append(plt.bar(left,histo,width,color=evcolours[3],bottom=prev,linewidth=linewidth))
        
        plt.xlim(xlims)
        if xlab: plt.xlabel(['Altitude (km)','Distance from tropopause (km)'][depth])
        if legend:
            plt.legend( [ p[0] for p in ph ], label, loc='upper right')
        if ylab: plt.ylabel('Event count')
    
    def plot_depth(sonde,xlab=False,ylab=False,legend=False):
        plot_altitude(sonde,depth=True,xlab=xlab,ylab=ylab,legend=legend)
    
    plot_funcs={0:plot_season, 1:plot_altitude, 2:plot_depth}
    for i in range(3):
        sonde=sondes[i]
        f,axes=plt.subplots(3,1,sharex=True, figsize=(14,13))
        
        # plot the sonde summaries
        for j, sonde in enumerate(sondes):
            plt.sca(axes[j])
            plot_funcs[i](sonde,ylab=(j==1),xlab=(j==2),legend=(j==1))
            plt.title(sonde.name)
        
        f.suptitle(titles[i],fontsize=24)
        plt.savefig(pltnames[i])
        print("Saved: %s"%pltnames[i])
        plt.close()

def plot_andrew_STT():
    # read data
    snames, stts=fio.read_ANDREW()
    
    linewidth=0.0
    
    X = np.arange(12)    # the x locations for the barchart
    left=X-0.5 # left side for barplot bins
    width=0.9 # bar widths can be array
    
    plt.figure(figsize=[14,11])
    # loop over the three sites we want to plot. (4 if we want laverton)
    for i in range(3):
        plt.subplot(311+i)
        plt.bar(left, stts[i,:], width, color='blue',linewidth=linewidth)
        if i==1: plt.ylabel('Event frequency [%]')
        plt.xlim([-0.5, 11.5])
        plt.xticks(X,['J','F','M','A','M','J','J','A','S','O','N','D'],fontsize=18)
        plt.ylim([0, 18.0])
        plt.yticks(np.arange(0,18.1,3))
        plt.title(snames[i],fontsize=24)
    plt.xlabel('month',fontsize=20)
    plt.suptitle('STT Proxy using $Brunt-Vi\\"{a}s\\"{a}l\\"{a}$ frequency',fontsize=28)
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    pltname='images/AndrewProxySTT.png'
    plt.savefig(pltname)
    print('saved '+pltname)

def check_weird_tp(year=2006):
    ''' show profiles where tropopause reading is weird '''
    sondes = [ fio.read_sonde(site=j) for j in range(3) ]
    for si,sonde in enumerate(sondes):
        name=sonde.name        
        print(name)
        tp = np.array(sonde.tp)
        tplr = np.array(sonde.tplr)
        tpo3 = np.array(sonde.tpo3)
        diffs=[]
        lows=[]
        for zangl in [True, False]:
            sonde._set_tps(zangl=zangl)
            for i,t in enumerate(tp):
                date=sonde.dates[i]
                if date.year == year:
                    dstr=date.strftime('%Y%m%d')
                    pname='images/eventprofiles/temp/'+['','zangl_'][zangl]
                    extra=''
                    if t < 5:
                        lows.append("%s: tp,tplr,tpo3: %5.3e %5.3e %5.3e"%(dstr,tp[i],tplr[i],tpo3[i]))
                        extra='_low'
                    elif np.abs(tplr[i]-tpo3[i]) > 3:
                        diffs.append("%s: tp,tplr,tpo3: %5.3e %5.3e %5.3e"%(dstr,tp[i],tplr[i],tpo3[i]))
                        extra='_dif'
                    pname=pname+'%s%s%s.png'%(name,dstr,extra)
                    fig=sonde.plot_profile(date=date, ytop=18, rh=True, alltps=True)
                    plt.savefig(pname)
                    plt.close(fig)
        print ("tp < 5km altitude")
        for low in lows: print(low)
        print ("tp difference > 3km")
        for diff in diffs: print(diff)
    
def seasonal_tropopause(show_event_tropopauses=False, shading=False):
    ''' Plot seasonal tropopause heights for each station '''
    #sonde data
    sondes = [ fio.read_sonde(site=j) for j in range(3) ]
    
    # interpolate up to 20km - every 100 metres
    f=plt.figure(figsize=(12,10))
    plt.tick_params(axis='x', which='major', labelsize=18)
    plt.tick_params(axis='y', which='major', labelsize=18)
    X=np.arange(12)
    Xstr=['J','F','M','A','M','J','J','A','S','O','N','D']
    colours=['k','chocolate','magenta']
    
    tp_m=np.zeros([12,3]) # median tropopause
    tp_u=np.zeros([12,3]) # upper percentile
    tp_d=np.zeros([12,3]) # downward percentile
    tp_e=np.zeros([12,3]) # event median tropopause
    for si,sonde in enumerate(sondes):
        for month in range(12):
            minds   = sonde.month_indices(month+1)
            tps     = np.array(sonde.tp)[minds]
            
            if show_event_tropopauses:
                eminds  = np.array(list(set.intersection(set(minds),set(sonde.einds))))
                # if no events in this month
                if len(eminds)==0: 
                    tp_e[month,si]=np.NaN
                else:
                    tp_e[month,si]  = np.nanmedian(np.array(sonde.tp)[eminds])
            tp_m[month,si]      = np.nanmedian(tps)
            tp_u[month,si]      = np.nanpercentile(tps, 90)
            tp_d[month,si]      = np.nanpercentile(tps, 10)
        # plot the median and shade the 80% percentile range
        plt.plot(X,tp_m[:,si],color=colours[si],label=sonde.name,linewidth=3)
        if shading:
            plt.fill_between(X, tp_d[:,si], tp_u[:,si], color=colours[si],alpha=0.4)
        else:
            plt.plot(X,tp_d[:,si], color=colours[si], linewidth=2, linestyle='--')
            plt.plot(X,tp_u[:,si], color=colours[si], linewidth=2, linestyle='--')
        if show_event_tropopauses:
            # dotted line on event tpheights
            plt.plot(X, tp_e[:,si], '--',color=colours[si])
        print(sonde.name)
        print("Max TP: %5.2e"%np.nanmax(sonde.tp))
        print("Min TP: %5.2e"%np.nanmin(sonde.tp))
    plt.xticks(X,Xstr)
    plt.xlim([-0.5,11.5])
    lbound,ubound=np.floor(np.nanmin(tp_d)), np.ceil(np.nanmax(tp_u))
    print([lbound,ubound])
    ylimits=[lbound,ubound]
    plt.ylim(ylimits)
    plt.ylabel('Altitude [km]',fontsize=20)
    plt.xlabel('Month',fontsize=20)
    plt.legend(fontsize=22,loc=0,frameon=False)
    ax=plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    pname='images/tpheights.png'
    plt.suptitle('Multi-year tropopause altitude',fontsize=28)
    plt.savefig(pname)
    print("Saved "+pname)
    plt.close(f)
    
def seasonal_tropozone():
    ''' Seasonal tropospheric ozone contour plots '''
    
    #sonde data
    sondes = [ fio.read_sonde(site=j) for j in range(3) ]
    
    # interpolate up to 20km - every 100 metres
    levels=np.arange(.1,20.001,.1)
    f=plt.figure(figsize=(15,15))
    plt.tick_params(axis='x', which='major', labelsize=24)
    plt.tick_params(axis='y', which='major', labelsize=18)
    X=np.arange(12)
    Xstr=['J','F','M','A','M','J','J','A','S','O','N','D']
    for sind,sonde in enumerate(sondes):    
        ppbv_m  = np.zeros([12,len(levels)])
        tp_m    = np.zeros([12,2])
        ppbvs   = sonde.o3ppbv
        alts    = sonde.gph/1000.0 # m to km
        tplr    = np.array(sonde.tplr) # in km
        tpo3    = np.array(sonde.tpo3) # in km
        
        for i in range(12):
            month_inds=sonde.month_indices(i+1)
            ppbv_month=np.zeros([len(month_inds),len(levels)])
            for j,k in enumerate(month_inds):
                ppbv_month[j,:]=np.interp(levels,alts[k,:],ppbvs[k,:],left=np.NaN,right=np.NaN)
            ppbv_m[i,:] = np.nanmean(ppbv_month,axis=0)
            
            tp_m[i,0] = np.nanmedian(tplr[month_inds])
            tp_m[i,1] = np.nanmedian(tpo3[month_inds])
        # plot each site
        plt.subplot(311 + sind)
        cf=plt.contourf(X,levels,ppbv_m.T,
                        levels=np.logspace(1,2.5,50),
                        cmap=plt.cm.jet,norm = LogNorm())
        plt.plot(X,tp_m[:,0],color='red') # plot lapse rate tropopause 
        plt.plot(X,tp_m[:,1],color='k') # plot ozone tropopause
        plt.title(sonde.name,fontsize=24)
        # plot axes lables
        plt.xlim([-0.5, 11.5])
        plt.xticks(X, Xstr)
        plt.ylim([0,17])
        if sind==1:
            plt.ylabel('GPH [km]')
    
    plt.xlabel('Month',fontsize=24)
    # set colour bar to the right
    f.subplots_adjust(right=0.8)
    cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar=f.colorbar(cf, cax=cbar_ax)
    cbar.set_label(label="Ozone [ppbv]",size=20)
    cbar_ticks=np.logspace(1,2.5,6)
    cbar_ticklabels=['%5.1f'%tick for tick in cbar_ticks]
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels(cbar_ticklabels)
    plt.suptitle("Multi-year average ozone",fontsize=28)
    pname='images/seasonaltropozone.png'
    plt.savefig(pname)
    print('Saved '+pname)
    plt.close(f)

def STT_extrapolation(Region):
    '''
    Estimate STT flux in a particular region
    Region: [S ,W ,N ,E]
    Returns: ILMT, flux, flux_range, ILMT_stdev
        ILMT= numpy array(4,12),
        flux=product on axis 0,
        flux_range is 2x12 lower to upper stddev for flux
        ILMT_stdev is 4x12 std deviations for each month in the multi year dataset
        I, L, M are mean measured monthly (I)mpact, (L)ikelihood, and (M)easurement count
        T is GEOS-Chem tropospheric ozone column, Flux is these 4 factors multiplied
    '''
    # Read sondes data
    all_sonde_files=[fio.read_sonde(s) for s in range(3)]
    
    # we only use sondes within extrapolation region...
    def in_region(lat,lon):
        return ((lat<Region[2])*(lat>Region[0])*(lon>Region[1])*(lon<Region[3]))
    sonde_files = [sf for sf in all_sonde_files if in_region(sf.lat,sf.lon)]
    if __DEBUG__:
        print("STT_Extrapolation")
        print(Region)
        print("Contains sites:")
        print([s.name for s in sonde_files])
    
    # model SO tropospheric O3 Column:
    GCData=fio.read_GC_global()
    TropO3, TropO3_stdev = GCData.averagedTVC(Region)
    
    # How many sondes have we:
    n_s=len(sonde_files)
    
    # monthly (I)mpact (event O3 / trop O3):
    # Events per measurement or (L)ikelihood:
    # (M)easurements per month
    # (T)ropO3
    ILMT = np.zeros([4,12,n_s])
    
    # Std_deviation of each factor
    ILMT_stdev=np.zeros([4,12,n_s])
    counts = np.zeros([12,n_s])
    
    # seasonal calculations
    ILMT_s = np.zeros([4,4,n_s])
    ILMT_stdev_s=np.zeros([4,4,n_s])
    counts_s = np.zeros([4,n_s])
    seasons=[[-1,0,1],[2,3,4],[5,6,7],[8,9,10]]
    
    # for each site:
    for ii, os in enumerate(sonde_files):
        # indices of all events and all sondes within this month
        allyears=np.array([ d.year for d in os.dates ])
        n_years=len(set(allyears))
        allmonths=np.array([ d.month for d in os.dates ])
        alleventsmonths=np.array([ d.month for d in os.edates ])
        alleventsyears=np.array([ d.year for d in os.edates ])
        # How many of each month do the measurements span?
        n_months_check = np.array([len(set([date.year for date in os.dates if date.month == month])) for month in np.arange(1,12.1)])
        
        ILMTi=np.ndarray([4,12,n_years]) +np.NaN# factors for each month each year
        
        # for each year
        for yi,year in enumerate(set(allyears)):
            # for each month
            for mi in range(12):
                inds=(allmonths == mi+1) * (allyears == year)
                n_m=np.sum(inds) # number of measurements
                einds=(alleventsmonths == mi+1) * (alleventsyears == year)
                n_e=np.sum(einds) # number of event detections
                counts[mi,ii]=counts[mi,ii]+n_e
                if n_m==0: 
                    if __DEBUG__: 
                        print("%s has no measurements on %d-%d"%(os.name,year,mi+1))
                    continue # no measuremets this month & year
                ILMTi[2,mi,yi] = n_m # count of this months measurements
                ILMTi[0,mi,yi] = np.NaN # Impact
                ILMTi[1,mi,yi] = 0.0 # Likelihood
                if n_e != 0:
                    ILMTi[0,mi,yi] = np.mean(np.array(os.eflux)[einds]/np.array(os.etropvc)[einds]) # Impact
                    ILMTi[1,mi,yi] = n_e / float(n_m) # Likelihood of event per measurement
            ILMTi[3,:,yi]=TropO3
        # end of year loop
        
        # take multiyear monthly average
        ILMT[:,:,ii]=np.nanmean(ILMTi,axis=2) 
        # get the stdev
        ILMT_stdev[:,:,ii]=np.nanstd(ILMTi,axis=2)
        ILMT_stdev[3,:,ii]=TropO3_stdev
        # When just one measurement used, set stdev to 100%
        ILMT_stdev[ILMT_stdev==0]=ILMT[ILMT_stdev==0] 
        
        # seasonal averages:
        for kk in range(4):
            for ll in range(4):
                ILMT_s[ll,kk,ii] = np.nanmean(ILMTi[ll,seasons[kk],:])
                ILMT_stdev_s[ll,kk,ii] = np.nanstd(ILMTi[ll,seasons[kk],:])
            counts_s[kk,ii] = np.sum(counts[seasons[kk],ii])
        
        if __DEBUG__:
            print("%s spans this many of each month:"%(os.name))
            print(n_months_check)
            print("With this many measurements per month:")
            print(ILMT[2,:,ii])
            print("With these likelihoods of event detection:")
            print(ILMT[1,:,ii])
    # end of sites loop
    
    ## Now we take the average between the sites
    #
    ILMT=np.mean(ILMT, axis=2) # Average between sondes
    ILMT_s=np.mean(ILMT_s, axis=2)
    
    flux = np.nanprod(ILMT,axis=0)
    flux_s = np.nanprod(ILMT_s, axis=0) # molecules/month
    
    # uncertainty reduced by root of independent measurement count
    # assume each sonde site is independent
    ILMT_stdev=np.nanmean(ILMT_stdev,axis=2) / np.sqrt(n_s) 
    ILMT_stdev_s=np.nanmean(ILMT_stdev_s,axis=2) / np.sqrt(n_s) 
    
    # proportional uncertainty
    uncertainty = ILMT_stdev / ILMT # [4,12]
    uncertainty_s = ILMT_stdev_s / ILMT_s # [4, 4]
    
    flux_uncertainty = np.nansum(uncertainty,axis=0) # sum of factor proportional uncertainties
    flux_range = np.zeros([2,12])
    flux_range[0,:] = flux * (1-flux_uncertainty)
    flux_range[1,:] = flux * (1+flux_uncertainty)
    
    if __DEBUG__: 
        print("ILMT_stdev %s:"%str(ILMT_stdev.shape))
        print (ILMT_stdev)
        print("flux uncertainty:")
        print (flux_uncertainty)
        print("flux_range:")
        for jj in range(12):
            print("%d: %.2e"%(jj+1,flux_range[1,jj]-flux_range[0,jj]))
    
    return {"ILMT":ILMT, "flux":flux, "flux_range":flux_range, 
            "ILMT_stdev":ILMT_stdev, "ILMT_s":ILMT_s, 
            "ILMT_stdev_s":ILMT_stdev_s,"flux_s":flux_s}

def STT_extrapolation_bracketed(Region, event_lifetime=2.5):
    '''
    Estimate STT flux in a particular region
    Region: [S ,W ,N ,E]
    Using ILMT from STT_extrapolation function, replacing M with a range 
        from 4-30 to give us an STT flux bracket
    Returns: ILMT, flux, flux_range, ILMT_stdev
    '''
    extrap = STT_extrapolation(Region)
    ILMT=extrap["ILMT"]
    ILMT_s=extrap["ILMT_s"]
    flux=extrap["flux"]
    flux_s=extrap["flux_s"]
    flux_range=extrap["flux_range"]
    ILMT_stdev=extrap["ILMT_stdev"]
    ILMT_stdev_s=extrap["ILMT_stdev_s"]
    
    # Reset M parameter using assumed event lifetime
    M=30.0/float(event_lifetime)
    ILMT[2,:] = M
    ILMT_s[2,:] = M
    ILMT_stdev[2,:] = 0.1*M # assume 10% error
    ILMT_stdev_s[2,:] = 0.1*M 
    
    # recalculate flux with new parameter
    flux[...] = np.nanprod(ILMT,axis=0)
    flux_s[...] = np.nanprod(ILMT_s,axis=0)

    # Set up bracketed range based on event lifetime of 1 day to 0.25 months
    prod=np.nanprod(ILMT[[0,1,3],:],axis=0)
    flux_range[0,:] = 4.0 * prod
    flux_range[1,:] = 30.0 * prod
    
    return extrap

def plot_extrapolation(Region, pltname='images/STT_extrapolation.png', Bracketed=False):
    '''
    Plot estimate of STT flux in some particular region
    Region: [S ,W ,N ,E]
    '''
    
    extrap_fn=[STT_extrapolation,STT_extrapolation_bracketed][Bracketed]
    extrap = extrap_fn(Region)
    ILMT=extrap["ILMT"]
    ILMT_s=extrap["ILMT_s"]
    flux=extrap["flux"]
    flux_s=extrap["flux_s"]
    flux_range=extrap["flux_range"]
    ILMT_stdev=extrap["ILMT_stdev"]
    ILMT_stdev_s=extrap["ILMT_stdev_s"] 
    
    I=ILMT[0,:]
    L=ILMT[1,:]
    GCTropO3=ILMT[3,:]
    
    # conversion to Tg/yr:
    gca=fio.get_GC_area()
    so_area=gca.region_area(Region) # m2
    g_per_mol=48 # g/Mol
    molecs=np.sum(flux) # molec/cm2/yr
    
    # [molec/cm2/month] * [cm2/m2] * [m2] * [Mol/molec] * [g/Mol] * Tg/g
    Tg_per_month= flux*1e4 * so_area * 1/N_A * g_per_mol * 1e-12 # = Tg/month
    
    # set ylimits for plot
    ylim0, ylim1=0.925*np.min(flux_range), 1.075*np.max(flux_range)
    ylim0tg, ylim1tg = 0.925*np.min(Tg_per_month), 1.075*np.max(Tg_per_month)
    
    # Plot the flux and the factors which are used to calculate it
    # set up plot and axes
    plt.clf()
    fig, axes=plt.subplots(nrows=2,ncols=1,sharex=True,figsize=(14,13))
    plt.sca(axes[1])
    
    # plot the flux line:
    X=range(12)
    plt.plot(X,flux, color='black', linewidth=3, label="STT Flux")
    # plot errorbars
    #plt.errorbar(X, flux, yerr=flux_error, linestyle="None", marker="None", color="k")
    plt.fill_between(X, flux_range[0,:], flux_range[1,:], color='grey', alpha='0.5')
    
    # plot limits, title, and labels
    plt.title("Ozone flux from STTs")
    plt.ylabel('Ozone flux [molecules cm$^{-2}$ month$^{-1}$]')
    plt.ylim([ylim0, ylim1])
    # second plot axes for other dimension
    rax=plt.twinx()
    rax.set_ylim([ylim0tg,ylim1tg])
    rax.set_ylabel('[Tg month$^{-1}$]')    
    # chuck formula onto plot
    plt.text(0.525,0.895,r'$\Delta \Omega_{trop O_3} = \Omega_{trop O_3} * I * L * M$', fontsize=28, transform = rax.transAxes)
    
    # plot factors on seperate subplot
    ax=axes[0]
    plt.sca(ax)
    l1=plt.plot(X, GCTropO3, 'k', linewidth=2, label="$\Omega_{trop O_3}$")
    #plt.title("Tropospheric ozone VCs in sub region (GEOS-Chem)")
    plt.xlim([-0.5, 11.5])
    plt.xlabel('Month')
    plt.xticks(X,['J','F','M','A','M','J','J','A','S','O','N','D'])
    plt.ylabel('Ozone [molecules cm$^{-2}$]')
    
    # plot percentages
    newax=plt.twinx()
    # likelihood pct * pct contribution
    l2=newax.plot(X,I*4, color='teal', label='I*4')
    l3=newax.plot(X,L, color='magenta', label='L')
    # l4=newax.plot(X,M/20.0, color='brown', label='M/20')
    
    # axes and legends
    #newax.legend(loc=1)
    newax.set_ylabel('factors')
    newax.set_ylim([0,0.45])
    lns=l1+l2+l3#+l4
    plt.legend(lns,[ln.get_label() for ln in lns], loc=0)
    sreg="[%dN, %dE, %dN, %dE]"%(Region[0],Region[1],Region[2],Region[3])
    plt.suptitle('Tropospheric ozone due to STT over %s'%sreg,fontsize=26)
    
    # save image
    plt.savefig(pltname)
    print("Created image at "+pltname)    
    plt.close(fig)

    ## Print Seasons
    
    # print both values
    print("%5.3e molecules/cm2/yr STT ozone contribution "%molecs)
    print("(%5.3e Tg/yr)"%np.sum(Tg_per_month))
    print("This occurs over %5.3e km2"%(so_area/1e6))
    
    # calculate monthly flux in kg/km2/month (to compare with skerlak2014)
    # flux : [molecules O3/ cm2 / month]
    #         x  * cm2/km2 * g/mol * kg/g * mol/molecule]
    flux_kg = flux * 1e10 * g_per_mol* 1e-3 * (1/N_A) # kg / km2 / month
    fks = flux_s * 1e10 * g_per_mol* 1e-3 * (1/N_A) # kg / km2 / month
    fus=np.sum(ILMT_stdev_s/ILMT_s * 100.0,axis=0) # sum along the parameters
    print ("seasonal [kg/km2/month]: ")
    seasons2_s=(fks[0],fus[0],fks[1],fus[1],fks[2],fus[2],fks[3],fus[3])
    seasons2=(np.sum(flux_kg[[-1,0,1]])/3.0, np.sum(flux_kg[[2,3,4]])/3.0,
               np.sum(flux_kg[[5,6,7]])/3.0, np.sum(flux_kg[[8,9,10]])/3.0)
    print("Summer(DJF) = %5.3f, Autumn(MAM) = %5.3f, Winter(JJA) = %5.3f, Sprint(SON) = %5.3f"%seasons2)
    print("Summer(DJF) = %5.3f(%5.1f%%), Autumn(MAM) = %5.3f(%5.1f%%), Winter(JJA) = %5.3f(%5.1f%%), Sprint(SON) = %5.3f(%5.1f%%)"%seasons2_s)
    print("monthly uncertainty:")
    flux_uncertainty=ILMT_stdev/ILMT*100.0
    tot_unc=np.nansum(flux_uncertainty,axis=0)
    for jj in range(12):
        outset=(jj+1,flux_uncertainty[0,jj],flux_uncertainty[1,jj],
                flux_uncertainty[2,jj],flux_uncertainty[3,jj],tot_unc[jj])
        print("%d : I:%.2f%%, L:%.2f%%, M:%.2f%%, T:%.2f%%, TOTAL: %3.2f"%outset)
    
    

def SO_extrapolation(north=-35,south=-75):
    '''
    Rough estimate of flux over southern ocean
    A wrapper for STT_extrapolation
    '''
    region=[south,-179.9,north,179.9] #Region: [S ,W ,N ,E]
    return STT_extrapolation(region)
    

def plot_SO_extrapolation(north=-35,south=-75):
    '''
    plot estimate of Southern Oceanic STT flux
    '''
    region=[south,-179.9,north,179.9] #Region: [S ,W ,N ,E]
    pltname='images/SO_extrapolation.png'
    plot_extrapolation(region,pltname=pltname)

def check_extrapolation():
    '''
    Calculate the extrapolation, and it's sensitivity to changes in range changes
    '''
    # model SO tropospheric O3 Column:
    GCData=fio.read_GC_global()
    
    # extrapolation using our lat range of -50 to -75
    print('calculating default extrapolation over -35 to -75:')
    sotropo3=np.mean(GCData.southernOceanTVC(north=-35,south=-75))
    
    # extrapolation with southern lat of -70 and -80
    sotropo3s1=np.mean(GCData.southernOceanTVC(north=-35,south=-70)) # South+5
    sotropo3s2=np.mean(GCData.southernOceanTVC(north=-35,south=-80)) # South-5
    
    # extrapolation with northern lat of -55 and -45
    sotropo3n1=np.mean(GCData.southernOceanTVC(north=-40,south=-75)) # North-5
    sotropo3n2=np.mean(GCData.southernOceanTVC(north=-30,south=-75)) # North+5
    
    dl=[] # [ds1, ds2, dn1, dn2]
    for d in [sotropo3s1,sotropo3s2,sotropo3n1,sotropo3n2]:
        dl.append( 100 * (sotropo3 - d) / sotropo3 )
    # differences
    print('relative differences with changing lat bounds: ')
    print('    %6s %6s '%('+5','-5') )
    print(' N :%6.3f %6.3f '%(dl[3],dl[2]) )
    print(' S :%6.3f %6.3f '%(dl[0],dl[1]) )

def seasonal_profiles(hour=0, degradesondes=False, pctl=10):
    '''
    Profile mean and std deviation for each season for each site
    If you only want to consider a particular hour then set the hour parameter
        hour can be one of [0, 6, 12, 18]
    set degradesondes=True to degrade the sondes to matching GEOS-Chem resolution before finding the average..
    '''
    
    # read site data
    sites = [ fio.read_GC_station(p) for p in range(3) ]
    o3sondes = [ fio.read_sonde(p) for p in range(3) ]
    
    # Model - Obs BIAS, TODO: record and print biases
    
    # some plot setups stuff
    months=np.array([[11,0,1],[2,3,4],[5,6,7],[8,9,10]])
    seasonstr=['Summer (DJF)','Autumn (MAM)','Winter (JJA)','Spring (SON)']
    seasonalcolours=seasonal_cmap.colors
    
    # Set up plot axes
    f, axes = plt.subplots(4,3, sharex=True, sharey=True, figsize=(16,16))
    axes[3,1].set_xlabel('Ozone (ppb)', fontsize=20)
    axes[3,0].set_ylabel('Altitude (km)', fontsize=20)
    xlim=[0,125]
    axes[3,1].set_xlim(xlim)
    ylim=[0,14]
    axes[1,0].set_ylim(ylim)
    
    #for each station do this
    # site,sonde=sites[1],o3sondes[1]
    for site, sonde, j in zip(sites, o3sondes, range(3)):
        
        # degrade sonde if we are doing that
        if degradesondes:
            sonde.degrade_vertically(site)
        
        # Grab Ozone
        O3 = site['O3']
        s_O3 = np.array(sonde.o3ppbv)
        
        # metres to kilometres
        s_TP = np.array(sonde.tp) # sonde TP is already in km
        TP = site['TropopauseAltitude'] / 1000.0
        Z  = site['Altitudes']/1000.0 
        s_Z  = np.array(sonde.gph) / 1000.0 
        # Interpolate to grid up to 15km
        Znewlen = 200
        Znew= np.linspace(0,15,Znewlen)
        # need to vertically bin the O3 profiles,
        # interpolated at 100 points up to 14km
        means=np.zeros([4,Znewlen])
        medians=np.zeros([4,Znewlen])
        pcta = np.zeros([4,Znewlen]) # ath percentile
        pctb = np.zeros([4,Znewlen]) # bth percentile
        stds =np.zeros([4,Znewlen])
        TPm = np.zeros(4)
        TPstd = np.zeros(4)
        counts = np.zeros(4)
        s_means=np.zeros([4,Znewlen])
        s_medians=np.zeros([4,Znewlen])
        s_pcta = np.zeros([4,Znewlen]) # ath percentile
        s_pctb = np.zeros([4,Znewlen]) # bth percentile
        s_stds =np.zeros([4,Znewlen])
        s_TPm = np.zeros(4)
        s_TPstd = np.zeros(4)
        s_counts=np.zeros(4)
        
        # bin data into 4 seasons
        allmonths=np.array([ d.month for d in site['Date'] ])
        s_allmonths=np.array( [ d.month for d in sonde.dates ])
        for season in range(4):
            # determine coincident indices for GEOS-Chem vs sondes at particular hour
            sondeAtHour=np.array([ datetime(d.year,d.month,d.day,hour) for d in sonde.dates ])
            hourmatch = np.in1d(site['Date'], sondeAtHour, assume_unique=True)
            s_hourmatch = np.in1d(sondeAtHour, site['Date'], assume_unique=False)
            
            # find coincident indices matching the current season
            inds = np.where( ((allmonths == months[season,0]+1) + 
                (allmonths == months[season,1]+1) +
                (allmonths == months[season,2]+1))  * hourmatch )[0]
            s_inds=np.where( ((s_allmonths == months[season,0]+1) + 
                (s_allmonths == months[season,1]+1) +
                (s_allmonths == months[season,2]+1)) * s_hourmatch )[0]
            n, s_n = len(inds), len(s_inds)
            counts[season]=n
            s_counts[season]=s_n
            
            if n != s_n:
                print(np.sum(hourmatch),np.sum(s_hourmatch))
                print(site['Date'][hourmatch][0])
                print(site['Date'][hourmatch][1])
                print(sondeAtHour[s_hourmatch][0])
                print(sondeAtHour[s_hourmatch][1])
                #for i,d in enumerate(site['Date'][hourmatch]):
                #    print(i,d.strftime('%Y%m%d'))
                #for i,d in enumerate(sondeAtHour[s_hourmatch]):
                #    print(i,d.strftime('%Y%m%d'))
                assert False, 'HNNNNGGGGGG'
            
            # each profile needs to be interpolated up to 14km
            profs=np.zeros([n,Znewlen])
            s_profs=np.zeros([s_n,Znewlen])
            for i in range(n):
                profs[i,:] = np.interp(Znew, Z[inds[i],:], O3[inds[i],:],left=np.NaN,right=np.NaN)
            for i in range(s_n):
                s_profs[i,:] = np.interp(Znew, s_Z[s_inds[i],:], s_O3[s_inds[i],:],left=np.NaN,right=np.NaN)
            means[season,:]=np.nanmean(profs,axis=0)
            medians[season,:]=np.nanmedian(profs,axis=0)
            pcta[season,:] = np.nanpercentile(profs,pctl,axis=0)
            pctb[season,:] = np.nanpercentile(profs,100-pctl,axis=0)
            stds[season,:] =np.nanstd(profs,axis=0)
            TPm[season] = np.nanmean(TP[inds])
            TPstd[season] = np.nanstd(TP[inds])
            s_means[season,:]=np.nanmean(s_profs,axis=0)
            s_medians[season,:]=np.nanmedian(s_profs,axis=0)
            s_pcta[season,:] = np.nanpercentile(s_profs,pctl,axis=0)
            s_pctb[season,:] = np.nanpercentile(s_profs,100-pctl,axis=0)
            s_stds[season,:] =np.nanstd(s_profs,axis=0)
            s_TPm[season] = np.nanmean(s_TP[s_inds])
            s_TPstd[season] = np.nanstd(s_TP[s_inds])
            
        stn_name=site['Station'].split(' ')[0]
        
        # plot the median profiles and shade the area of xth-(100-x)th percentile
        for i in range(4):
            plt.sca(axes[i,j]) # set current axis
            X=medians[i,:]                
            Xl=pcta[i,:]
            Xr=pctb[i,:]
            s_X=s_medians[i,:]
            s_Xl=s_pcta[i,:]
            s_Xr=s_pctb[i,:]
            
            ## plot averaged profiles + std from the model
            #
            plt.plot(X, Znew , linewidth=3, color=col['GEOS'])
            #plt.fill_betweenx(Znew, Xl, Xr, alpha=0.5, color=seasonalcolours[i])
            plt.plot(Xl, Znew, linewidth=2, color=col['GEOS'], linestyle='--')
            plt.plot(Xr, Znew, linewidth=2, color=col['GEOS'], linestyle='--')
            
            ## plot averaged profiles + std from the sondes
            #
            plt.plot(s_X, Znew , linewidth=3, color=col['Sonde'])
            plt.fill_betweenx(Znew, s_Xl, s_Xr, alpha=0.5, color=seasonalcolours[i])
            
            # plot tropopauses
            Y=TPm[i]
            s_Y=s_TPm[i]
            plt.plot(xlim, [Y, Y], '--', linewidth=2, color=col['GEOS'])
            plt.plot(xlim, [s_Y, s_Y], '--', linewidth=2, color=col['Sonde'])
            
            # add count text to lower right corner
            plt.text(.65*xlim[1], 2, "%d simulated"%counts[i],color=col['GEOS'])
            plt.text(.65*xlim[1], .5, "%d sondes"%s_counts[i])
            if i == 0:
                plt.title(stn_name, fontsize=26)
                print("Tropopause heights:")
                print("Station(season), GEOS-Chem(std), Sonde(std)")
            if j == 2:
                twinax=plt.twinx()
                twinax.set_yticks([]) # turn off ticks
                twinax.set_ylabel(seasonstr[i],fontsize=26,color=seasonalcolours[i])
                    
            
            # Print mean bias between model and obs tropopause heights.
            # 
            print ("%s(%s), TP:%6.3f(%6.3f), s_TP:%6.3f(%6.3f)"%(stn_name,seasonstr[i], TPm[i],TPstd[i],s_TPm[i],s_TPstd[i]))
        
    
    ## set title, and layout, then save figure
    #
    f.suptitle("Seasonal profiles",fontsize=32)
    outfile='images/eventprofiles/seasonalprofiles.png'
    if hour is not None: outfile='images/eventprofiles/seasonalprofiles%02d.png'%hour
    if degradesondes:
        outfile='images/eventprofiles/seasonalprofilesdegraded.png'%stn_name
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig(outfile)
    print("Image saved to %s"%outfile)
    plt.close(f)

def monthly_profiles(hour=0, degradesondes=False):
    '''
    Profile mean and std deviation for each month for each site
    If you only want to consider a particular hour then set the hour parameter
        hour can be one of [0, 6, 12, 18, None]
    set degradesondes=True to degrade the sondes to matching GEOS-Chem resolution before finding the average..
    '''
    
    # read site data
    sites = [ fio.read_GC_station(p) for p in range(3) ]
    o3sondes = [ fio.read_sonde(p) for p in range(3) ]
    
    # some plot setups stuff
    titles=np.array([['Dec','Jan','Feb'],['Mar','Apr','May'],['Jun','Jul','Aug'],['Sep','Oct','Nov']])
    months=np.array([[11,0,1],[2,3,4],[5,6,7],[8,9,10]])
    seasonalcolours=seasonal_cmap.colors
    
    #for each station do this
    # site,sonde=sites[1],o3sondes[1]
    for site, sonde in zip(sites, o3sondes):
        
        # degrade sonde if we are doing that
        if degradesondes:
            sonde.degrade_vertically(site)
        
        # Set up plot axes
        f, axes = plt.subplots(4,3, sharex=True, sharey=True, figsize=(16,16))
        axes[3,1].set_xlabel('Ozone (ppb)')
        xlim=[0,125]
        axes[3,1].set_xlim(xlim)
        axes[1,0].set_ylabel('Altitude (km)')
        ylim=[0,14]
        axes[1,0].set_ylim(ylim)
        
        # Grab Ozone
        O3 = site['O3']
        s_O3 = np.array(sonde.o3ppbv)
        
        # metres to kilometres
        s_TP = np.array(sonde.tp) # sonde TP is already in km
        TP = site['TropopauseAltitude'] / 1000.0
        Z  = site['Altitudes']/1000.0 
        s_Z  = np.array(sonde.gph) / 1000.0 
        Znew= np.linspace(0,14,100)
        # need to vertically bin the O3 profiles,
        # interpolated at 100 points up to 14km
        means=np.zeros([12,100])
        medians=np.zeros([12,100])
        pct5 = np.zeros([12,100]) # 5th percentile
        pct95 = np.zeros([12,100]) # 95th percentile
        stds =np.zeros([12,100])
        TPm = np.zeros(12)
        s_means=np.zeros([12,100])
        s_medians=np.zeros([12,100])
        s_pct5 = np.zeros([12,100]) # 5th percentile
        s_pct95 = np.zeros([12,100]) # 95th percentile
        s_stds =np.zeros([12,100])
        s_TPm = np.zeros(12)
        s_counts=np.zeros(12)
        
        # bin data into 12 months
        allmonths=np.array([ d.month for d in site['Date'] ])
        s_allmonths=np.array( [ d.month for d in sonde.dates ])
        for month in range(12):
            # find indices matching the current month
            inds=np.where(allmonths == month+1)[0]
            if hour is not None:
                allhours =np.array([ d.hour for d in site['Date'] ])
                inds = np.where( (allmonths == month+1) * (allhours==hour) )[0]
            s_inds=np.where(s_allmonths == month+1)[0]
            n, s_n = len(inds), len(s_inds)
            s_counts[month]=s_n
            
            # each profile needs to be interpolated up to 14km
            profs=np.zeros([n,100])
            s_profs=np.zeros([s_n,100])
            for i in range(n):
                profs[i,:] = np.interp(Znew, Z[inds[i],:], O3[inds[i],:],left=np.NaN,right=np.NaN)
            for i in range(s_n):
                s_profs[i,:] = np.interp(Znew, s_Z[s_inds[i],:], s_O3[s_inds[i],:],left=np.NaN,right=np.NaN)
            means[month,:]=np.nanmean(profs,axis=0)
            medians[month,:]=np.nanmedian(profs,axis=0)
            stds[month,:] =np.nanstd(profs,axis=0)
            pct5[month,:] = np.nanpercentile(profs,5,axis=0)
            pct95[month,:] = np.nanpercentile(profs,95,axis=0)
            TPm[month] = np.nanmean(TP[inds])
            s_means[month,:]=np.nanmean(s_profs,axis=0)
            s_medians[month,:]=np.nanmedian(s_profs,axis=0)
            s_stds[month,:] =np.nanstd(s_profs,axis=0)
            s_pct5[month,:] = np.nanpercentile(s_profs,5,axis=0)
            s_pct95[month,:] = np.nanpercentile(s_profs,95,axis=0)
            s_TPm[month] = np.nanmean(s_TP[s_inds])
            
        # plot the median profiles and shade the area of 5th-95th percentiles
        for i in range(4):
            for j in range(3):
                plt.sca(axes[i,j]) # set current axis
                mind=months[i,j]
                X=medians[mind,:]                
                Xl=X-pct5[mind,:]
                Xr=X+pct95[mind,:]
                s_X=s_medians[mind,:]
                s_Xl=s_X-s_pct5[mind,:]
                s_Xr=s_X+s_pct95[mind,:]
                
                # plot averaged profiles + std
                plt.plot(X, Znew , linewidth=3, color=col['GEOS'])
                plt.fill_betweenx(Znew, Xl, Xr, alpha=0.5, color=seasonalcolours[i])
                plt.plot(s_X, Znew , linewidth=3, color=col['Sonde'])
                plt.fill_betweenx(Znew, s_Xl, s_Xr, alpha=0.5, color=seasonalcolours[i])
                # plot tropopause
                Y=TPm[mind]
                s_Y=s_TPm[mind]
                plt.plot(xlim, [Y, Y], '--', linewidth=2, color=col['GEOS'])
                plt.plot(xlim, [s_Y, s_Y], '--', linewidth=2, color=col['Sonde'])
                # plot title and text
                plt.title(titles[i,j])
                # add count text to upper corner
                plt.text(.72*xlim[1], .5, "%d sondes"%s_counts[mind])
        
        # set title, and layout, then save figure
        stn_name=site['Station'].split(' ')[0]
        if hour is not None: stn_name+='_H%02d'%hour
        f.suptitle("Monthly Median Profiles over "+stn_name)
        outfile='images/eventprofiles/%s_monthprofiles.png'%stn_name
        if degradesondes:
            outfile='images/eventprofiles/%s_monthprofilesdegraded.png'%stn_name
        plt.tight_layout()
        plt.subplots_adjust(top=0.94)
        plt.savefig(outfile)
        print("Image saved to %s"%outfile)
        plt.close(f)

def plot_a_profile(date,sondes,GC,legend=False):
    ## First read the station data
    
    stn_name=GC['Station'].split(' ')[0]
    title=stn_name + date.strftime(" %Y %m %d")
    
    SO3s  = sondes.o3ppbv # [time, altitude]
    SAlts = sondes.gph/1000 # [time, altitude] GPH in km
    
    # find matching model profile
    ymd = np.array([ d.strftime('%Y%m%d') for d in GC['Date']])
    gcind = np.where( ymd=='%4d%02d%02d'%(date.year,date.month,date.day) )[0]
    
    # Using first hour index, since this is 7AM, 11AM, 11AM for our sites
    #
    ind=gcind[0]
    
    O3 = GC['O3'][ind,:]
    Altitude=GC['Altitudes'][ind,:] *1e-3 # m to km
    
    plt.plot(O3,Altitude,color=model_colour, 
             linewidth=3, label='GEOS-Chem')
    
    # plot the model pressure levels
    ms=30; mew=1.5    #marker size and edge width
    plt.plot(O3,Altitude,'_',color=model_colour, ms=ms,mew=mew)
    #plt.ylabel('Altitude(km)') #plt.xlabel('O3 (ppbv)') #plt.xlim([5,120]) #plt.ylim([0,14])
    plt.title(title, fontsize=20)
    
    ## event profile on same plot
    Sind  = sondes.get_index(date)
    plt.plot(SO3s[Sind,:], SAlts[Sind,:], color=data_colour,
             linewidth=3, label='Sonde')
    # plot the Sonde pressure levels
    plt.plot(SO3s[Sind,:], SAlts[Sind,:],'_', color=data_colour, ms=ms,mew=mew*0.8)
    if legend:
        plt.legend(loc='upper left',frameon=False)
    
def event_profiles_best():
    ''' plot the  three best profiles side by side '''
    
    dates=[datetime(2006,9,22),datetime(2012,12,19),datetime(2010,2,3)]
    f, axes=plt.subplots(1,3,sharey=True,sharex=True,figsize=(16,13))
    for i,date in enumerate(dates):
        #print('plotting %s at %s'%(date.strftime('%Y%m%d'),['dav','mac','melb'][i]))
        plt.sca(axes[i])
        sondes=fio.read_sonde(i)
        gc=fio.read_GC_station(i)
        plot_a_profile(date,sondes,gc, legend=(i==1))
    
    axes[0].set_ylabel('Altitude(km)')
    plt.sca(axes[1])
    plt.xlabel('O3 (ppbv)')
    plt.xlim([5,120])
    plt.ylim([0,14])
    pltname='images/eventprofiles/event_profile_comparison.png'
    plt.savefig(pltname)
    print (pltname+' saved!')
    
def event_profiles(station=2, data=None, legend=False):
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
    sondes= fio.read_sonde(station)
    SO3s  = sondes.o3ppbv # [time, altitude]
    SAlts = sondes.gph/1000 # [time, altitude] GPH in km
    eventdates=sondes.edates # this info is read from a csv I create in IDL
    print("saving profiles for station %s"%stn_name)
    ## At each event date plot the profiles side by side.
    for date in eventdates:
        title=stn_name + date.strftime(" %Y %m %d")
        outfile='images/eventprofiles/%s/%s_%s.png'%(stn_name,stn_name, date.strftime("%Y%m%d"))
        
        # find matching model profile
        ymd = np.array([ d.strftime('%Y%m%d') for d in dates])
        ind=np.where( ymd=='%4d%02d%02d'%(date.year,date.month,date.day) )[0]
        
        if len(ind) == 0:
            outfile='images/eventprofiles/%s/missing_%s_%s.png'%(stn_name,stn_name, date.strftime("%Y%m%d"))
            outf=open(outfile,'w')
            print('Missing %s'%date.strftime('%Y%m%d'))
            outf.close()
            continue
        
        # Using first hour index, since this is 7AM, 11AM, 11AM for our sites
        #
        ind=ind[0]
        
        O3 = data['O3'][ind,:]
        Altitude=data['Altitudes'][ind,:] *1e-3 # m to km
        
        # plot the modelled event
        f=plt.figure(1,figsize=(6,10))
        plt.plot(O3,Altitude,color=model_colour, 
                 linewidth=3, label='GEOS-Chem')
        # plot the model pressure levels
        plt.plot(O3,Altitude,'_',color=model_colour, markersize=12)
        plt.ylabel('Altitude(km)')
        plt.xlabel('O3 (ppbv)')
        plt.xlim([5,120])
        plt.ylim([0,14])
        plt.title(title, fontsize=24)
        
        ## event profile on same plot
        Sind  = sondes.get_index(date)
        plt.plot(SO3s[Sind,:], SAlts[Sind,:], color=data_colour,
                 linewidth=3, label='Sonde')
        if legend:
            plt.legend(loc='upper left',frameon=False)
        plt.savefig(outfile)
        print("plotted %s"%date.strftime('%Y%m%d'))
        plt.close(f)
    

def time_series(outfile='images/StationSeries.png'):
    '''
    Plot timeseries for each station, also shows sonde measurements
    '''
    f3, f3axes = plt.subplots(3, 1, sharex=True, figsize=(14,10))
    f3axes[2].set_xlabel('Date')
    f3axes[1].set_ylabel('$\Omega_{O_3}$ [molec cm$^{-2}$]')
    
    # read the three files into a list
    files=[ fio.read_GC_station(f) for f in range(3) ]
    o3sondes = [fio.read_sonde(s) for s in range(3) ]
    
    # for each station do this
    # gc, os, f3ax, i = files[0], o3sondes[0], f3axes[0], range(len(files))[0]
    for gc, os, f3ax, i in zip(files, o3sondes, f3axes, range(len(files))):
        ## grab variables
        
        # Ozone data array [time] in molecules/cm2
        data=gc['O3TropColumn']
        sdata=np.array(os.tropvc)
        
        # create string title and turn tau's into matplotlib date numbers
        station =gc['Station']
        dates   = gc['Date']
        hours   = np.array([ d.hour for d in dates ])
        h0inds  = np.where(hours==0)[0] # use UTC+0 times to grab local morning estimates only
        dates   = date2num(dates[h0inds])
        data    = data[h0inds]
        sdates  = date2num(np.array(os.dates))
        # also need to grab the model points which are on the same days as sondes
        s0dates = date2num(np.array([ datetime(sdate.year,sdate.month,sdate.day,0) for sdate in os.dates ]))
        matching= np.in1d(dates, s0dates, assume_unique=True)
        mdates  = dates[matching]
        mdata   = data[matching]
        
        # plot time series for each station
        print(station)
        print(dates)
        #f3ax.plot_date(dates, data, '-b', label='GEOS-Chem')
        #f3ax.plot_date(sdates, sdata, 'k*', label='sondes')
        
        f3ax.plot_date(dates, data, linestyle='None', marker='.', color='pink', alpha=0.5) 
        f3ax.plot_date(mdates, mdata, linestyle='None', marker='.',
                       color=model_colour, label='GEOS-Chem')
        f3ax.plot_date(sdates, sdata, linestyle='None', marker='*',
                       color=data_colour, label='Ozonesonde')
        
        f3ax.set_title(station)
        if i == 0: 
            f3ax.legend(numpoints=1)
            f3ax.set_xlim(date2num([datetime(2003,9,1),datetime(2014,3,1)]))
        
        # add dobson units
        newax=f3ax.twinx()
        if i == 1: newax.set_ylabel('Dobson Units')
        newylim= [1/2.69e16 * oldlim for oldlim in f3ax.get_ylim()]
        newax.set_ylim(newylim)
    
    # set plot titles
    f3.suptitle('Tropospheric ozone column ($\Omega_{O_3}$)',fontsize=21)
    
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
    plt.ylabel('ozone [molec cm$^{-2}$]')
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



###########################################################################
#####################    Finished Functions                ################
###########################################################################

def degraded_profile_comparison(station=0,date=datetime(2007,1,1)):
    '''
    Profile comparison at a site and date with sonde degraded to geos chem resolution
    '''
    
    # read site data
    site = fio.read_GC_station(station)
    os = fio.read_sonde(station)
    
    # some plot setup stuff
    xlim=[0,125]
    ylim=[0,14]
    col={'GEOS':model_colour,'Sonde':data_colour,'SondeDeg':'fuchsia'} # colours for model vs sondes
            
    # Set up plot axes
    f = plt.figure(figsize=(10,14))
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.ylabel('Altitude(km)')
    plt.xlabel('O3 (ppbv)')
    
    # plot original and geos chem for specific date:
    # find Closest matching sonde date
    sind=os.get_index(date)
    olddate=date
    date=os.dates[sind]
    yyyymmdd=date.strftime("%Y%m%d")
    print(olddate-date)
    
    dates=site['Date']
    ymd = np.array([ d.strftime('%Y%m%d') for d in dates])
    ind=np.where( ymd=='%4d%02d%02d'%(date.year,date.month,date.day) )[0]
    
    if len(ind) == 0:
        print('Missing %s from model output'%date.strftime('%Y%m%d'))
    ind=ind[0]
    
    # get GEOS-Chem profile and plot it
    GC = site['O3'][ind,:]
    GCZ=site['Altitudes'][ind,:] *1e-3 # m to km
    TP = site['TropopauseAltitude'][ind] / 1000.0
    plt.plot(GC,GCZ,marker='x',color=col['GEOS'],linewidth=3,label='GEOS-Chem')
    plt.plot(xlim,[TP,TP],'--',color=col['GEOS'])
    
    # get sonde profile and plot it
    prof=os.o3ppbv[sind,:]
    osz= os.gph[sind,:]*1e-3 # m to km
    ostp= os.tp[sind] # sonde TP is in km
    plt.plot(prof,osz,color=col['Sonde'],linewidth=3,label='Sonde(orig)')
    plt.plot(xlim,[ostp,ostp],'--',color=col['Sonde'])
    
    # Degrade sonde to GEOS vertical resolution
    # first what are the GEOS edges
    GCedges=np.zeros(73)
    GCedges[1:]=np.cumsum(site['BXHEIGHT'][ind,:])*1e-3
    prof2 = np.zeros(72)
    #degrade prof to average within GC box
    for i in range(72):
        prof2[i]=np.mean(prof[ (GCedges[i] < osz) * (osz < GCedges[i+1]) ])
    
    plt.plot(prof2,GCZ,marker='x',linewidth=3,color=col['SondeDeg'], label='Sonde (degraded)')
    
    # set title, and layout, then save figure
    plt.legend()
    stn_name=site['Station'].split(' ')[0]
    f.suptitle("Profile over %s on %s"%(stn_name,yyyymmdd))
    outfile='images/eventprofiles/%s_%s_profiles.png'%(stn_name,yyyymmdd)
    plt.savefig(outfile)
    print("Image saved to %s"%outfile)
    plt.close(f)

def anomaly_correlation(outfile='images/correlations_anomalies.png'):
    '''
    plot correlation of monthly average anomalies between sondes and GC where both occur on the same day.
    '''
    
    # read the three files into a list
    files=[ fio.read_GC_station(f) for f in range(3) ]
    o3sondes = [fio.read_sonde(s) for s in range(3) ]
    
    f3, f3axes = plt.subplots(3, 1, figsize=(12,16))
    f3axes[2].set_xlabel('Sonde tropospheric O3 relative anomaly')
    f3axes[1].set_ylabel('GEOS-Chem tropospheric O3 relative anomaly')
    ssnmap= {1:0,2:0,3:1,4:1,5:1,6:2,7:2,8:2,9:3,10:3,11:3,12:0} # map month to season:
    
    # Use seasonal colourmap (I defined at start)
    cmap=seasonal_cmap
    
    # for each station do this
    # gc, os, f3ax, i = files[0], o3sondes[0], f3axes[0], range(len(files))[0]
    for gc, os, f3ax, i in zip(files, o3sondes, f3axes, range(len(files))):
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
        season=[]
        for si, sdate in enumerate(osh0dates):
            if sdate in dates:
                ind=np.where(dates == sdate)[0]
                if len(ind) == 0 or np.isnan(sdata[si]):
                    continue
                else: 
                    ind=ind[0]
                #Work out anomaly and append to list
                mind=sdate.month-1
                Yos.append((sdata[si]-smean[mind])/smean[mind])
                Ygc.append((data[ind]-mean[mind])/mean[mind])
                season.append(ssnmap[sdate.month])
        
        Yos,Ygc = np.array(Yos),np.array(Ygc)
        # plot correlation coefficient
        m,b0,r,ci1,ci2= RMA(Yos,Ygc)
        colors=cmap(season)
        f3ax.scatter(Yos, Ygc, color=colors, cmap=cmap, label='Trop O3')
        f3ax.plot(Yos, b0+m*Yos, 'k-', label='Regression')
        f3ax.plot([-1,1], [-1, 1], 'k--', label='1-1 line')
        f3ax.set_title(station)
        f3ax.set_ylim(ylim)
        f3ax.set_xlim(xlim)
        txty=ylim[0]+0.1*(ylim[1]-ylim[0])
        txtx=xlim[0]+0.81*(xlim[1]-xlim[0])
        txty2=ylim[0]+0.17*(ylim[1]-ylim[0])
        txty3=ylim[0]+.24*(ylim[1]-ylim[0])
        plt.text(txtx,txty,"N=%d"%len(Yos))
        plt.text(txtx,txty2,"r=%5.3f"%r)
        plt.text(txtx,txty3,"slope=%5.3f"%m)
        if i==1: plt.legend()
        print("%s: r^2=%5.3f"%(station,r**2))
    # set plot titles
    f3.suptitle('Relative anomaly from monthly mean',fontsize=21)
    
    # add colourbar space to the right
    f3.subplots_adjust(right=0.85)
    cbar_ax = f3.add_axes([0.9, 0.25, 0.04, 0.5])
    
    # add colourbar, force the stupid thing to be the same as the one used in plotting...
    sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(vmin=1,vmax=4))
    sm._A=[]
    cb=f3.colorbar(sm,cax=cbar_ax)
    
    # set the 12 ticks nicely and roughly centred
    cb.set_ticks(np.linspace(1,4,4) + (2.5-np.linspace(1,4,4))/4.)
    cb.ax.set_yticklabels(['Summer','Autumn','Winter','Spring'],rotation=90)
    #cb.set_label('month')
    
    # save then close plots
    #plt.show()
    f3.savefig(outfile)
    print("Image saved to %s"%outfile)
    plt.close(f3)

def correlation(outfile='images/correlations.png'):
    '''
    plot correlation between sondes and GC where both occur on the same day.
    '''
    
    # read the three files into a list
    files=[ fio.read_GC_station(f) for f in range(3) ]
    o3sondes = [fio.read_sonde(s) for s in range(3) ]
    
    f3, f3axes = plt.subplots(3, 1, figsize=(12,16))
    f3axes[2].set_xlabel('Sonde $\Omega_{O_3}$ [molec cm$^{-2}$]',fontsize=20)
    f3axes[1].set_ylabel('GEOS-Chem $\Omega_{O_3}$ [molec cm$^{-2}$]',fontsize=20)
    xlims=[[1e17,1e18], [1e17,1e18], [1e17,1.5e18]]
    ylims=[[1e17,1.5e18], [1e17,1.5e18], [1e17, 2e18]]
    ssnmap= {1:0,2:0,3:1,4:1,5:1,6:2,7:2,8:2,9:3,10:3,11:3,12:0} # map month to season:
    
    # set up colorbar
    #cmap=plt.get_cmap('gist_rainbow', 4)
    cmap=seasonal_cmap
    
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
        season=[]
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
                season.append(ssnmap[sdate.month])
        
        Yos,Ygc = np.array(Yos),np.array(Ygc)
        # plot correlation coefficient
        m,b0,r,ci1,ci2= RMA(Yos,Ygc)
        
        colors=cmap(season)
        f3ax.scatter(Yos, Ygc, color=colors, cmap=cmap, label='Trop O3')
        f3ax.plot(Yos, b0+m*Yos, 'k-', label='Regression')
        f3ax.plot(xlim, xlim, 'k--', label='1-1 line')
        f3ax.set_title(station,fontsize=22)
        f3ax.set_ylim(ylim)
        f3ax.set_xlim(xlim)
        txty=ylim[0]+0.1*(ylim[1]-ylim[0])
        txtx=xlim[0]+0.81*(xlim[1]-xlim[0])
        txty2=ylim[0]+0.17*(ylim[1]-ylim[0])
        txty3=ylim[0]+.24*(ylim[1]-ylim[0])
        plt.text(txtx,txty,"N=%d"%len(Yos))
        plt.text(txtx,txty2,"r=%5.3f"%r)
        plt.text(txtx,txty3,"slope=%5.3f"%m)
        if i==1: plt.legend()
        print("%s: r^2=%5.3f"%(station,r**2))
    # set plot titles
    f3.suptitle('Correlation',fontsize=26)
    
    # add colourbar space to the right
    f3.subplots_adjust(right=0.85)
    cbar_ax = f3.add_axes([0.9, 0.25, 0.04, 0.5])
    
    # add colourbar, force the stupid thing to be the same as the one used in plotting...
    sm = plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(vmin=1,vmax=4))
    sm._A=[]
    cb=f3.colorbar(sm,cax=cbar_ax)
    
    # set the 12 ticks nicely and roughly centred
    cb.set_ticks(np.linspace(1,4,4) + (2.5-np.linspace(1,4,4))/4.)
    cb.ax.set_yticklabels(['Summer','Autumn','Winter','Spring'],rotation=90)
    #cb.set_label('month')
    
    # save then close plots
    #plt.show()
    f3.savefig(outfile)
    print("Image saved to %s"%outfile)
    plt.close(f3)

def check_GC_output():
    
    # read in GC data
    GCData=fio.read_GC_global()
    
    # Check the surface ozone
    data=GCData.O3density[0,0] # take first month, surface layer
    
    plt.figure(figsize=(14,10))
    GCData.mapData(data,label="Molecules cm$^{-3}$")
    plt.title("GEOS-Chem simulated surface ozone (Jan 2004)")
    pltname="images/GEOS-Chem_surface_ozone_example.png"
    plt.savefig(pltname)
    print(pltname+ " saved")
    plt.close()
    
    # Now check the Tropospheric VC
    data=GCData.O3tropVC[0] # first time slice
    
    plt.figure(figsize=(14,10))
    GCData.mapData(data,label="Molecules cm$^{-2}$")
    plt.title("GEOS-Chem simulated tropospheric ozone Column (Jan 2004)")
    plt.xlabel("$\Sigma_{z(troposphere)}($ Ozone ppb x $10^{-9}$x boxheight x$ N_{Air})$",fontsize=20)
    pltname="images/GEOS-Chem_tropospheric_VC_example.png"
    plt.savefig(pltname)
    print(pltname+ " saved")
    plt.close()

###########################################################################
#####################    Run Section                       ################
###########################################################################

if __name__ == "__main__":
    print ("Running")
    #brief_summary()
    #summary_plots()
    #event_profiles_best()
    #plot_andrew_STT()
    #check_extrapolation()
    #plot_SO_extrapolation()
    # N W S E regions:
    #Region1=[-60, 140, -35, 165] # region for Melb and Macca
    #Region2=[-70, 60, -55, 90] # region for Davis
    #plot_extrapolation(Region1,pltname='images/STT_extrapolation_MelbMac.png')
    #plot_extrapolation(Region2,pltname='images/STT_extrapolation_Dav.png')
    Reg_Melb=[-48,135,-28,155]
    Reg_Mac=[-65,149,-45, 169]
    Reg_Dav=[-79,68,-59,88]
    plot_extrapolation(Reg_Melb,pltname='images/STT_extrapolation_Melb_B.png',Bracketed=True)
    plot_extrapolation(Reg_Mac,pltname='images/STT_extrapolation_Mac_B.png',Bracketed=True)
    plot_extrapolation(Reg_Dav,pltname='images/STT_extrapolation_Dav_B.png',Bracketed=True)
    #check_weird_tp(2006)# look at profile of low tp sondes
    #seasonal_tropopause(shading=False) # plot tpheights.png
    #seasonal_tropozone() # plot seasonaltropozone.png
    #check_GC_output()
    #[event_profiles(s,legend = (s==1)) for s in [0,1,2]]
    #time_series()
    #seasonal_profiles(hour=0,degradesondes=False)
    #monthly_profiles(hour=0,degradesondes=False)
    #anomaly_correlation()
    #correlation()
    #yearly_cycle()
    
