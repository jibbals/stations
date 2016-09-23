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
from matplotlib.dates import DateFormatter, date2num
import numpy as np
from datetime import datetime

# Local module for reading sonde dataset
import fio

# GLOBALS:
#

# seasonal colormap
seasonal_cmap=matplotlib.colors.ListedColormap(['fuchsia','chocolate','cyan','darkgreen'])
# colours for model vs sondes
model_colour='red'
data_colour='k'
col={'GEOS':model_colour,'Sonde':data_colour}
###########################################################################
#####################    Functions                         ################
###########################################################################

def summary_plots():
    '''
    Summary of seasonality, STT altitude, STT depth for each site.
    '''
    # read sonde event data
    sondes=[fio.read_sonde(s) for s in range(3)]
    
    # set up plot limits, labels etc..
    # three seperate summary bar plots.
    titles=['Event '+s for s in ['seasonality','altitude','depth']]
    pltnames=[ 'summary_'+s+'.png' for s in ['season','altitude','depth'] ]
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
        prev=np.zeros(12)
        for k in krange:
            inds=list(set(np.where(etypes == k)[0]) - set(fireinds))
            mons=[edates[i].month-1 for i in inds]
            histo=np.histogram(mons,bins)[0]
            ph.append(plt.bar(left, histo, width, color=evcolours[k],bottom=prev,linewidth=linewidth))
            prev=histo+prev
        # add fires 
        firem=[edates[i].month-1 for i in fireinds]
        histo=np.histogram(firem,bins)[0]
        ph.append(plt.bar(left,histo,width,color=evcolours[3],bottom=prev,linewidth=linewidth))
        
        plt.xlim([-0.5, 11.5])
        if xlab: plt.xlabel('month')
        plt.xticks(X,['J','F','M','A','M','J','J','A','S','O','N','D'])
        if legend:
            plt.legend( [ p[0] for p in ph ], label, loc='upper center')
        if ylab: plt.ylabel('occurences')
    
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
        if xlab: plt.xlabel('Altitude (km)')
        if legend:
            plt.legend( [ p[0] for p in ph ], label, loc='upper right')
        if ylab: plt.ylabel('occurences')
    
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
        

def SO_extrapolation(north=-35,south=-75):
    '''
    Rough estimate of extrapolation
    '''
    # Read sondes data
    sonde_files=[fio.read_sonde(s) for s in range(2)]
    # monthly likelihood = occurrences/sondecount each month, for davis and macquarie
    n_s=len(sonde_files)
    likelihood=np.zeros([12,n_s])
    # monthly flux percentage
    fluxpct=np.zeros([12,n_s])
    
    # for davis and macquarie:
    for i, os in enumerate(sonde_files):
        # indices of all events and all sondes within this month
        allmonths=np.array([ d.month for d in os.dates ])
        alleventsmonths=np.array([ d.month for d in os.edates ])
        # for each month
        fluxpc=np.array(os.eflux)/np.array(os.etropvc)
        for month in range(12):
            inds=np.where(allmonths == month+1)[0]
            einds=np.where(alleventsmonths == month+1)[0]
            likelihood[month,i] = len(einds)/ float(len(inds))
            # Flux pct per month
            if len(einds) != 0:
                fluxpct[month,i] = np.mean(fluxpc[einds])
    
    # model SO tropospheric O3 Column:
    GCData=fio.read_GC_global()
    SOTropO3=GCData.southernOceanTVC(north=north,south=south)
    
    # plot estimated flux on left hand axis
    # plot likelihoods and flux pct on the right hand axis
    f=np.mean(fluxpct, axis=1)
    l=np.mean(likelihood, axis=1)
    flux = SOTropO3 * f * l
    return f,l,flux, SOTropO3

def plot_SO_extrapolation(north=-35,south=-75):
    '''
    plot estimate of Southern Oceanic STT flux
    '''
    X=range(12)
    f,l,flux, SOTropO3=SO_extrapolation(north=north,south=south)
    
    print("%4e molecules/cm2/yr STT ozone contribution to the southern high latitudes"%np.sum(flux))
    
    plt.clf()
    fig, axes=plt.subplots(nrows=2,ncols=1,sharex=True,figsize=(14,13))
    plt.sca(axes[0])
    plt.plot(X,flux, color='black', linewidth=3, label="STT Flux")
    plt.title("Ozone flux from STTs")
    plt.ylabel('Ozone flux [molec/cm$^2$/yr]')
    plt.ylim([1e14, 6e15])
    
    ax=axes[1]
    plt.sca(ax)
    l1=plt.plot(X, SOTropO3, 'k', linewidth=2, label="$\Omega_{SO_{O_3}}$")
    #plt.title("Tropospheric ozone VCs in SO (GEOS-Chem)")
    plt.xlim([-0.5, 11.5])
    plt.xlabel('Month')
    plt.xticks(X,['J','F','M','A','M','J','J','A','S','O','N','D'])
    plt.ylabel('Molecules / cm$^3$')
    
    # plot percentages
    newax=plt.twinx()
    # likelihood pct * pct contribution
    l2=newax.plot(X,f*100*3, color='teal', label='f*3')
    l3=newax.plot(X,l*100, color='magenta', label='l')
    
    # axes and legends
    #newax.legend(loc=1)
    newax.set_ylabel('percent')
    newax.set_ylim([0,35])
    lns=l1+l2+l3
    plt.legend(lns,[ln.get_label() for ln in lns], loc=0)
    ax.set_ylabel('Ozone flux [molec/cm2]')
    plt.xlim([-0.5, 11.5])
    ax.set_xlabel('Month')
    plt.xticks(X,['J','F','M','A','M','J','J','A','S','O','N','D'])
    plt.suptitle('Tropospheric ozone due to STT to over the Southern Ocean',fontsize=26)
    print("Created image at image/SO_extrapolation.png")
    plt.savefig('images/SO_extrapolation.png')

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
    Profile mean and std deviation for each month for each site
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
        # Interpolate to grid up to 14km
        Znewlen = 200
        Znew= np.linspace(0,14,Znewlen)
        # need to vertically bin the O3 profiles,
        # interpolated at 100 points up to 14km
        means=np.zeros([4,Znewlen])
        medians=np.zeros([4,Znewlen])
        pcta = np.zeros([4,Znewlen]) # ath percentile
        pctb = np.zeros([4,Znewlen]) # bth percentile
        stds =np.zeros([4,Znewlen])
        TPm = np.zeros(4)
        counts = np.zeros(4)
        s_means=np.zeros([4,Znewlen])
        s_medians=np.zeros([4,Znewlen])
        s_pcta = np.zeros([4,Znewlen]) # 5th percentile
        s_pctb = np.zeros([4,Znewlen]) # 95th percentile
        s_stds =np.zeros([4,Znewlen])
        s_TPm = np.zeros(4)
        s_counts=np.zeros(4)
        
        # bin data into 4 seasons
        allmonths=np.array([ d.month for d in site['Date'] ])
        s_allmonths=np.array( [ d.month for d in sonde.dates ])
        for season in range(4):
            # determine coincident indices for GEOS-Chem vs sondes at particular hour
            sondeAtHour=np.array([ datetime(d.year,d.month,d.day,hour) for d in sonde.dates ])
            hourmatch = np.in1d(site['Date'], sondeAtHour, assume_unique=True)
            
            # find coincident indices matching the current month
            inds = np.where( ((allmonths == months[season,0]+1) + 
                (allmonths == months[season,1]+1) +
                (allmonths == months[season,2]+1))  * hourmatch )[0]
            s_inds=np.where((s_allmonths == months[season,0]+1) + 
                (s_allmonths == months[season,1]+1) +
                (s_allmonths == months[season,2]+1) )[0]
            n, s_n = len(inds), len(s_inds)
            counts[season]=n
            s_counts[season]=s_n
            
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
            s_means[season,:]=np.nanmean(s_profs,axis=0)
            s_medians[season,:]=np.nanmedian(s_profs,axis=0)
            s_pcta[season,:] = np.nanpercentile(s_profs,pctl,axis=0)
            s_pctb[season,:] = np.nanpercentile(s_profs,100-pctl,axis=0)
            s_stds[season,:] =np.nanstd(s_profs,axis=0)
            s_TPm[season] = np.nanmean(s_TP[s_inds])
            
        stn_name=site['Station'].split(' ')[0]
        
        # plot the median profiles and shade the area of 5th-95th percentile
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
            if j == 2:
                twinax=plt.twinx()
                twinax.set_yticks([]) # turn off ticks
                twinax.set_ylabel(seasonstr[i],fontsize=26)
    
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
        
        # pick out midday hour TODO:
        #
        ind=ind[0]
        
        O3 = data['O3'][ind,:]
        Altitude=data['Altitudes'][ind,:] *1e-3 # m to km
        
        # plot the modelled event
        f=plt.figure(1,figsize=(6,10))
        plt.plot(O3,Altitude,color=model_colour, 
                 linewidth=3, label='Modelled Profile')
        plt.ylabel('Altitude(km)')
        plt.xlabel('O3 (ppbv)')
        plt.xlim([5,120])
        plt.ylim([0,14])
        
        ## event profile on same plot
        Sind  = sondes.get_index(date)
        plt.plot(SO3s[Sind,:], SAlts[Sind,:], color=data_colour,
                 linewidth=3, label='Sonde Profile')
        if legend:
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
    f3axes[1].set_ylabel('Ozone (molec/cm2)')
    
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
        f3ax.plot_date(dates, data, linestyle='None', marker='.', color='grey') 
        f3ax.plot_date(sdates, sdata, linestyle='None', marker='*',
                       color=data_colour, label='Ozonesonde')
        f3ax.plot_date(mdates, mdata, linestyle='None', marker='.',
                       color=model_colour, label='GEOS-Chem')
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
    from scipy import stats
    
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
        slope,intercept,r_value,p_value,std_err= stats.linregress(Yos,Ygc)
        colors=cmap(season)
        f3ax.scatter(Yos, Ygc, color=colors, cmap=cmap, label='Trop O3')
        f3ax.plot(Yos, intercept+slope*Yos, 'k-', label='Regression')
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
    from scipy import stats
    
    # read the three files into a list
    files=[ fio.read_GC_station(f) for f in range(3) ]
    o3sondes = [fio.read_sonde(s) for s in range(3) ]
    
    f3, f3axes = plt.subplots(3, 1, figsize=(12,16))
    f3axes[2].set_xlabel('Sonde tropospheric O3 (molecules/cm2)')
    f3axes[1].set_ylabel('GEOS-Chem tropospheric O3 (molecules/cm2)')
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
        slope,intercept,r_value,p_value,std_err= stats.linregress(Yos,Ygc)
        
        colors=cmap(season)
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
    f3.suptitle('Correlation',fontsize=21)
    
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
    #check_extrapolation()
    plot_SO_extrapolation()
    #check_GC_output()
    #[event_profiles(s) for s in [0,1,2]]
    time_series()
    seasonal_profiles(hour=0,degradesondes=False)
    #summary_plots()
    #monthly_profiles(hour=0,degradesondes=False)
    #anomaly_correlation()
    #correlation()
    #yearly_cycle()
    #monthly_GC_profiles()
    #monthly_sonde_profiles()
    
    #[monthly_GC_profiles(hour=h) for h in [0,6,12,18] ]
    
