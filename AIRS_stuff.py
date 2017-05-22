# Create images from AIRS of CO total column on Event days
#

# plot library
import matplotlib
# don't display stuff, we are just saving to file:
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# local libraries
import fio as fio

# datetime library
from datetime import datetime

import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
from glob import glob

def read_AIRS_day(date):
    '''
    Read an AIRS day
    '''
    airsfolder = "/hpc/data/chemistry/CAC/GEOS_Chem/Satellite/AIRS/CO/daily/"

    # find single file match, eg: "*AIRS.2004.01.01.*"
    pattern=date.strftime("*AIRS.%Y.%m.%d.*")
    filename = glob(airsfolder+pattern)
    if len(filename) == 0:
        print("MISSING DATE: "+pattern)
        return([-1],[-1],[-1])
    filename = filename[0]
    
    # open file
    fh = nc.Dataset(filename, mode='r')
    
    # pull out variables of interest
    lons=fh.variables['Longitude'][:]
    lats=fh.variables['Latitude'][:]
    # Ascending is during the day at 130pm local time    
    totco = fh.variables['TotCO_A'][:]
    
    # close file
    fh.close()
    return (lats,lons,totco)
    
def plot_AIRS_day(date):
    '''
    Create a plot from the AIRS dataset for one day
    '''
    lats,lons,totco = read_AIRS_day(date)
    if lats[0]==-1: return (-1)
    # plot stuff
    lon0=lons.mean()
    lat0=lats.mean()
    
    # width, height in meters, 
    #lon = -137.5, 172.5
    #lat = 15.5, -75.5
    m=Basemap(llcrnrlat=-80,  urcrnrlat=20,
              llcrnrlon=-140, urcrnrlon=175,
              resolution='l',projection='merc',
              lat_0=lat0, lon_0=lon0)
    
    # lat lon are 1D, basemap uses 2D mesh
    lon,lat = np.meshgrid(lons,lats)
    xi, yi = m(lon,lat)
    
    # draw the CO total column onto the map
    cs = m.pcolor(xi,yi,np.squeeze(totco)) # squeeze removes any 1 length dimensions
    
    # set up consistent colour map (the colour bar)
    cmap = plt.cm.jet # blue to red
    plt.set_cmap(cmap)
    plt.clim(1e18, 3.5e18) # bounds for cmap    
    
    #add coastlines and equator
    m.drawcoastlines()
    m.drawparallels([0], labels=[0,0,0,0])
    
    #add title, colorbar
    cb=m.colorbar(cs,"right",size="5%", pad="2%")
    cb.set_label('CO')
    
    #ax.set_title('Total Column Ascending '+str(date))
    plt.title('Total Column CO'+date.strftime("%Y%m%d"))
    return(m)

def plot_all_events_AIRS(test=False):
    '''
    loop through sites and plot all events
    '''

    # Get the site data for each site
    all_sonde_files=[fio.read_sonde(s) for s in range(3)]
    
    # saving to here:
    imagesfolder="images/AIRS/"
    for sonde in all_sonde_files:
        for i in sonde.einds:
            date=sonde.dates[i]
            dstr=date.strftime("%Y%m%d")
            outf="%s%s_%s.png"%(imagesfolder,sonde.name,dstr)
            
            # Set up figure window
            fig=plt.figure(figsize=(12,6))
        
            # Plotted in seperate function:    
            m=plot_AIRS_day(date)
            # if no airs data on this day then save an empty plot
            if (m == -1) :
                plt.savefig(outf+'.missing.png')
                plt.close(fig)
                continue
            # add site marker
            x,y = m(sonde.lon, sonde.lat)
            m.plot(x, y, 'mo', markersize=6 )
        
            #save
            plt.savefig(outf, bbox_inches='tight')
            print ("Saved "+outf)
            plt.close(fig)
            if test: 
                return ()

def check_high_CO(date, site ,radius, threshold=2e18):
    '''
    If there is a column with greather than threshold [molecules/cm2] of CO, 
        within radius of site then return True
    '''
    lat=site.lat
    lon=site.lon

if __name__ == "__main__":
    plot_all_events_AIRS(test=True)
