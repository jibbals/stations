# -*- coding: utf-8 -*-
"""
Created on Thu May 18 11:52:54 2017

@author: jesse
"""


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
import copy
import os
import errno

# Local module for reading sonde dataset
import fio
from JesseRegression import RMA

# GLOBALS:
#

# seasonal colormap
seasonal_cmap=matplotlib.colors.ListedColormap(['fuchsia','chocolate','cyan','darkgreen'])

# colours for model vs sondes
model_colour='red'
data_colour='k'
col={'GEOS':model_colour,'Sonde':data_colour}

def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def read_diff_events(sites=[0,1,2]):
    '''
    '''
    nam={"Davis":"dav","Macquarie":"mac","Melbourne":"melb"}
    all_sonde_files=[fio.read_sonde(s) for s in sites]
    slist=[]
    snames=[]
    for sonde in all_sonde_files:
        sn=sonde.name
        for ftype in ['_orig','_tpo3','_co985','_co98','_co97','_co96','_co95']:
            name=nam[sn]+ftype
            fname='../ozoneIDL/Data/'+name+".csv"
            sonde._set_events(csvname=fname)
            slist.append(copy.deepcopy(sonde))
            snames.append(name)
    return (slist, snames)

def summary():
    # read sonde data
    for sites in [[0],[1],[2]]:
        slist,snames=read_diff_events(sites=sites)
        ecount = [len(s.einds) for s in slist]
        mintp = [np.nanmin(s.tp) for s in slist]
        meantp = [np.nanmean(s.tp) for s in slist]
        maxtp = [np.nanmax(s.tp) for s in slist]
        
        head="%9s"%slist[0].name
        ecount = "events   "
        meantp = "mean tph "
        minmax = "tph bound"
        for sonde, sname in zip(slist,snames):
            
            head=head+'| %16s'%sname
            ecount=ecount+'| %16d'%len(sonde.einds)
            meantp=meantp+'| %16.2f'%np.nanmean(sonde.tp)
            minmax=minmax+'| %7.2f,%7.2f '%(np.nanmin(sonde.tp),np.nanmax(sonde.tp))
            
        print("")
        print(head)
        print(ecount)
        print(meantp)
        print(minmax)

def plot_all_events():
    '''
    Look at events after changes to the thresh and tp def
    '''
    # read sonde data
    slist,snames=read_diff_events()
    for sonde,sname in zip(slist,snames):
        sn=sonde.name
        for eind in sonde.einds:
            top = [16,13][sn == 'Davis']
            sonde.plot_profile(sonde.dates[eind],ytop=top)
            dstr=sonde.dates[eind].strftime('%Y%m%d')
            path='images/eventprofiles/%s/%s'%(sn,sname)
            fn=path+'/%s_%s.png'%(sname,dstr)
            make_sure_path_exists(path)
            print (fn)
            plt.savefig(fn)
            plt.close()
    
if __name__ == "__main__":
    summary()
    plot_all_events()