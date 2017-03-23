# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 16:29:05 2017

@author: jesse
"""

import matplotlib
# Stop python from displaying images, faster to save images this way
#matplotlib.use('Agg')

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

sondes=fio.sondes(site=0)
        