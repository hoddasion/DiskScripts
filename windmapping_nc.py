# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 16:16:00 2019

@author: Wilhelm Hodder

To compute cross-sections and resolved maps of wind fields from model .nc output data.
"""
#%% Imports

# personal modules
import datareading2 as dar # open and load data from .nc files
import auxiliary_functions as aux
import personal_classes as percy # variable and data handling

# anaconda and standard library
import xarray as xr
import sys # for exit()
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from time import strftime # for figure name formatting
# installed
import metpy.calc as mpcalc # to calculate theta and theta_e and such like
from metpy.units import units #  to specificy metpy.calc units
print('Cell: Imports successful\n')

#%%  open file

# open and extract data from file
try:
    dataset = dar.load_file(f'../Model_Data/u-bk574/nc/20180312T0000Z_IcelandGreenlandSeas_4p4km_RA1M_pj006.nc')
    
except:
    print ('\nCell: File failed to open.\n\n\n')
    sys.exit()
print ('\nCell: File successfully opened.\n\n')
datakeys = dar.load_keys(dataset)
print (datakeys)