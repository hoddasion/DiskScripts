# to read the atmospheric flight data (flight 301 and more) and process it
# data comes in netcdf format
# program created 29/10/2018

import numpy as np
from netCDF4 import Dataset
#import pandas as pd


# =============================================================================
# files to look at: 'core_masin_20180312_r002_flight301_50hz.nc'
#                   'core_masin_20180312_r002_flight301_1hz.nc'
# =============================================================================

def load_file(datafile):
    # quick function to open and load netCDF datafiles
    dataset = Dataset(datafile, 'r', format = "NETCDF3_CLASSIC")
    return dataset

def load_keys(dataset):
    # extract and return array of variable names (keys)
    datakeys = dataset.variables.keys ()
    return datakeys



def show_variable_info(dataset, datakeys):
    for j in range(len(datakeys)):
        # to print out informations about variables in a long list
        if j == 0:
            print (dataset[datakeys[j]], '\n')
        if j % 2 == 0: #even number
            pass
        else: #uneven number
            print (dataset[datakeys[j]], '\n')
    return
 
    
def show_fillvalues(dataset, datakeys):
    # print a list of the fill values for each variable for reference
    for k in range(len(datakeys)):
        try:
            dummy = dataset[datakeys[k]]._FillValue
            print (datakeys[k], '=', dummy)
        except AttributeError:
            pass
    return










