# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 15:21:35 2019

@author: Wilhelm Hodder
"""
#%% Import modules
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import iris
import iris.coord_categorisation
import iris.analysis.cartography as cairis
#import CubeCrossSectioner_UK as ccs
import iris.quickplot as qplt
import cartopy
import cartopy.feature as cfeat
import cartopy.crs as ccrs
import model_functions as modf # home-made functions
import itertools
import datetime # datetime formatting for figure names and titles
import xarray # for loading netcdf data
import pandas as pd # for loading observation data from text files in csv format


#%% functions and testing

def load_database_as_pandasdataset(path = '../Obvs_Data/Databases/', flightno = 301, 
                                   run_type = '60s', data_source = 'obs', file_extn = 'txt',
                                   datafile_name = None,
                                   delimiter = ' ', study_variable = None, 
                                   pressure = True, altitude = True, time = True, coord = True, legs = 'all'):
    """Function uses pandas read_csv to load selected variables from IGP database files, 
    which returns data structures as pandas DataFrame objects.
    Path to and name of data file and shape of delimiter are preset, but can be overridden.
    Study_variable can be stand-alone string or list/array-like of strings.
    Returns datastructures for
    Variables, Coordniates, Altitudes, Pressure, Times
    in that order."""
    ## load entire database first
    if datafile_name == None:
        datafile_name = f'IGP_flights_database_{data_source}_{run_type}_{flightno}.{file_extn}'
    database = pd.read_csv(f'{path}{datafile_name}', delimiter = delimiter)

    ## make list of variables to load out of database
    if altitude == True:
        altitude_names = ['altgps', 'altrad', 'altradSTD'] # names of altitude variables
    if pressure == True:
        pressure_name = 'prs' # name of pressure quantity; conditional
    if coord == True:
        coord_names = ['lon', 'lat', 'legno'] # names of coordinate variables and of leg designation
    if time == True:
        time_names = ['meantime', 'starttime', 'endtime'] # names of time quantitities
    
    legno = np.array(database['legno'])
    ## load variables out of database
    variable = database[study_variable] # DataFrame with the variable/variables of interest, eg. theta or w
    coordinates = database[coord_names] # DataFrame with longitude, latitude, and leg number
    altitudes = database[altitude_names] # DataFrame with gps and radar altitude, as well as radar standard deviation (m)
    pressure = database[pressure_name] # DataFrame with pressure in HPa
    times = database[time_names] # DataFrame with mean time, start time of run, and end time of run (seconds since midnight)
    
    ## replace time seconds with Datetime objects
    for i in range(len(times['meantime'])):
        # loop through each element and apply datetime.timedelta() 
        # this flags up a LOT of warnings, but it seems to work fine
        times['meantime'][i] = datetime.timedelta(0, int(times['meantime'][i]), 0)
        times['starttime'][i] = datetime.timedelta(0, int(times['starttime'][i]), 0)
        times['endtime'][i] = datetime.timedelta(0, int(times['endtime'][i]), 0)
    
    if legs == 'all':
        ## return all DataFrames
        return variable, coordinates, altitudes, pressure, times
    # initialise numpy arrays of variables with selected legs
    try:
        # if legs keyword has list of integers given, this will succeed
        # extract specified legs from data
        for i in range(len(legs)):
            if i == 0: # first step; need something to concatenate future iterations with
                conditions = legno == legs[0] # pre-generate list of booleans from condition
                # extract each leg in turn and also designate frst iteration of concatenated DataFrames
                conc_var = variable[conditions]
                conc_coor = coordinates[conditions]
                conc_alt = altitudes[conditions]
                conc_press = pressure[conditions]
                conc_time = times[conditions]
            else: # all steps after i = 0
                conditions = legno == legs[i] # pre-generate list of booleans from condition
                # extract each leg in turn
                temp_variable = variable[conditions]
                temp_coor = coordinates[conditions]
                temp_alt = altitudes[conditions]
                temp_press = pressure[conditions]
                temp_time = times[conditions]
                
                # concatenate DataFrames
                conc_var = pd.concat([conc_var, temp_variable])
                conc_coor = pd.concat([conc_coor, temp_coor])
                conc_alt = pd.concat([conc_alt, temp_alt])
                conc_press = pd.concat([conc_press, temp_press])
                conc_time = pd.concat([conc_time, temp_time])
        ## after loop is complete, return all concatenated DataFrames
        return conc_var, conc_coor, conc_alt, conc_press, conc_time
                
    except:
        print('Legs could not be extracted; falling back onto full flight.')
        ## return all DataFrames
        return variable, coordinates, altitudes, pressure, times


#Var, Coors, Alts, Prs, Times = load_database_as_pandasdataset(study_variable = 'w', legs = np.array([1]))

