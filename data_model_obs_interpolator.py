# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 13:57:10 2020

@author: Wilhelm Hodder

Script to interpolate UM model data onto new (linear) set of coordinates and 
output interpolated values into text file. Example use is to match model data points to the same
posiitons as observational data points
"""
#%% imports

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import iris
import iris.coord_categorisation
import model_functions as modf
import itertools
import datetime
import pandas as pd
import iris.analysis.trajectory as trajectory
# homemade modules
import model_foundry as foundry

#%% global configurations 

import warnings
warnings.filterwarnings("ignore") ##BAD! If something doesn't work, re-enable warnings to check
font = {'size' : 18}
matplotlib.rc('font', **font)

flight = 301
leg = 1
#%% load 60s dataset from database

# using pandas
database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_60s_301.txt', delimiter = ' ')
#print(database['meantime'])
#print(database[['altgps', 'lon', 'lat']])
db_leg = np.array(database['legno'])
db_alt = np.array(database['altgps'])[np.where(db_leg == leg)]
db_lon = np.array(database['lon'])[np.where(db_leg == leg)]
db_lat = np.array(database['lat'])[np.where(db_leg == leg)]
db_time = np.array(database['meantime'])[np.where(db_leg == leg)]/3600 # mean time coordinate in hours

#%% load model cube

variable = 'air_potential_temperature'
nc_path = '../../Model_Data/u-bk574/nc/Control/'
file_name = 'air_potential_temperature/0p5km_air_potential_temperature_24hrs.nc'
# load model data cube

cube = iris.load_cube(f'{nc_path}{file_name}', variable);
print(cube)

model_time = cube.coord('time').points
#print(model_time - model_time[0])
#%% Rotate coordinates to model frame
# extract rotated north pole coordinates
NPoleLon = cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
NPoleLat = cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
# rotate database coordinates
db_rot_lon, db_rot_lat = iris.analysis.cartography.rotate_pole(db_lon, db_lat, NPoleLon, NPoleLat)

#%% Temporal isolation and leg altitude meaning
if flight == 301:
    Year = 2018; Mon = 3; Day = 12
elif flight == 306:
    Year = 2018; Mon = 3; Day = 19
# seperate decimal time marker into hours and minuts, and cast as integers



# split decimal time value into hours and minutes, then cast elements as integers
Hours = (db_time//1).astype(int) 
Mins = ((db_time - Hours)*60//1).astype(int)
# initialise list to attach datetimes to
timeframe = [] # list for datetime-formatted times from database
for j in range(len(db_time)):
    # prepare list of values for datetime reformatting
    instance = [Year,Mon,Day,Hours[j],Mins[j]]
    # reformat list into datetime object
    instance = datetime.datetime(instance[0], instance[1], instance[2], instance[3], instance[4])
    timeframe.append(instance)
timeframe = np.array(timeframe)

#%% Nearest-neighbour interpolation

# define interpolation constraints
lon_constr = ('grid_longitude', db_rot_lon)
lat_constr = ('grid_latitude', db_rot_lat)
time_constr = ('time', np.mean(db_time + model_time[0]))
alt_constr = ('altitude', np.mean(db_alt))
sample_points = [time_constr]
#sample_points = [lon_constr, lat_constr, alt_constr, time_constr]
#cube_at_mean_time = cube.interpolate([time_constr], iris.analysis.Nearest())
coloumn = cube.interpolate([lon_constr, lat_constr, time_constr], iris.analysis.Nearest())
#model_points = coloumn.interpolate([alt_constr], iris.analysis.Linear())

print(coloumn)

fig = plt.figure()
q = iris.quickplot.scatter(coloumn.coord('grid_longitude'), coloumn.coord('grid_latitude'))

col_data = coloumn.data
col_alts = coloumn.coord('altitude').points
col_lat = coloumn.coord('grid_latitude').points
col_lon = coloumn.coord('grid_longitude').points 
temp_data = []
print(coloumn.data)
#print(col_data)
for i in range(len(col_data)):
    temp_single_coloumn = np.array([])
    
    for j in range(len(col_data[i])):
        temp_single_coloumn = np.concatenate((temp_single_coloumn, col_data[i][j]))
       
    temp_data.append(temp_single_coloumn)
col_data = np.array(temp_data)
print(col_data[15])    
