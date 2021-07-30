#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 14:21:49 2020

@author: Wilhelm Hodder

Script for generating a series of detailed maps based on model land mask and orography
with observational data and detailed flight path overlay
"""

#%% Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import iris.plot as iplt
import iris.quickplot as qplt
import iris
import iris.coord_categorisation
import model_foundry as foundry
#%% set global variables
leg = 1
res = '0p5km'
#%% load 60s dataset from database

# using pandas
database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_sr_301.txt', delimiter = ' ')
#print(database['meantime'])
#print(database[['altgps', 'lon', 'lat']])
db_leg = np.array(database['legno'])
db_alt = np.array(database['altgps'])#[np.where(db_leg == leg)]
db_lon = np.array(database['lon'])#[np.where(db_leg == leg)]
db_lat = np.array(database['lat'])#[np.where(db_leg == leg)]
db_time = np.array(database['meantime'])/3600#[np.where(db_leg == leg)]/3600 # mean time coordinate in hours
#%% load land cubes

nc_path = '../../Model_Data/u-bk574/nc/Control/'

og_file_name =  f'surface_altitude/{res}_surface_altitude_301.nc'
lm_file_name = f'land_binary_mask/{res}_land_binary_mask_301.nc'
og_cube = iris.load_cube(f'{nc_path}{og_file_name}', 'surface_altitude') # load orography cube
lm_cube = iris.load_cube(f'{nc_path}{lm_file_name}', 'land_binary_mask')

#%% extract data from lm_cube
lm_data = lm_cube.data
lm_x = lm_cube.coord('grid_longitude').points
lm_y = lm_cube.coord('grid_latitude').points

#%% rotate db coordinates

# extract rotated north pole coordinates
NPoleLon = lm_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
NPoleLat = lm_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
db_x, db_y = iris.analysis.cartography.rotate_pole(db_lon, db_lat, NPoleLon, NPoleLat)
print(lm_x,lm_y)
#%% unrotate model coordinates for grid overlay

lons, lats = foundry.modf.unrotate_coords(lm_x, lm_y, NPoleLon, NPoleLat)

#%% generate textfile with start and end coords of legs
## put coordinate data into pandas dataframe
#create lists first
start_lats = []
start_lons = []
end_lats = []
end_lons = []
legs = []
start_time = []
end_time = []
for i in range(np.max(db_leg) ):
    start_lats.append(db_lat[np.where(db_leg == i + 1)][0])
    start_lons.append(db_lon[np.where(db_leg == i + 1)][0])
    end_lats.append(db_lat[np.where(db_leg == i + 1)][-1])
    end_lons.append(db_lon[np.where(db_leg == i + 1)][-1])
    legs.append(i+1)
    start_time.append(db_time[np.where(db_leg == i + 1)][0])
    end_time.append(db_time[np.where(db_leg == i + 1)][-1])
coors = {'leg' : legs,
         'start_lats' : start_lats,
         'start_lons' : start_lons,
         'end_lats' : end_lats,
         'end_lons' : end_lons,
         'start_time' : start_time,
         'end_time' : end_time}

dataframe = pd.DataFrame(coors, columns = ['leg','start_lats','start_lons','end_lats','end_lons', 'start_time', 'end_time'])
print(dataframe)
dataframe.to_csv('../../Obvs_Data/Databases/301_leg_start_end_positions_sr.csv', index = False, header = True)
#%% plot full domain map
#fig0, ax0 = plt.subplots(1,1, figsize = (10,10))
#q0 = ax0.contour(lm_x, lm_y, lm_data)
#qlons = plt.contour(lm_x, lm_y, lons)
#qlats = plt.contour(lm_x, lm_y, lats)
#q1 = plt.scatter(db_x+360, db_y)
#plt.show()