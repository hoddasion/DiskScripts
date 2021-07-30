# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 13:12:51 2019

@author: Wilhelm Hodder

Script for closer inspection of variables at model level and 
model validation with observational data
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
import database_functions as database

import warnings
warnings.filterwarnings("ignore") ##BAD! If something doesn't work, re-enable warnings to check
#%% load land data
land_file_name = '20180312T0000Z_IcelandGreenlandSeas_1p5km_RA1M_pf000.nc'
land_path = f'../Model_Data/u-bk574/nc/{land_file_name}'

modf.file_info(land_path)

stash = 'surface_altitude'
orography_set, land_latitude, land_longitude = modf.return_cube_components(land_path, stash, 
                                                        lat_name = 'grid_latitude',
                                                        lon_name = 'grid_longitude')
land_mask, land_latitude, land_longitude = modf.return_cube_components(land_path, 'land_binary_mask', 
                                                        lat_name = 'grid_latitude',
                                                        lon_name = 'grid_longitude')
orography = orography_set[0][0]
coastline = land_mask[0][0]

orocube = iris.load_cube(land_path, stash)

#%% load var from model output

file_name = '1p5km_upward_air_velocity_1.nc'
sub_folder = 'upward_air_velocity'
var_path = f'../Model_Data/u-bk574/nc/{sub_folder}/'
stash = 'upward_air_velocity'
modf.file_info(f'{var_path}{file_name}')#, stash1)
modf.cube_info(f'{var_path}{file_name}', stash)
modf.level_info(f'{var_path}{file_name}', stash, 3)
#model_height = modf.get_model_levels(f'{path1}{file_name}', stash2)
#print(model_height)
arguments = {'Hybrid height' : True}
mvar, mvar_lat, mvar_lon = modf.return_concatenated_cube_components(var_path, 'h', '1p5km',
                                                          stash,
                                                          file_number = 4,
                                                          lat_name = 'grid_latitude',
                                                          lon_name = 'grid_longitude')

#%%
# load the cube itself for good measure 

var_cube = iris.load_cube(f'{var_path}{file_name}', stash)
#print(var_cube.coord('altitude'))

#%% load obs from Database

oVar, oCoors, oAlts, oPrs, oTimes = database.load_database_as_pandasdataset(study_variable = 'w', legs = np.array([1]))

#%% interpolate cube onto flight altitude


mean_gps = np.nanmean(np.array(oAlts['altgps']))
print(mean_gps)
std_gps = np.nanstd(np.array(oAlts['altgps']))
print(std_gps)
print(var_cube.shape)
#x = len(var_cube.coords('grid_longitude')); y = len(var_cube.coords('grid_latitude'))
x = 420; y =400
print(np.ones((x,y))*mean_gps)
#%%
x = len(var_cube.coords('grid_longitude')); y = len(var_cube.coords('grid_latitude'))
sample_points =[('altitude', np.ones((x,y))*mean_gps)]
inter_cube = var_cube.interpolate(sample_points, iris.analysis.Linear())
