# -*- coding: utf-8 -*-
"""
Created on Fri Jan  3 14:47:31 2020

@author: Wilhelm Hodder

Script to produce vertical cross-sections for vertical velocity
"""

#%% imports

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import iris.plot as iplt
import iris
import iris.coord_categorisation

import itertools
import datetime
import pandas as pd
from scipy.interpolate import interp2d

# homemade modules
import model_foundry as foundry

#%% global configurations 

import warnings
warnings.filterwarnings("ignore") ##BAD! If something doesn't work, re-enable warnings to check
font = {'size' : 18}
matplotlib.rc('font', **font)
variable = 'Vertical_velocity'
figure_path_pdf = f'../../Figures/PDF/301/Model_levels/{variable}/Vertical'
figure_path_png = f'../../Figures/PNG/301/Model_levels/{variable}/Vertical'
figure_path_mp4 = f'../../Figures/MP4/301/Model_levels/{variable}/Vertical'

res = '0p5km'




#%% load the cube anyway and extract altitude and surface/land data
nc_path = '../../Model_Data/u-bk574/nc/Control/'
w_file_name = f'upward_air_velocity/{res}_upward_air_velocity_24hrs.nc'
th_file_name = f'air_potential_temperature/{res}_air_potential_temperature_24hrs.nc'
og_file_name =  f'surface_altitude/{res}_surface_altitude_301.nc'
w_cube = iris.load_cube(f'{nc_path}{w_file_name}', 'upward_air_velocity') # load vertical velocity cube
th_cube = iris.load_cube(f'{nc_path}{th_file_name}', 'air_potential_temperature') # load theta cube
og_cube = iris.load_cube(f'{nc_path}{og_file_name}', 'surface_altitude') # load orography cube
#lm_cube = iris.load_cube(f'{nc_path}{og_file_name}', 'land_binary_mask')
print(w_cube)
#print(w_cube.coords('grid_longitude'),w_cube.coords('grid_latitude'))
#print(w_cube.coords('grid_latitude'))

#%%

xy_start = (-1.0,359.0)
xy_end = (0.5,360.5)

x_points = np.linspace(xy_start[0],xy_end[0],100)
y_points = np.linspace(xy_start[1],xy_end[1],100)
print(x_points,y_points)
sample_points = [('grid_latitude',x_points),('grid_longitude',y_points)]
time_index = 23
w_slice = w_cube[time_index].interpolate(sample_points, iris.analysis.Nearest())
print(w_slice)

#%%
fig = plt.figure( figsize = (18,15))
ax = iplt.contourf(w_slice, cmap = 'seismic')
plt.show()
