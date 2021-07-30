"""
Created on Sun Nov 10 14:39:23 2019

@author: Wilhelm Hodder

Script to create various horizontal cross-sections from model data on model levels
"""
#%% module imports
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import iris
import iris.coord_categorisation
#import CubeCrossSectioner_UK as ccs
import iris.quickplot as qplt
import cartopy
import cartopy.feature as cfeat
import cartopy.crs as ccrs
import model_functions as modf
import itertools
import datetime

#%% global configurations 

import warnings
warnings.filterwarnings("ignore") ##BAD! If something doesn't work, re-enable warnings to check
font = {'size' : 18}
matplotlib.rc('font', **font)

figure_path_pdf = '../Figures/PDF/301/Model_levels/Vertical_velocity/'
figure_path_png = '../Figures/PNG/301/Model_levels/Vertical_velocity/'

# Turn interactive plotting off
plt.ioff()

#%% load land data
land_file_name = '20180312T0000Z_IcelandGreenlandSeas_1p5km_RA1M_pf000.nc'
land_path = f'../Model_Data/u-bk574/nc/{land_file_name}'

modf.file_info(land_path)

stash = 'land_binary_mask'
land_mask, land_latitude, land_longitude = modf.return_cube_components(land_path, stash, 
                                                        lat_name = 'grid_latitude',
                                                        lon_name = 'grid_longitude')

#%% load first variable to analyse: Theta on theta levels
file_name = '20180312T0000Z_IcelandGreenlandSeas_1p5km_RA1M_pi000.nc'
path1 = f'../Model_Data/u-bk574/nc/'
stash1 = 'Moisture flux'
modf.file_info(f'{path1}{file_name}')#, stash1)
modf.cube_info(f'{path1}{file_name}', stash1)
modf.level_info(f'{path1}{file_name}', stash1, 3)
model_height = modf.get_model_levels(f'{path1}{file_name}', stash1)
print(model_height)
arguments = {'Hybrid height' : True}
data, lat, lon = modf.return_concatenated_cube_components(path1, 'i', '1p5km',
                                                          stash1,
                                                          file_number =4,
                                                          lat_name = 'grid_latitude',
                                                          lon_name = 'grid_longitude')

#%% unortate coordinate grid
cube = iris.load_cube(f'{path1}{file_name}', stash1)

lons, lats = modf.unrotate_coords(lon, lat, 
                                  cube.coord('grid_longitude').coord_system.grid_north_pole_longitude, 
                                  cube.coord('grid_latitude').coord_system.grid_north_pole_latitude)

#%% full domain plot
counter = 0
for level in range(10):
    for time in range(48):
        counter += 1; print(counter)
        fig = modf.make_horizontal_model_contourf_linear(data, land_mask[0][0], lon, lat, lons, lats, 
                                                         data_label = r'Moisture flux, [$kgm^{-2}s^{-1}$]', domain = 'fulldomain', 
                                                         flight = 301, savefig = True, res = '1p5km',
                                                         variable_name = 'Moisture flux', level = level, time = time,
                                                         variable_in_file_name = 'moisture_flux',
                                                         figure_path_pdf = '../Figures/PDF/301/Model_levels/Moisture_flux',
                                                         figure_path_png = '../Figures/PNG/301/Model_levels/Moisture_flux')
        plt.close()
