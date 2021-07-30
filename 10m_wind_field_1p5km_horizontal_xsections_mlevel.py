# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 14:20:40 2020

@author: Wilhelm Hodder

maiing horizontal xsections of 10m wind field
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

#%% global settings#

import warnings
warnings.filterwarnings("ignore") ##BAD! If something doesn't work, re-enable warnings to check
font = {'size' : 18}
matplotlib.rc('font', **font)

variable = '10m_wind_field'
figure_path_pdf = f'../Figures/PDF/301/Model_levels/{variable}/'
figure_path_png = f'../Figures/PNG/301/Model_levels/{variable}/'
figure_path_mp4 = f'../Figures/MP4/301/Model_levels/{variable}/'
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

#%%

datafile_path = '../Model_Data/u-bk574/nc'
datafile_name = '20180312T0000Z_IcelandGreenlandSeas_1p5km_RA1M_pg000.nc'
modf.file_info(datafile_path + '/' + datafile_name)

#%%
stash = 'x wind component (with respect to grid)'
cube = iris.load_cube(f'{datafile_path}/{datafile_name}', stash)
print(cube)

#%% load both wind components
data_path = f'{datafile_path}/{datafile_name}'
stash1 = 'y wind component (with respect to grid)'
u, u_lat, u_lon = modf.return_concatenated_cube_components(data_path, 'g', '1p5km',
                                                          stash1,
                                                          file_number =4,
                                                          lat_name = 'grid_latitude',
                                                          lon_name = 'grid_longitude')
stash2 = 'x wind component (with respect to grid)'
v, v_lat, v_lon = modf.return_concatenated_cube_components(data_path, 'g', '1p5km',
                                                          stash2,
                                                          file_number =4,
                                                          lat_name = 'grid_latitude',
                                                          lon_name = 'grid_longitude')