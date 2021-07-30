# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 15:57:20 2020

@author: kse18nru
"""

#%% Import modules and scripts
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import iris
import iris.coord_categorisation
import iris.plot as iplt
import model_foundry as foundry
import pandas as pd

plt.ioff()
#%% global definitions
Case = '0p5km'
variable = 'air_temperature'

#%% first load and extract data
# main data cube 
cube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/{variable}/{Case}_{variable}_24hrs_g_301.nc', variable)
data = cube.data
data[np.where(data == 0)] = np.nan
# land mask data cube
lmcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/land_binary_mask/{Case}_land_binary_mask_301.nc')
lmdata = lmcube.data
lmlon = lmcube.coord('grid_longitude').points
lmlat = lmcube.coord('grid_latitude').points
# mean-sealevel-pressure cube
mspcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/air_pressure_at_sea_level/{Case}_air_pressure_at_sea_level_24hrs_301.nc')
mspdata = mspcube.data
print(len(mspdata))
msplon = mspcube.coord('grid_longitude').points
msplat = mspcube.coord('grid_latitude').points
print(len(msplon), len(msplat))
gphcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/geopotential_height/{Case}_geopotential_height_24hrs_301.nc')
gphdata = gphcube.data
gphlon = gphcube.coord('grid_longitude').points
gphlat = gphcube.coord('grid_latitude').points
gphcube = None
model_lat = cube.coord('grid_latitude').points
model_lon = cube.coord('grid_longitude').points
polelat = cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
polelon = cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
## free up memory/ unpoint cubes
#thcube = None, mspcube = None
print(cube)
pressure_labels = ['200hPa','300hPa', '500hPa','800hPa','850hPa','950hPa']

#%% mask out data above orography/land
print(np.shape(lmdata))
print(np.shape(data))
for i in range(48):
    data[i][np.where(lmdata == 1)] = np.nan
#%%

domain = 'fulldomain'
if domain == 'fulldomain':
    lons, lats = foundry.modf.unrotate_coords(model_lon, model_lat, polelon, polelat)
    msp_levels = np.arange(600, 1200, 0.5)
    gph_levels = np.arange(0,12000,10)
    
    counter = 0
    for level in range(1):
        for time in range(48):
           
            if (time) % 8 == 0:
                mspdata_at_t = mspdata[time]
                mspdata_at_t[np.where(lmdata == 1)] = np.nan
                counter += 1; print('counter:',counter, 'time:', time, 'level:', level)
                fig0,ax0,norm0 = foundry.plot_hoz_xsec_mass(data, lmdata, model_lon, model_lat, lons, lats, lmlon, lmlat, time, level,
                                                                 data_label = r'Temperature, K', domain = 'fulldomain', 
                                                                 flight = 301, savefig = True, res = Case, unrotated = True,
                                                                 variable_name = '1.5m air temperature', time_norm = True,
                                                                 variable_in_file_name = '1p5m_air_temperature_mancscale_wol',
                                                                 figure_path_pdf = f'../../Figures/PDF/301/P_levels/{variable}/{Case}/{domain}/surface_1p5m',
                                                                 figure_path_png = f'../../Figures/PNG/301/P_levels/{variable}/{Case}/{domain}/surface_1p5m',
                                                                 contour_type = 'pcolor', colourmap = 'plasma',sampling_rate = 'halfhourly',
                                                                 msp_contour = True, msp_data = mspdata_at_t/100, msp_lon = msplon, msp_lat = msplat,
                                                                 surface_variable = True, msp_levels = msp_levels, pressure_levels = True,
                                                                 gph_contour = False, gph_data = gphdata, gph_lon = gphlon, gph_lat = gphlat, gph_levels = gph_levels,
                                                                 gph_style = 'solid',figsize = (15,18),
                                                                 cscale_override = True, cscale_min = 268, cscale_max = 278)

    
                plt.close()