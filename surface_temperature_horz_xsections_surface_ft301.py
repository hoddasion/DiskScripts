# -*- coding: utf-8 -*-
"""
Created on Sun Mar 22 18:41:52 2020

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
Case = '4p4km'
variable = 'surface_temperature'

#%% first load and extract data
# main data cube 
cube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/{variable}/{Case}_{variable}_24hrs_301.nc', variable)
data = cube.data
# land mask data cube
lmcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/land_binary_mask/{Case}_land_binary_mask_301.nc')
lmdata = lmcube.data
lmlon = lmcube.coord('grid_longitude').points
lmlat = lmcube.coord('grid_latitude').points
# mean-sealevel-pressure cube
mspcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/surface_air_pressure/{Case}_surface_air_pressure_24hrs_301.nc')
mspdata = mspcube.data
print(len(mspdata))
msplon = mspcube.coord('grid_longitude').points
msplat = mspcube.coord('grid_latitude').points
print(len(msplon), len(msplat))
model_lat = cube.coord('grid_latitude').points
model_lon = cube.coord('grid_longitude').points
polelat = cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
polelon = cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
## free up memory/ unpoint cubes
#thcube = None, mspcube = None
print(cube)

#%%

domain = 'fulldomain'
if domain == 'fulldomain':
    lons, lats = foundry.modf.unrotate_coords(model_lon, model_lat, polelon, polelat)
    msp_levels = np.arange(600, 1200, 5)
    
    

    counter = 0
    for level in range(1):
        for time in range(48):
            if (time )%8==0:
                counter += 1; print('counter:',counter, 'time:', time, 'level:', level)
                mspdata_at_t = mspdata[time]
                # replace msp instances above land with nans/ or prdduce masked array
                mspdata_at_t[np.where(lmdata == 1)] = np.nan
                fig0,ax0,norm0 = foundry.plot_hoz_xsec_mass(data, lmdata, model_lon, model_lat, lons, lats, lmlon, lmlat, time, level,
                                                                 data_label = r'Temperature, K', domain = 'fulldomain', 
                                                                 flight = 301, savefig = True, res = Case, unrotated = True,
                                                                 variable_name = 'Surface temperature', time_norm = True,
                                                                 variable_in_file_name = 'surface_temperature',
                                                                 figure_path_pdf = f'../../Figures/PDF/301/Surface/{variable}/{Case}/{domain}',
                                                                 figure_path_png = f'../../Figures/PNG/301/Surface/{variable}/{Case}/{domain}',
                                                                 contour_type = 'pcolor', colourmap = 'plasma',
                                                                 msp_contour = True, msp_data = mspdata_at_t/100, msp_lon = msplon, msp_lat = msplat,
                                                                 surface_variable = True, msp_levels = msp_levels)
    
                plt.close()