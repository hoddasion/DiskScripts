# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 17:15:25 2020

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
suite = 'u-bu807'
stream = 'verc'

#%% first load and extract data
# main data cube 
xcube = iris.load_cube(f'../../Model_Data/{suite}/nc/Control/windvector_pl/{Case}_x_wind_24hrs_{stream}_301.nc', 'x_wind')
xdata = xcube.data
ycube = iris.load_cube(f'../../Model_Data/{suite}/nc/Control/windvector_pl/{Case}_y_wind_24hrs_{stream}_301.nc', 'y_wind')
ydata = ycube.data
magcube = iris.load_cube(f'../../Model_Data/{suite}/nc/Control/windvector_pl/{Case}_mag_wind_24hrs_{stream}_301.nc', 'mag_wind')
magdata = magcube.data
magdata[np.where(magdata == 0)] = np.nan


# land mask data cube
lmcube = iris.load_cube(f'../../Model_Data/{suite}/nc/Control/land_binary_mask/{Case}_land_binary_mask_flt301.nc')
lmdata = lmcube.data
lmlon = lmcube.coord('grid_longitude').points
lmlat = lmcube.coord('grid_latitude').points
# mean-sealevel-pressure cube
mspcube = iris.load_cube(f'../../Model_Data/{suite}/nc/Control/air_pressure_at_sea_level/{Case}_air_pressure_at_sea_level_24hrs_verb_301.nc')
mspdata = mspcube.data
print(np.shape(xdata), np.shape(ydata))
msplon = mspcube.coord('grid_longitude').points
msplat = mspcube.coord('grid_latitude').points
print(len(msplon), len(msplat))
gphcube = iris.load_cube(f'../../Model_Data/{suite}/nc/Control/geopotential_height/{Case}_geopotential_height_24hrs_verd_301.nc')
gphdata = gphcube.data
gphlon = gphcube.coord('grid_longitude').points
gphlat = gphcube.coord('grid_latitude').points
gphcube = None
model_lat = xcube.coord('grid_latitude').points
model_lon = xcube.coord('grid_longitude').points
polelat = xcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
polelon = xcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
## free up memory/ unpoint cubes
#thcube = None, mspcube = None
print(xcube)
pressure_labels = ['100hPa','150hPa','200hPa','250hPa','300hPa','400hPa','500hPa','600hPa','650hPa','700hPa','750hPa','800hPa','850hPa','925hPa','950hPa','1000hPa']
print(pressure_labels[0])


#%%

print(xcube,'\n', ycube,'\n', magcube)
print(xcube.coord('pressure'))
#%%
for i in range(16):
    timemeans = []
    for j in range(24):
        try:
            timemean = np.nanmean(gphdata[j][i])
            timemeans.append(timemean)
        except:
            break
    mean_at_level = np.nanmean(np.array(timemeans))
    print(mean_at_level)   
#(gph_levelmeans)
#%%
domain = 'peninsulas'
if domain == 'peninsulas':
    ## boundaries (in points):
    if Case == '4p4km':
        South = 115; North = 200; West = 150; East = 185
    if Case == '1p5km':
        South = 110; North = 350; West = 90; East = 185
    if Case == '0p5km':
        South = 80; North = -1; West = 50; East = 350
    # subset 4d data
    sub_x = foundry.modf.subset_4d(xdata, South, North, West, East)
    sub_y = foundry.modf.subset_4d(ydata, South, North, West, East)
    sub_m = foundry.modf.subset_4d(magdata, South, North, West, East)
    sub_msp = foundry.modf.subset_2d(mspdata, South, North, West, East) 
    sub_gph = foundry.modf.subset_4d(gphdata, South, North, West, East)      
    # two-dimensional subsetting of land-mask
    lm = foundry.modf.subset_2d(lmdata, South, North, West, East )
    
    ## filter out sea points
    #orog[np.where(lm != 1)] = np.nan
    # subset coordinates and unrotate
    sub_lon = foundry.modf.subset_1d(model_lon, West, East)
    sub_lat = foundry.modf.subset_1d(model_lat, South, North)
    sub_lmlon = foundry.modf.subset_1d(lmlon, West, East)
    sub_lmlat = foundry.modf.subset_1d(lmlat, South, North)
    gphlon = foundry.modf.subset_1d(gphlon, West, East)
    gphlat = foundry.modf.subset_1d(gphlat, South, North)
    msplon = foundry.modf.subset_1d(msplon, West, East)
    msplat = foundry.modf.subset_1d(msplat, South, North)
    #%% unrotate pole
    lons, lats = foundry.modf.unrotate_coords(sub_lon, sub_lat, polelon, polelat)
    # make list for lons lats levels to be plotted
    lats_levels = np.arange(50,70,0.5); print(lats_levels)
    lons_levels = np.sort(-np.arange(10,30,0.5)); print(lons_levels)
    msp_levels = np.arange(600, 1200, 1)
    gph_levels = np.arange(0,50000,10)
    print('#############\n',np.shape(sub_x), np.shape(sub_y), '\n###############\n')
    counter = 0
    level_idx = np.arange(16); print(level_idx)
    for level in level_idx:
        
            for time in range(9):
               
                
                try:
                    #mspdata_at_t[np.where(lmdata == 1)] = np.nan
                    counter += 1; print('counter:',counter, 'time:', time, 'level:', level)
                    fig0,ax0,norm0 = foundry.plot_hoz_xsec_mass(sub_m, lm, sub_lon, sub_lat, lons, lats, sub_lmlon, sub_lmlat, time, level,
                                                                     data_label = r'Windspeed, $ms^{-1}$', domain = domain, 
                                                                     flight = 301, savefig = True, res = Case, unrotated = True,
                                                                     variable_name = 'Horizontal windspeed', time_norm = True,
                                                                     variable_in_file_name = 'windvector_pl',
                                                                     figure_path_pdf = f'../../Figures/PDF/301/u-bu807/P_levels/windvector_pl/{Case}',
                                                                     figure_path_png = f'../../Figures/PNG/301/u-bu807/P_levels/windvector_pl/{Case}',
                                                                     contour_type = 'contourf', colourmap = 'cool',sampling_rate = 'threehourly',
                                                                     msp_contour = False ,
                                                                     surface_variable = False, msp_levels = msp_levels, pressure_levels = True,
                                                                     gph_contour = True, gph_data = sub_gph, gph_lon = gphlon, gph_lat = gphlat, gph_levels = gph_levels,
                                                                     quiver = True, quiver_u = sub_x, quiver_v = sub_y,
                                                                     quiver_n = 1, gph_style = 'solid',figsize=(15,18),
                                                                     cscale_override = False, cscale_min = 0, cscale_max = 17.5, verstash = True,
                                                                     pressure_labels = pressure_labels)
        
        
                    plt.close()
                except:
                    break
        
