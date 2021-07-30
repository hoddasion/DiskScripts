#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 16:13:14 2020

@author: Wilhelm Hodder
"""


#%% Import modules and scripts
import sys
import numpy
import matplotlib.pyplot as plt
import matplotlib
import iris
import iris.coord_categorisation
import iris.plot as iplt
import model_foundry as foundry
import numpy as np

#%% global definitions
Case = '1p5km'
variable = 'wind_components_ml'
eastwind = 'm01s15i002'
northwind = 'm01s15i003'
if Case == '0p5km':
    #%% first load and extract data
    ucube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/{variable}/{Case}_{eastwind}_24hrs_301.nc', eastwind)
    u = ucube.data
    vcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/{variable}/{Case}_{northwind}_24hrs_301.nc', northwind)
    v = vcube.data
    lmcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/land_binary_mask/{Case}_land_binary_mask_301.nc')
    lmdata = lmcube.data
    model_lat = vcube.coord('grid_latitude').points
    model_lon = vcube.coord('grid_longitude').points
    land_lat = lmcube.coord('grid_latitude').points
    land_lon = lmcube.coord('grid_longitude').points
    polelat = vcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = vcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    ## free up memory/ unpoint cubes
    ecube = None; ncube = None; lmcube = None
    magdata = np.sqrt(u**2 + v**2)
    
    #%% unrotate pole
    lons, lats = foundry.modf.unrotate_coords(model_lon, model_lat, polelon, polelat)
    counter = 0
    for level in range(1):
        for time in range(1):
            counter += 1; print('counter:',counter, 'time:', time, 'level:', level)
            print(counter, 'model-level:', level + 1, 'time-step:', time + 1)
            fig = foundry.plot_hoz_xsec_mass(magdata, lmdata, model_lon, model_lat,  lons, lats, land_lon, land_lat,
                                                      level = level, time = time, colourmap = 'viridis',
                                                      quiver = True, quiver_u = u, quiver_v = v, quiver_n = 20,
                                                      variable_in_file_name = 'horizontal_velocity',
                                                      figure_path_pdf = f'../../Figures/PDF/301/Model_levels/Horizontal_velocity/Horizontal/{Case}',
                                                      figure_path_png = f'../../Figures/PNG/301/Model_levels/Horizontal_velocity/Horizontal/{Case}',
                                                      data_label = r'ms$^{-1}$', domain = 'fulldomain', 
                                                      flight = 301, savefig = True, res = '0p5km',
                                                      variable_name = 'Horizontal velocity',
                                                      contour_type = 'contourf')
            plt.close()
            
    
    plt.close()
    
if Case == '1p5km':
    #%% first load and extract data
    ucube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/{variable}/{Case}_{eastwind}_24hrs_301.nc', eastwind)
    u = ucube.data
    vcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/{variable}/{Case}_{northwind}_24hrs_301.nc', northwind)
    v = vcube.data
    lmcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/land_binary_mask/{Case}_land_binary_mask_301.nc')
    lmdata = lmcube.data
    model_lat = vcube.coord('grid_latitude').points
    model_lon = vcube.coord('grid_longitude').points
    land_lat = lmcube.coord('grid_latitude').points
    land_lon = lmcube.coord('grid_longitude').points
    polelat = vcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = vcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    ## free up memory/ unpoint cubes
    ecube = None; ncube = None; lmcube = None
    magdata = np.sqrt(u**2 + v**2)
    
    
    #%% Subset to fit 0p5km domain
    ## boundaries (in points):
    South = 90; North = 320; West = 75; East = 215
    
    # subset 4d data
    sub_v = foundry.modf.subset_4d(v, South, North, West, East)
    sub_u = foundry.modf.subset_4d(u, South, North, West, East)
    sub_mag = foundry.modf.subset_4d(magdata, South, North, West, East)
           
    # two-dimensional subsetting of land-mask
    sub_land = foundry.modf.subset_2d(lmdata, South, North, West, East )
    
    # subset coordinates and unrotate
    sub_lon = foundry.modf.subset_1d(model_lon, West, East)
    sub_lat = foundry.modf.subset_1d(model_lat, South, North)
    sub_lmlon = foundry.modf.subset_1d(land_lon, West, East)
    sub_lmlat = foundry.modf.subset_1d(land_lat, South, North)
    #%% unrotate pole
    lons, lats = foundry.modf.unrotate_coords(sub_lon, sub_lat, polelon, polelat)
    
    #%% unrotate pole
    lons, lats = foundry.modf.unrotate_coords(model_lon, model_lat, polelon, polelat)
    counter = 0
    for level in range(20):
        for time in range(48):
            counter += 1; print('counter:',counter, 'time:', time, 'level:', level)
            print(counter, 'model-level:', level + 1, 'time-step:', time + 1)
            fig = foundry.plot_hoz_xsec_mass(magdata, lmdata, model_lon, model_lat,  lons, lats, land_lon, land_lat,
                                                      level = level, time = time, colourmap = 'viridis',
                                                      quiver = True, quiver_u = u, quiver_v = v, quiver_n = 11,
                                                      variable_in_file_name = 'horizontal_velocity',
                                                      figure_path_pdf = f'../../Figures/PDF/301/Model_levels/Horizontal_velocity/Horizontal/{Case}',
                                                      figure_path_png = f'../../Figures/PNG/301/Model_levels/Horizontal_velocity/Horizontal/{Case}',
                                                      data_label = r'ms$^{-1}$', domain = 'fulldomain', 
                                                      flight = 301, savefig = True, res = '1p5km',
                                                      variable_name = 'Horizontal velocity',
                                                      contour_type = 'contourf')
            plt.close()
            
    
<<<<<<< HEAD
=======
    
    
if Case == '4p4km':
    #%% first load and extract data
    ucube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/{variable}/{Case}_{eastwind}_24hrs_301.nc', eastwind)
    u = ucube.data
    vcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/{variable}/{Case}_{northwind}_24hrs_301.nc', northwind)
    v = vcube.data
    lmcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/land_binary_mask/{Case}_land_binary_mask_301.nc')
    lmdata = lmcube.data
    model_lat = vcube.coord('grid_latitude').points
    model_lon = vcube.coord('grid_longitude').points
    land_lat = lmcube.coord('grid_latitude').points
    land_lon = lmcube.coord('grid_longitude').points
    polelat = vcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = vcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    ## free up memory/ unpoint cubes
    ecube = None; ncube = None; lmcube = None
    magdata = np.sqrt(u**2 + v**2)
    
    
    #%% Subset to fit 0p5km domain
    ## boundaries (in points):
    South = 110; North = 190; West = 140; East = 200
    
    # subset 4d data
    sub_v = foundry.modf.subset_4d(v, South, North, West, East)
    sub_u = foundry.modf.subset_4d(u, South, North, West, East)
    sub_mag = foundry.modf.subset_4d(magdata, South, North, West, East)
           
    # two-dimensional subsetting of land-mask
    sub_land = foundry.modf.subset_2d(lmdata, South, North, West, East )
    
    # subset coordinates and unrotate
    sub_lon = foundry.modf.subset_1d(model_lon, West, East)
    sub_lat = foundry.modf.subset_1d(model_lat, South, North)
    sub_lmlon = foundry.modf.subset_1d(land_lon, West, East)
    sub_lmlat = foundry.modf.subset_1d(land_lat, South, North)
    #%% unrotate pole
    lons, lats = foundry.modf.unrotate_coords(sub_lon, sub_lat, polelon, polelat)
    counter = 0
    for level in range(20):
        for time in range(48):
            counter += 1; print('counter:',counter, 'time:', time, 'level:', level)
            print(counter, 'model-level:', level + 1, 'time-step:', time + 1)
            fig = foundry.plot_hoz_xsec_mass(sub_mag, sub_land, sub_lon, sub_lat,  lons, lats, sub_lmlon, sub_lmlat,
                                                      level = level, time = time, colourmap = 'viridis',
                                                      quiver = True, quiver_u = sub_u, quiver_v = sub_v, quiver_n = 5,
                                                      variable_in_file_name = 'horizontal_velocity',
                                                      figure_path_pdf = f'../../Figures/PDF/301/Model_levels/Horizontal_velocity/{Case}',
                                                      figure_path_png = f'../../Figures/PNG/301/Model_levels/Horizontal_velocity/{Case}',
                                                      data_label = r'ms$^{-1}$', domain = 'fulldomain', 
                                                      flight = 301, savefig = True, res = '4p4km',
                                                      variable_name = 'Horizontal velocity',
                                                      contour_type = 'pcolor')
            plt.close()
            
    
>>>>>>> 5fa36f46c6430e5b5c7db522a2b4c2763ecc96f0
    plt.close()
    