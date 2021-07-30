# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 14:39:04 2020

@author: Wilhelm Hodder

Script to produce horizontal cross-sections of upward_air_velocity on model levels for flight 301
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
import multiprocessing


#%% global definitions
Case = '0p5km'
variable = 'upward_air_velocity'
Domain = 'fulldomain'
if Domain == 'fulldomain':
    #%% first load and extract data
    wcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/{variable}/{Case}_{variable}_24hrs.nc', variable)
    wdata = wcube.data
    lmcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/land_binary_mask/{Case}_land_binary_mask_301.nc')
    lmdata = lmcube.data
    model_lat = wcube.coord('grid_latitude').points
    model_lon = wcube.coord('grid_longitude').points
    polelat = wcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = wcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    ## free up memory/ unpoint cubes
    wcube = None; lmcube = None
    
    #%% unrotate pole
    lons, lats = foundry.modf.unrotate_coords(model_lon, model_lat, polelon, polelat)
    
    #%% loop to plot
    counter = 0
    select_levels = [11,20,23,39,51,56]
    for level in range(70):
        if level in (np.array(select_levels) - 1):
            for time in range(48):
                if time % 8 == 0:
                    counter += 1; print('counter:',counter, 'time:', time, 'level:', level)
                    fig0,ax0, norm0 = foundry.plot_hoz_xsec_mass(wdata, lmdata, model_lon, model_lat, lons, lats, model_lon, model_lat, time, level,
                                                                     data_label = r'w, [ms$^{-1}$]', domain = Domain, 
                                                                     flight = 301, savefig = True, res = Case,
                                                                     variable_name = 'Vertical velocity', 
                                                                     variable_in_file_name = 'upward_air_velocity_mcscale_ne5po5',
                                                                     figure_path_pdf = f'../../Figures/PDF/301/Model_levels/{variable}/Horizontal/{Case}',
                                                                     figure_path_png = f'../../Figures/PNG/301/Model_levels/{variable}/Horizontal/{Case}',
                                                                     contour_type = 'pcolor', colourmap = 'seismic', sym_cmap = True, figsize = (15,18),
                                                                     cscale_override = True, cscale_min = -5, cscale_max = 5)
if Case == 0:
    #%% first load and extract data
    wcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/{variable}/{Case}_{variable}_24hrs.nc', variable)
    wdata = wcube.data
    lmcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/land_binary_mask/{Case}_land_binary_mask_301.nc')
    lmdata = lmcube.data
    model_lat = wcube.coord('grid_latitude').points
    model_lon = wcube.coord('grid_longitude').points
    polelat = wcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = wcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    ## free up memory/ unpoint cubes
    wcube = None; lmcube = None
    
    #%% Subset to fit 0p5km domain
    ## boundaries (in points):
    South = 90; North = 320; West = 75; East = 215
    
    # subset 4d data
    sub_data_region1 = foundry.modf.subset_4d(wdata, South, North, West, East)
           
    # two-dimensional subsetting of land-mask
    sub_land_region1 = foundry.modf.subset_2d(lmdata, South, North, West, East )
    
    # subset coordinates and unrotate
    sub_lon1 = foundry.modf.subset_1d(model_lon, West, East)
    sub_lat1 = foundry.modf.subset_1d(model_lat, South, North)
    #%% unrotate pole
    lons, lats = foundry.modf.unrotate_coords(sub_lon1, sub_lat1, polelon, polelat)
    
    
    
    #%% loop over plots
    
    for level in range(20):
        for time in range(48):
            counter = level*time + 1 + time; print('counter:',counter, 'time:', time, 'level:', level)
            fig0,ax0 = foundry.plot_hoz_xsec_mass(sub_data_region1, sub_land_region1, sub_lon1, sub_lat1, lons, lats, sub_lon1, sub_lat1, time, level,
                                                             data_label = r'w, [ms$^{-1}$]', domain = 'fulldomain', 
                                                             flight = 301, savefig = True, res = '1p5km',
                                                             variable_name = 'Vertical velocity', 
                                                             variable_in_file_name = 'upward_air_velocity',
                                                             figure_path_pdf = f'../../Figures/PDF/301/Model_levels/{variable}/Horizontal/{Case}',
                                                             figure_path_png = f'../../Figures/PNG/301/Model_levels/{variable}/Horizontal/{Case}',
                                                             contour_type = 'pcolor', colourmap = 'seismic', sym_cmap = True)

            print('Figure saved.')
            plt.close()

if Case == 0:
    #%% first load and extract data
    wcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/{variable}/{Case}_{variable}_24hrs.nc', variable)
    wdata = wcube.data
    lmcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/land_binary_mask/{Case}_land_binary_mask_301.nc')
    lmdata = lmcube.data
    model_lat = wcube.coord('grid_latitude').points
    model_lon = wcube.coord('grid_longitude').points
    polelat = wcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = wcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    ## free up memory/ unpoint cubes
    wcube = None; lmcube = None
    
    #%% Subset to fit 0p5km domain
    ## boundaries (in points):
    South = 110; North = 190; West = 140; East = 200
    
    # subset 4d data
    sub_data_region1 = foundry.modf.subset_4d(wdata, South, North, West, East)
           
    # two-dimensional subsetting of land-mask
    sub_land_region1 = foundry.modf.subset_2d(lmdata, South, North, West, East )
    
    # subset coordinates and unrotate
    sub_lon1 = foundry.modf.subset_1d(model_lon, West, East)
    sub_lat1 = foundry.modf.subset_1d(model_lat, South, North)
    #%% unrotate pole
    lons, lats = foundry.modf.unrotate_coords(sub_lon1, sub_lat1, polelon, polelat)
    
    
    
    #%% loop over plots
    
    for level in range(20):
        for time in range(48):
            counter = level*time + 1 + time; print('counter:',counter, 'time:', time, 'level:', level)
            fig0,ax0 = foundry.plot_hoz_xsec_mass(sub_data_region1, sub_land_region1, sub_lon1, sub_lat1, lons, lats, sub_lon1, sub_lat1, time, level,
                                                             data_label = r'w, [ms$^{-1}$]', domain = 'fulldomain', 
                                                             flight = 301, savefig = True, res = '4p4km',
                                                             variable_name = 'Vertical velocity', 
                                                             variable_in_file_name = 'upward_air_velocity',
                                                             figure_path_pdf = f'../../Figures/PDF/301/Model_levels/{variable}/Horizontal/{Case}',
                                                             figure_path_png = f'../../Figures/PNG/301/Model_levels/{variable}/Horizontal/{Case}',
                                                             contour_type = 'pcolor', colourmap = 'seismic', sym_cmap = True)

            print('Figure saved.')
            plt.close()  