# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 14:56:19 2020

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


#%% global definitions
Case = '0p5km'
variable = 'air_potential_temperature'

if Case == '0p5km':
    #%% first load and extract data
    thcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/{variable}/{Case}_{variable}_24hrs.nc', variable)
    thdata = thcube.data
    lmcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/land_binary_mask/{Case}_land_binary_mask_301.nc')
    lmdata = lmcube.data
    model_lat = thcube.coord('grid_latitude').points
    model_lon = thcube.coord('grid_longitude').points
    polelat = thcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = thcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    ## free up memory/ unpoint cubes
    thcube = None
    
    #%% unrotate pole
    lons, lats = foundry.modf.unrotate_coords(model_lon, model_lat, polelon, polelat)
    
    #%% make plots
    counter = 0
    for level in range(20):
        for time in range(48):
            if level == 19:
                if time == 30:
                    counter = level*time + 1 + time; print('counter:',counter, 'time:', time, 'level:', level)
                    fig0,ax0, my_norm = foundry.plot_hoz_xsec_mass(thdata, lmdata, model_lon, model_lat, lons, lats, model_lon, model_lat, time, level,
                                                                     data_label = r'w, [ms$^{-1}$]', domain = 'fulldomain', 
                                                                     flight = 301, savefig = False, res = '0p5km', unrotated = False,
                                                                     variable_name = 'Potential temperature', time_norm = False,
                                                                     variable_in_file_name = 'air_potential_temperature',
                                                                     figure_path_pdf = f'../../Figures/PDF/301/Model_levels/{variable}/Horizontal/{Case}',
                                                                     figure_path_png = f'../../Figures/PNG/301/Model_levels/{variable}/Horizontal/{Case}',
                                                                     contour_type = 'pcolor', colourmap = 'plasma')
        
                    #%% load 60s dataset from database
        
                    # using pandas
                    database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_60s_301.txt', delimiter = ' ')
                    #print(database['meantime'])
                    #print(database[['altgps', 'lon', 'lat']])
                    leg = 6
                    db_leg = np.array(database['legno'])
                    db_alt = np.array(database['altgps'])[np.where(db_leg == leg)]
                    db_theta = np.array(database['theta'])[np.where(db_leg == leg)]
                    db_lon = np.array(database['lon'])[np.where(db_leg == leg)]
                    db_lat = np.array(database['lat'])[np.where(db_leg == leg)]
                    db_time = np.array(database['meantime'])[np.where(db_leg == leg)]/3600 # mean time coordinate in hours
                    print(np.mean(db_alt), 'm')
                    print(db_theta)
                    #%% Rotate coordinates to model frame
                    # extract rotated north pole coordinates
                    NPoleLon = lmcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
                    NPoleLat = lmcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
                    # rotate database coordinates
                    db_rot_lon, db_rot_lat = iris.analysis.cartography.rotate_pole(db_lon, db_lat, NPoleLon, NPoleLat)
                    
                    
                    q = ax0.scatter(db_rot_lon + 360, db_rot_lat, db_theta, norm = my_norm, cmap = 'plasma')
                    
                    plt.show()
                    
 #%%           
if Case == '1p5km':
    #%% first load and extract data
    thcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/{variable}/{Case}_{variable}_24hrs.nc', variable)
    thdata = thcube.data
    lmcube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/land_binary_mask/{Case}_land_binary_mask_301.nc')
    lmdata = lmcube.data
    model_lat = thcube.coord('grid_latitude').points
    model_lon = thcube.coord('grid_longitude').points
    polelat = thcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = thcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    ## free up memory/ unpoint cubes
    thcube = None; lmcube = None
    
    #%% Subset to fit 0p5km domain
    ## boundaries (in points):
    South = 90; North = 320; West = 75; East = 215
    
    # subset 4d data
    sub_data_region1 = foundry.modf.subset_4d(thdata, South, North, West, East)
           
    # two-dimensional subsetting of land-mask
    sub_land_region1 = foundry.modf.subset_2d(lmdata, South, North, West, East )
    
    # subset coordinates and unrotate
    sub_lon1 = foundry.modf.subset_1d(model_lon, West, East)
    sub_lat1 = foundry.modf.subset_1d(model_lat, South, North)
    #%% unrotate pole
    lons, lats = foundry.modf.unrotate_coords(sub_lon1, sub_lat1, polelon, polelat)
    
    #%% make plots
    counter = 0
    for level in range(20):
        for time in range(48):
            counter = level*time + 1 + time; print('counter:',counter, 'time:', time, 'level:', level)
            fig0,ax0 = foundry.plot_hoz_xsec_mass(sub_data_region1, sub_land_region1, sub_lon1, sub_lat1, lons, lats, sub_lon1, sub_lat1, time, level,
                                                             data_label = r'w, [ms$^{-1}$]', domain = 'fulldomain', 
                                                             flight = 301, savefig = True, res = '1p5km', unrotated = False,
                                                             variable_name = 'Potential temperature', time_norm = False,
                                                             variable_in_file_name = 'air_potential_temperature',
                                                             figure_path_pdf = f'../../Figures/PDF/301/Model_levels/{variable}/Horizontal/{Case}',
                                                             figure_path_png = f'../../Figures/PNG/301/Model_levels/{variable}/Horizontal/{Case}',
                                                             contour_type = 'pcolor', colourmap = 'plasma')

            plt.close()