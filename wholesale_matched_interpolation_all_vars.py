# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 08:42:24 2021

@author: kse18nru
"""

#%% module imports
import iris
from iris.analysis.cartography import rotate_pole
import data_management as dmna
import pandas as pd
import numpy as np
import persinterpolate as pert
import matplotlib
import matplotlib.pyplot as plt
import xsection_workshop as workshop
import iris.coord_categorisation
import iris.plot as iplt
import matplotlib.gridspec as gridspec
from iris.analysis import trajectory

#%% 
res = '0p5km'
flight = '306'
suite = 'u-cf117'
config = 'LONGTAIL'
alt_idx = 40
xwind_name = 'x_wind'
ywind_name = 'y_wind'
stream_prefix = 'p'
tstart_idx = 24; tend_idx = 40
raw_all_points_sr = False
raw_all_points_60s = True
nc_path ='/'
#%%
if raw_all_points_sr:
    ## validate model with all obs points from short run averaged database (no distinction in legs applied)
    #%% load obs data
    obs_filename = f'IGP_flights_database_obs_sr_{flight}.txt'
    path = 'D:/Project/Obvs_Data/Databases/'
    
    database = pd.read_csv(f'{path}{obs_filename}', delimiter = ' ')
    db_coors = database[['lon','lat','altgps','starttime']]
    db_atmost = database[['wsp','wd','w','theta','q_BUCK']]
    db_fluxes = database[['sh','lh','tke','windstress']]
    
    
    #%% load and compute model data - atmostate
    theta_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_air_potential_temperature_24hrs_{stream_prefix}i_{flight}.nc', 'air_potential_temperature')[tstart_idx:tend_idx,:alt_idx]
    q_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_specific_humidity_24hrs_{stream_prefix}i_{flight}.nc', 'specific_humidity')[tstart_idx:tend_idx,:alt_idx]
    xwind_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_{xwind_name}_24hrs_{stream_prefix}h_{flight}.nc', xwind_name)[tstart_idx:tend_idx,:alt_idx]
    ywind_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_{ywind_name}_24hrs_{stream_prefix}h_{flight}.nc', ywind_name)[tstart_idx:tend_idx,:alt_idx]
    w_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_upward_air_velocity_24hrs_{stream_prefix}h_{flight}.nc', 'upward_air_velocity')[tstart_idx:tend_idx,:alt_idx]
    wsp_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_upward_air_velocity_24hrs_{stream_prefix}h_{flight}.nc', 'upward_air_velocity')[tstart_idx:tend_idx,:alt_idx]
    hoz_wind_data = (xwind_cube.data**2+ ywind_cube.data[:,:,:-1]**2)**0.5
    wsp_cube.data = hoz_wind_data
    #print(wsp_cube.coords('Time'))
    #%%
    wsp_cube.standard_name = None
    wsp_cube.long_name = 'windspeed'
    
    #%% load and compute model data - windstress
    
    ## windstress (two cumponents, need to be combined)
    # wse = windstrss eastward, wsn = windstress northward
    wse_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_atmosphere_downward_eastward_stress_24hrs_{stream_prefix}h_{flight}.nc', 'atmosphere_downward_eastward_stress')[tstart_idx:tend_idx,:alt_idx]
    wsn_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_atmosphere_downward_northward_stress_24hrs_{stream_prefix}h_{flight}.nc','atmosphere_downward_northward_stress')[tstart_idx:tend_idx,:alt_idx]
    # regrid wsn cube and combine cubes for magnitude
    #wsn_cube_reg = wsn_cube.regrid(wse_cube, iris.analysis.Linear())
    temp_cube =  iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_atmosphere_downward_eastward_stress_24hrs_{stream_prefix}h_{flight}.nc', 'atmosphere_downward_eastward_stress')[tstart_idx:tend_idx,:alt_idx]
    temp_cube.data = wsn_cube.data[:,:,:-1]
    
    ws_cube = (wse_cube**2 + temp_cube**2)**0.5
    ws_cube.standard_name = None
    ws_cube.long_name = 'windstress'
    #%%load and compute model data - heat fluxes
   
    # lh = latent heat, sh = sensible heat
    
    sh_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_upward_heat_flux_in_air_24hrs_{stream_prefix}i_{flight}.nc', 'upward_heat_flux_in_air')[tstart_idx:tend_idx,:alt_idx]
    lh_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_upward_latent_heat_flux_in_air_24hrs_{stream_prefix}i_{flight}.nc', 'upward_latent_heat_flux_in_air')[tstart_idx:tend_idx,:alt_idx]
    
    
    #%% rotate db coords
    
    polelat = w_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = w_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    rot_db_lon, rot_db_lat = rotate_pole(np.array(db_coors['lon']), np.array(db_coors['lat']), polelon, polelat)
    rot_db_lon = rot_db_lon + 360
    
    #%% adjust time coordinates
    
    ## db time is since 00:00 of same day, model time is since 1970
    ## convert db from econds in hours first
    db_time = np.array(db_coors['starttime'])/(60*60) + w_cube.coords('time')[0].points[0] - 0.5*(tstart_idx+1)
   
    time_coor = theta_cube.coords('time')[0]
    dates = time_coor.units.num2date(time_coor.points)
    print(dates)
    #%% interpolate coloumns
    db_alt = np.array(db_coors['altgps'])
    sample_points = [('grid_longitude', rot_db_lon), ('grid_latitude',rot_db_lat),('time', db_time) ]
    theta_cols = trajectory.interpolate(theta_cube, sample_points, method = 'nearest')
    q_cols = trajectory.interpolate(q_cube, sample_points, method = 'nearest')
    wsp_cols = trajectory.interpolate(wsp_cube, sample_points, method = 'nearest')
    w_cols = trajectory.interpolate(w_cube, sample_points, method = 'nearest')
    ws_cols = trajectory.interpolate(ws_cube, sample_points, method = 'nearest')
    lh_cols = trajectory.interpolate(lh_cube, sample_points, method = 'nearest')
    sh_cols = trajectory.interpolate(sh_cube, sample_points, method = 'nearest')
    print(theta_cols)
    #%%
    plt.plot(theta_cols.coords('grid_longitude')[0].points, theta_cols.coords('grid_latitude')[0].points)
    
    #%% interpolations
    theta_results = pert.coloumn_alt_interpolation(theta_cols, db_alt)
    q_results = pert.coloumn_alt_interpolation(q_cols, db_alt)
    wsp_results = pert.coloumn_alt_interpolation(wsp_cols, db_alt)
    w_results = pert.coloumn_alt_interpolation(w_cols, db_alt)
    ws_results = pert.coloumn_alt_interpolation(ws_cols, db_alt)
    lh_results = pert.coloumn_alt_interpolation(lh_cols, db_alt)
    sh_results = pert.coloumn_alt_interpolation(sh_cols, db_alt)
    
    #%% put data into pndas dataframe
    all_results = [theta_results,q_results,wsp_results,w_results,ws_results,lh_results,sh_results]
    print(np.shape(np.transpose(all_results)))
    #%%
    column_names = ['theta','q','wsp','w','ws','lh','sh']
    results_dataframe = pd.DataFrame(data = np.transpose(all_results),columns = column_names)
    print(results_dataframe)
    #%%
    results_dataframe.to_csv(f'D:/Project/Model_Data/{suite}/{config}_Matched_interpolated_values_{res}_sr_{flight}.txt')
    
    
#%%
if raw_all_points_60s:
    ## validate model with all obs points from 60s run averaged database (no distinction in legs applied)
    #%% load obs data
    obs_filename = f'IGP_flights_database_obs_60s_{flight}.txt'
    path = 'D:/Project/Obvs_Data/Databases/'
    
    database = pd.read_csv(f'{path}{obs_filename}', delimiter = ' ')
    db_coors = database[['lon','lat','altgps','starttime']]
    db_atmost = database[['wsp','wd','w','theta','q_BUCK']]
    db_fluxes = database[['sh','lh','tke','windstress']]
    
    
    #%% load and compute model data - atmostate
    theta_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_air_potential_temperature_24hrs_{stream_prefix}i_{flight}.nc', 'air_potential_temperature')[tstart_idx:tend_idx,:alt_idx]
    q_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_specific_humidity_24hrs_{stream_prefix}i_{flight}.nc', 'specific_humidity')[tstart_idx:tend_idx,:alt_idx]
    xwind_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_{xwind_name}_24hrs_{stream_prefix}h_{flight}.nc', xwind_name)[tstart_idx:tend_idx,:alt_idx]
    ywind_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_{ywind_name}_24hrs_{stream_prefix}h_{flight}.nc', ywind_name)[tstart_idx:tend_idx,:alt_idx]
    w_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_upward_air_velocity_24hrs_{stream_prefix}h_{flight}.nc', 'upward_air_velocity')[tstart_idx:tend_idx,:alt_idx]
    wsp_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_upward_air_velocity_24hrs_{stream_prefix}h_{flight}.nc', 'upward_air_velocity')[tstart_idx:tend_idx,:alt_idx]
    try:
        hoz_wind_data = (xwind_cube.data**2 + ywind_cube.data[:,:,:-1]**2)**0.5
    except:
        wsp_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_{xwind_name}_24hrs_{stream_prefix}h_{flight}.nc', xwind_name)[tstart_idx:tend_idx,:alt_idx]
        hoz_wind_data = (xwind_cube.data**2 + ywind_cube.data**2)**0.5
    wsp_cube.data = hoz_wind_data
    #%%
    wsp_cube.standard_name = None
    wsp_cube.long_name = 'windspeed'
    
    #%% load and compute model data - windstress
    
    ## windstress (two cumponents, need to be combined)
    # wse = windstrss eastward, wsn = windstress northward
    wse_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_atmosphere_downward_eastward_stress_24hrs_{stream_prefix}h_{flight}.nc', 'atmosphere_downward_eastward_stress')[tstart_idx:tend_idx,:alt_idx]
    wsn_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_atmosphere_downward_northward_stress_24hrs_{stream_prefix}h_{flight}.nc','atmosphere_downward_northward_stress')[tstart_idx:tend_idx,:alt_idx]
    # regrid wsn cube and combine cubes for magnitude
    #wsn_cube_reg = wsn_cube.regrid(wse_cube, iris.analysis.Linear())
    temp_cube =  iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_atmosphere_downward_eastward_stress_24hrs_{stream_prefix}h_{flight}.nc', 'atmosphere_downward_eastward_stress')[tstart_idx:tend_idx,:alt_idx]
    temp_cube.data = wsn_cube.data[:,:,:-1]
    
    ws_cube = (wse_cube**2 + temp_cube**2)**0.5
    ws_cube.standard_name = None
    ws_cube.long_name = 'windstress'
    #%%load and compute model data - heat fluxes
   
    # lh = latent heat, sh = sensible heat
    
    sh_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_upward_heat_flux_in_air_24hrs_{stream_prefix}i_{flight}.nc', 'upward_heat_flux_in_air')[tstart_idx:tend_idx,:alt_idx]
    lh_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_path}{config}_{res}_um_upward_latent_heat_flux_in_air_24hrs_{stream_prefix}i_{flight}.nc', 'upward_latent_heat_flux_in_air')[tstart_idx:tend_idx,:alt_idx]
    
    
    #%% rotate db coords
    
    polelat = w_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = w_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    rot_db_lon, rot_db_lat = rotate_pole(np.array(db_coors['lon']), np.array(db_coors['lat']), polelon, polelat)
    rot_db_lon = rot_db_lon + 360
    
    #%% adjust time coordinates
    
    ## db time is since 00:00 of same day, model time is since 1970
    ## convert db from econds in hours first
    db_time = np.array(db_coors['starttime'])/(60*60) + w_cube.coords('time')[0].points[0] - 0.5*(tstart_idx+1)
   
    time_coor = theta_cube.coords('time')[0]
    dates = time_coor.units.num2date(time_coor.points)
    print(dates)
    #%% interpolate coloumns
    db_alt = np.array(db_coors['altgps'])
    sample_points = [('grid_longitude', rot_db_lon), ('grid_latitude',rot_db_lat),('time', db_time) ]
    theta_cols = trajectory.interpolate(theta_cube, sample_points, method = 'nearest')
    q_cols = trajectory.interpolate(q_cube, sample_points, method = 'nearest')
    wsp_cols = trajectory.interpolate(wsp_cube, sample_points, method = 'nearest')
    w_cols = trajectory.interpolate(w_cube, sample_points, method = 'nearest')
    ws_cols = trajectory.interpolate(ws_cube, sample_points, method = 'nearest')
    lh_cols = trajectory.interpolate(lh_cube, sample_points, method = 'nearest')
    sh_cols = trajectory.interpolate(sh_cube, sample_points, method = 'nearest')
    print(theta_cols)
    #%%
    plt.plot(theta_cols.coords('grid_longitude')[0].points, theta_cols.coords('grid_latitude')[0].points)
    
    #%% interpolations
    theta_results = pert.coloumn_alt_interpolation(theta_cols, db_alt)
    q_results = pert.coloumn_alt_interpolation(q_cols, db_alt)
    wsp_results = pert.coloumn_alt_interpolation(wsp_cols, db_alt)
    w_results = pert.coloumn_alt_interpolation(w_cols, db_alt)
    ws_results = pert.coloumn_alt_interpolation(ws_cols, db_alt)
    lh_results = pert.coloumn_alt_interpolation(lh_cols, db_alt)
    sh_results = pert.coloumn_alt_interpolation(sh_cols, db_alt)
    
    #%% put data into pndas dataframe
    all_results = [theta_results,q_results,wsp_results,w_results,ws_results,lh_results,sh_results]
    print(np.shape(np.transpose(all_results)))
    #%%
    column_names = ['theta','q','wsp','w','ws','lh','sh']
    results_dataframe = pd.DataFrame(data = np.transpose(all_results),columns = column_names)
    print(results_dataframe)
    results_dataframe.to_csv(f'D:/Project/Model_Data/{suite}/{config}_Matched_interpolated_{res}_values_60s_{flight}.txt')