# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 13:35:09 2021

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
#%% 
res = '0p5km'

suite = 'u-cc134'
config = 'RA1M'
alt_idx = 30
tstart_idx = 24; tend_idx = 40
flight = 306



#%%
if flight == 306:
    batch = [1,4,6,11,13]
    #%% load obs data
    obs_filename = f'IGP_flights_database_obs_60s_{flight}.txt'
    path = 'D:/Project/Obvs_Data/Databases/'
    
    df_obs = pd.read_csv(f'{path}{obs_filename}', delimiter = ' ')
    db_legno = np.array(df_obs['legno'])
    db_theta = np.array(df_obs['theta'])
    db_lon = np.array(df_obs['lon'])
    db_lat = np.array(df_obs['lon'])
    db_ws = np.array(df_obs['windstress'])
    db_sh = np.array(df_obs['sh'])
    db_lh = np.array(df_obs['lh'])
    db_time = np.array(df_obs['starttime'])
    db_alt = np.array(df_obs['altgps'])
    #%%
    
    df_avrs = None
    for leg in batch:
        print(leg)
        if leg == 1:
            df_avrs = pd.DataFrame(df_obs[db_legno == leg].mean(axis = 0)).transpose()
            
        else:
            df_avrs_temp = pd.DataFrame(df_obs[db_legno == leg].mean(axis = 0)).transpose()
            #df_avrs.append(df_avrs_temp, ignore_index = True)
            df_avrs = pd.concat([df_avrs, df_avrs_temp], ignore_index = True)
    df_avrs.set_index('legno', inplace = True)
    print(df_avrs['theta'])
    #%%
    db_theta = np.array(df_avrs['theta'])
    db_lon = np.array(df_avrs['lon'])
    db_lat = np.array(df_avrs['lat'])
    db_ws = np.array(df_avrs['windstress'])
    db_sh = np.array(df_avrs['sh'])
    db_lh = np.array(df_avrs['lh'])
    db_time = np.array(df_avrs['starttime'])
    db_alt = np.array(df_avrs['altgps'])
    #%%
    theta_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/air_potential_temperature/{config}_vertslice_{res}_air_potential_temperature_flt306_leg7.nc', 'air_potential_temperature')[tstart_idx:tend_idx,:40]
    rot_theta_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/air_potential_temperature/{config}_vertslice_{res}_air_potential_temperature_flt306_leg7.nc', 'air_potential_temperature')[tstart_idx:tend_idx,:40]
    surf_alt_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/air_potential_temperature/{config}_vertslice_{res}_air_potential_temperature_flt306_leg7.nc', 'surface_altitude')
    print(surf_alt_cube)
    #%%
    w_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{config}_{res}_um_upward_air_velocity_24hrs_ph_{flight}.nc', 'upward_air_velocity')[tstart_idx:tend_idx,:2]
    polelat = w_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = w_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    grid_db_lon, grid_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
    grid_db_lon = grid_db_lon + 360
    
    #print(theta_um_mean)
    
    
    #%% coordinate system rotation
    ## need phi = angle between transect and y-axis (grid_latitude)
    um_x = theta_cube.coords('grid_longitude')[0].points
    um_y = theta_cube.coords('grid_latitude')[0].points
    print(um_x)
    dy = um_y[-1] - um_y[0]
    dx = um_x[-1] - um_x[0]
    phi = np.arctan(dx/dy)
    print(phi)
    
    ## peform for model
    rot_um_x = um_x*np.cos(-phi) + um_y*np.sin(-phi)
    rot_um_y = -um_x*np.sin(-phi) + um_y*np.cos(-phi)
    rot_theta_cube.coords('grid_longitude')[0].points = rot_um_x
    rot_theta_cube.coords('grid_latitude')[0].points = rot_um_y
    
    ## perform for obs
    rot_obs_x = grid_db_lon*np.cos(-phi) + grid_db_lat*np.sin(-phi)
    rot_obs_y = -grid_db_lon*np.sin(-phi) + grid_db_lat*np.cos(-phi)
    
    #%% collapse cube in time
    theta_um_mean = rot_theta_cube.collapsed('time',iris.analysis.MEAN)
    
    #%%
    tmax = 290 #np.max(db_theta)
    tmin =275 #np.min(db_theta)
    #%% plotting
    fig, ax  = plt.subplots(1,1, figsize = (10,4))
    iplt.contourf(theta_um_mean, coords = ['grid_latitude', 'altitude'], levels = 30, vmin = tmin, vmax = tmax, cmap = 'plasma')
    ax.plot(rot_um_y,surf_alt_cube.data, color = 'black' )
    ax.set_ylim(top = 2000, bottom = 0)
    ax.set_xlim(right = 246.6, left = 246.1)
    ax.scatter(rot_obs_y, db_alt, s = 1000, c = db_theta, vmin = tmin, vmax = tmax, cmap = 'plasma', edgecolor = 'white')
    
    