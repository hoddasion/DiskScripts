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
    theta_um_mean = theta_cube.collapsed('time',iris.analysis.MEAN)
    
    #%%
    w_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{config}_{res}_um_upward_air_velocity_24hrs_ph_{flight}.nc', 'upward_air_velocity')[tstart_idx:tend_idx,:2]
    polelat = w_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = w_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
    rot_db_lon = rot_db_lon + 360
    print(polelat, polelon)
    print(rot_db_lat, db_lat)
    #print(theta_um_mean)
    
    
    #%%
    tmax = 300
    tmin = 275
    #%% plotting
    fig, ax  = plt.subplots(1,1, figsize = (20,10))
    iplt.contourf(theta_um_mean, coords = ['grid_latitude', 'altitude'], levels = 30, vmin = tmin, vmax = tmax)
    ax.set_ylim(top = 5000)
    ax.scatter(rot_db_lat, db_alt, s = 1000, c = db_theta, vmin = tmin, vmax = tmax)
    
    