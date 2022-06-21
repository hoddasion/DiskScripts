# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 17:20:09 2021

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

#%% model land mask loading
um_path = 'D:/Project/Model_Data/u-cc134/'
land_cube = iris.load_cube(f'{um_path}RA1M_1p5km_um_land_binary_mask_24hrs_pi_306.nc')
land_data = land_cube.data
land_x = land_cube.coord('grid_longitude').points
land_y = land_cube.coord('grid_latitude').points

#%% load obs data
## gstas
gsta_306_df_all = pd.read_csv('D:/Project/Obvs_Data/groundstations/weather_stations_201803.csv')
gsta_301_df = pd.read_csv('D:/Project/Obvs_Data/groundstations/snaesfellsnes_stations_info.csv')

## sample correct 306 stations
station_codes= np.array([2862,2655,2692,2738])
gsta_condition = np.array(gsta_306_df_all.STATION) == 0
for scode in station_codes:
    temp_cond =  np.array(gsta_306_df_all.STATION) == scode
    gsta_condition = gsta_condition + temp_cond
gsta_306_df = gsta_306_df_all[gsta_condition]  
print(gsta_301_df,'\n', gsta_306_df)
## gt gsta snames
gsta_301_sname = gsta_301_df.SHORTNAME
gsta_306_sname = gsta_306_df.SHORTNAME
#%%
### flights
flight306_filename = 'IGP_flights_database_obs_60s_306.txt'
flight301_filename = 'IGP_flights_database_obs_60s_301.txt'
path = 'D:/Project/Obvs_Data/Databases/'
flight_306_df = pd.read_csv(f'{path}{flight306_filename}',delimiter = ' ')
flight_301_df = pd.read_csv(f'{path}{flight301_filename}',delimiter = ' ')

lon_306 = flight_306_df.lon
lat_306 = flight_306_df.lat
lon_301 = flight_301_df.lon
lat_301 = flight_301_df.lat

#%% down stream profile flt 306
ds_profile_df = pd.read_csv('D:/Project/Obvs_Data/profiles/MASIN_1hz_downstream_profile_ft306.csv')
ds_x = ds_profile_df.gridx
ds_y = ds_profile_df.gridy
#%% rotate obs coords onto model grid 
polelat = land_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
polelon = land_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
gsta_306_x, gsta_306_y = rotate_pole(np.array(gsta_306_df.LON), np.array(gsta_306_df.LAT), polelon, polelat)
gsta_306_x = gsta_306_x + 360
gsta_301_x, gsta_301_y = rotate_pole(np.array(gsta_301_df.LON), np.array(gsta_301_df.LAT), polelon, polelat)
gsta_301_x = gsta_301_x + 360
kefl_x, kefl_y = rotate_pole(np.array([-22.60]), np.array([63.96]), polelon, polelat)
kefl_x = kefl_x + 360
case2x, case2y = rotate_pole(np.array(lon_306), np.array(lat_306), polelon, polelat)
case2x = case2x + 360
case1x, case1y = rotate_pole(np.array(lon_301), np.array(lat_301), polelon, polelat)
case1x = case1x + 360
#%% plot map

fig, ax = plt.subplots(1,1, figsize = (8,12))

## set up land borders and frame
ax.contour(land_x, land_y, land_data, colors = 'grey', linewidths = 0.2)
ax.set_xlim(left = 359, right = 361.4)
ax.set_ylim(bottom = -2.5, top = 1.5)
ax.set_xticks([])
ax.set_yticks([])

ax.plot(case1x,case1y, c = 'blue', label = 'runaveraged aircraft')
ax.plot(case2x,case2y, c = 'blue')
ax.plot(ds_x,ds_y, c= 'olive', linewidth = 4, label = '1hz aircraft, for profile')
## add ground stations
ax.scatter(gsta_301_x,gsta_301_y, color = 'red', label = 'groundstations') # 301 stations
for i,txt in enumerate(np.array(gsta_301_sname)):
        ax.annotate(txt, (gsta_301_x[i],gsta_301_y[i]), xytext = (gsta_301_x[i] - 0.065,gsta_301_y[i]+0.05), fontsize = 12)
ax.scatter(gsta_306_x,gsta_306_y, color = 'red') # 306 stations
for i,txt in enumerate(np.array(gsta_306_sname)):
        ax.annotate(txt, (gsta_306_x[i],gsta_301_y[i]), xytext = (gsta_306_x[i] - 0.065,gsta_306_y[i]+0.03), fontsize = 12)
        
## plot keflavik location
ax.scatter(kefl_x, kefl_y, color = 'purple', marker = 'x', s = 100)
ax.annotate('Keflavik', (kefl_x[0], kefl_y[0]), xytext = (kefl_x[0] - 0.15, kefl_y[0] + 0.04))
ax.annotate('IGP 306', (case2x[0] +0.2,case2y[0]+0.5))
ax.annotate('IGP 301', (case1x[0]-0.8,case1y[0]+0.40))
plt.savefig('obs_source_map.png')