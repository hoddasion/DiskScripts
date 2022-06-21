# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 15:36:43 2021

@author: kse18nru
"""

#%% module imports
import sys
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
import model_foundry_toolkit as toolkit
import data_management as dmna
import time
import model_foundry as foundry
import matplotlib.gridspec as gridspec
import iris.analysis.cartography
import cartopy.crs as ccrs
import glm_functions as glmf
import datetime
from iris.analysis import trajectory

#%% global variables
plt.rcParams.update({'font.size': 20})
flight = 306
suite = 'u-cc134'
variable = ''
kconstraints = {'longitude' :(-25,-17), 'latitude' :(64,67)} # coordinate constraints for iris intersection function
#%%load glm data

cubes = iris.load(f'D:/Project/Model_Data/{suite}/RA1M_glm_1_flt{flight}.nc')
#%%
print(cubes)
#%%
cubes = None
#%%
temp19_cube = glmf.glm_load_and_concat('m01s03i236', kconstraints)
#rh_cube = glmf.glm_load_and_concat('m01s03i245', kconstraints)
lbm_cube = glmf.glm_load_and_concat('m01s00i030', kconstraints, fileno = 1)
sa_cube = glmf.glm_load_and_concat('m01s00i033', kconstraints, fileno = 1)



#%% load regional data

#%% load gsta obs
obs_path = 'D:/Project/Obvs_Data/groundstations/weather_obs_201803.csv'
stations_path = 'D:/Project/Obvs_Data/groundstations/weather_stations_201803.csv'

obs_df = pd.read_csv(obs_path)
obs_df['TIMESTAMP'] = pd.to_datetime((obs_df['TIMESTAMP'] - obs_df.TIMESTAMP[0])*24, unit = 'h', origin =pd.Timestamp('2018-03-01'))
obs_df.rename(columns={'T': 'TEMP'}, inplace=True)
gsta_df = pd.read_csv(stations_path) # gsta: groundstations
gsta_lon = np.array(gsta_df.LON)
gsta_lat = np.array(gsta_df.LAT)
gsta_sname = np.array(gsta_df.SHORTNAME)
obs_df.TIMESTAMP = obs_df.TIMESTAMP.round('H')


#print(gsta_df)
#%%
print(obs_df)
print(gsta_df)

#%%
lower19_time = obs_df[obs_df.TIMESTAMP <= datetime.datetime(2018,3,20,0,0)]
obs_19df = lower19_time[lower19_time.TIMESTAMP > datetime.datetime(2018,3,19,0,0)]
tstamp19 = obs_19df.TIMESTAMP#.round('H')#.replace(second = 0, microsecond = 0, minute = 0)
print(tstamp19[tstamp19 == datetime.datetime(2018,3,20,0)])
#obs_time = unit.date2num(np.array(obs_19df.TIMESTAMP)[0], 'hours since 1970-01-01 00:00:00',unit.CALENDAR_STANDARD)
#%% subset obs to model timesteps
obs19_0900 = obs_19df[tstamp19 == datetime.datetime(2018,3,19,9)]
obs19_1200 = obs_19df[tstamp19 == datetime.datetime(2018,3,19,12)]
obs19_1500 = obs_19df[tstamp19 == datetime.datetime(2018,3,19,15)]
obs19_1800 = obs_19df[tstamp19 == datetime.datetime(2018,3,19,18)]
obs19_2100 = obs_19df[tstamp19 == datetime.datetime(2018,3,19,21)]
obs19_0000 = obs_19df[tstamp19 == datetime.datetime(2018,3,20,0)]


#%% loop through and identify missing stations
print(gsta_df)
c, int_idx, uni_idx = np.intersect1d(np.array(gsta_df.STATION), np.array(obs19_0900.STATION), assume_unique = True, return_indices = True)


print(int_idx)
lon_reduced = [gsta_lon[i] for i in int_idx]
lat_reduced = [gsta_lat[i] for i in int_idx]
alt_reduced = [gsta_df.HEIGHT[i] for i in int_idx]
print(alt_reduced)
#%% interpolate model data by nearet neighbour
sample_points = [('longitude',lon_reduced),('latitude',lat_reduced)]
temp19_inter = trajectory.interpolate(temp19_cube,sample_points, method =  'nearest') # air temperature
#rh_inter = trajectory.interpolate(rh19_cube,sample_points, method =  'nearest') # relative humidity
sa_inter = trajectory.interpolate(sa_cube, sample_points, method = 'nearest') # surface altitude
print(sa_inter)


#%% T correction
model_alt = sa_inter.data
## delta = model - obs
delta_alt = model_alt - alt_reduced
print(delta_alt)
padia_rate = 6/1000 # pseudo-adiabatic lapse rate ~ 6K/km
dT = delta_alt*padia_rate
T19_model_cor = temp19_inter.data + dT
#%% plot boxplots etc
fig, (ax0,ax1) = plt.subplots(2,1, figsize = (14,10))
widths = 0.5
## temperature
ax0.boxplot((np.array(obs19_0900.TEMP)+273.15,T19_model_cor[0]), positions = [1.2,1.8], widths = widths)
ax0.boxplot((np.array(obs19_1200.TEMP)+273.15, T19_model_cor[1]), positions = [3.2,3.8], widths = widths)
ax0.boxplot((np.array(obs19_1500.TEMP)+273.15, T19_model_cor[2]), positions = [5.2,5.8], widths = widths)
ax0.boxplot((np.array(obs19_1800.TEMP)+273.15,T19_model_cor[3]), positions = [7.2,7.8], widths = widths)
ax0.boxplot((np.array(obs19_2100.TEMP)+273.15, T19_model_cor[4]), positions = [9.2,9.8], widths = widths)
ax0.boxplot((np.array(obs19_0000.TEMP)+273.15, T19_model_cor[5]), positions = [11.2,11.8], widths = widths)
## plot means as timeseries on taylored x-axis array
x_vals = [1.5,3.5,5.5,7.5,9.5,11.5]
x_labels = ['0900','1200','1500','1800','2100','0000']
ax0.set_xticks(x_vals); ax0.set_xticklabels([])#; ax0.set_xticklabels(x_labels)
ax0.set_ylabel('Air temperature, K')
ax0.grid(True)

## relative humidity
# ax1.boxplot((np.array(obs_0900.RH)[~np.isnan(obs_0900.RH)], rh_inter.data[0][~np.isnan(obs_0900.RH)]), positions = [1.2,1.8], widths = widths)
# ax1.boxplot((np.array(obs_1200.RH)[~np.isnan(obs_1200.RH)], rh_inter.data[1][~np.isnan(obs_1200.RH)]), positions = [3.2,3.8], widths = widths)
# ax1.boxplot((np.array(obs_1500.RH)[~np.isnan(obs_1500.RH)], rh_inter.data[2][~np.isnan(obs_1500.RH)]), positions = [5.2,5.8], widths = widths)
# ax1.boxplot((np.array(obs_1800.RH)[~np.isnan(obs_1800.RH)], rh_inter.data[3][~np.isnan(obs_1800.RH)]), positions = [7.2,7.8], widths = widths)
# ax1.boxplot((np.array(obs_2100.RH)[~np.isnan(obs_2100.RH)], rh_inter.data[4][~np.isnan(obs_2100.RH)]), positions = [9.2,9.8], widths = widths)
# ax1.boxplot((np.array(obs_0000.RH)[~np.isnan(obs_0000.RH)], rh_inter.data[5]), positions = [11.2,11.8], widths = widths)
ax1.set_xticks(x_vals); ax1.set_xticklabels(x_labels)
ax1.set_ylabel('Relative humidity, %')
ax1.set_xlabel('Time UTC')
ax1.grid(True)
fig.suptitle('Case 2 - Groundstation obs vs GA6.1 N768')

plt.tight_layout()

plt.savefig('D:/Project/Figures/PNG/306/u-cc134/stat_val/Case1_2_GA61_GSTAobs_temp_20180319_padiacor_boxplots.png')



#%% start plotting
data_crs = ccrs.PlateCarree()
gs = gridspec.GridSpec(nrows = 3, ncols = 4)
projtype = ccrs.NorthPolarStereo() 
splt_kwargs = {'projection' : projtype, 'frameon' : False} # subplot keyword arguments
lbm_kwargs = {'colors':'k','transform':data_crs} # land mask keyword arguments

fig = plt.figure(figsize = (18,13))
    
ax0 = fig.add_subplot(gs[0,0],**splt_kwargs)
ax0.contour(lbm_lon, lbm_lat, lbm_data, colors = 'k', levels = [1], linewidths = 0.5)#, **lbm_kwargs)
lon_mesh, lat_mesh = np.meshgrid(lbm_lon, lbm_lat)
#ax0.scatter(lon_mesh, lat_mesh, transform = data_crs)
ax0.pcolormesh(temp_lon, temp_lat, temp_data[0])#, transform = data_crs)
ax0.scatter(gsta_lon, gsta_lat, marker = 'x', color = 'k')
