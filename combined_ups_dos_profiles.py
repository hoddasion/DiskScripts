# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 11:52:44 2022

@author: kse18nru

"Combined upstream-downstream profiles"
"""

#%% module imports
import iris
from iris.analysis.cartography import rotate_pole
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import xarray as xr
import equations
from iris.analysis import trajectory
#%% miscellaneous settings

knots_conv = 0.5144 # knots to m/s conversion factor
matplotlib.rcParams.update({'font.size': 24}) # set font size in figures unviversally


#%% load keflavik obs data
kefl_obs_path = 'D:/Project/Obvs_Data/profiles/'
kefl_filename = 'keflavik_04018BIKF_profiles_upto500mb.csv' # keflavik obs data file name (independent of case and suite since it's obs)
df_kefl = pd.read_csv(f'{kefl_obs_path}{kefl_filename}') # load keflavik obs data

obs_time = '20180319T1800Z'


kefl_profiles = df_kefl[df_kefl.Time == obs_time]
print(kefl_profiles)               
kefl_h = np.array(kefl_profiles.HGHT).astype(np.float)

kefl_p = np.array(kefl_profiles.PRES).astype(np.float)
kefl_theta = np.array(kefl_profiles.THTA).astype(np.float)
kefl_speed =np.array(kefl_profiles.SKNT).astype(np.float)*knots_conv # convert from knots to m/s
kefl_q = equations.specific_humidity(np.array(kefl_profiles.MIXR).astype(np.float))
print(kefl_q)
print(kefl_speed)
#%% load downwind MASIN data

MASIN_dataset = xr.open_dataset('../Obvs_Data/UEA_qc_MASIN_data/MASIN_flightdata_306.nc')
print(MASIN_dataset)

idx_0 = 13820; idx_1 = 14740
MASIN_lat = np.array(MASIN_dataset.lat_50hz)[::50][idx_0:idx_1]
MASIN_lon = np.array(MASIN_dataset.lon_50hz)[::50][idx_0:idx_1]
MASIN_alt = np.array(MASIN_dataset.gps_alt_50hz)[::50][idx_0:idx_1]
MASIN_temp = np.array(MASIN_dataset.air_temp_1hz)[idx_0:idx_1]
MASIN_pres = np.array(MASIN_dataset.static_pressure_1hz)[idx_0:idx_1]
MASIN_time = np.array(MASIN_dataset.time_1hz)[idx_0:idx_1]
MASIN_U = (np.array(MASIN_dataset.u_50hz)[::50][idx_0:idx_1]**2 + np.array(MASIN_dataset.v_50hz)[::50][idx_0:idx_1]**2)**0.5
#MASIN_q = np.array(MASIN_dataset.q)[idx_0:idx_1]
print(len(MASIN_alt))
print(MASIN_time[0], MASIN_time[-1])
MASIN_theta = equations.potential_temperature(MASIN_temp, MASIN_pres*100)
#%%
print(MASIN_alt)
#%% load model data
suite = 'u-cc134'
res = '1p5km'
um_path = f'D:/Project/Model_Data/{suite}/'
## velocity diagnostics
nwind_cube = iris.load_cube(f'{um_path}RA1M_{res}_um_y_wind_24hrs_ph_306.nc', 'y_wind')[:,:40]
ewind_cube = iris.load_cube(f'{um_path}RA1M_{res}_um_x_wind_24hrs_ph_306.nc', 'x_wind')[:,:40]
print(nwind_cube, ewind_cube)
magwind_cube = iris.load_cube(f'{um_path}RA1M_{res}_um_x_wind_24hrs_ph_306.nc', 'x_wind')[:,:40]
mag_data = (nwind_cube.data[:,:,:-1,:]**2 + ewind_cube.data**2)**0.5 
magwind_cube.data = mag_data
## theta diagnostic
theta_cube = iris.load_cube(f'{um_path}RA1M_{res}_um_air_potential_temperature_24hrs_pi_306.nc', 'air_potential_temperature')[:,:40]
q_cube = iris.load_cube(f'{um_path}RA1M_{res}_um_specific_humidity_24hrs_pi_306.nc', 'specific_humidity')[:,:40]

#%% isolate coloumns in model data

## get keflavik grid coordinates
kef_lat = 63.96
kef_lon = -22.60
polelat = theta_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
polelon = theta_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
kef_x, kef_y = rotate_pole(np.array([kef_lon]), np.array([kef_lat]), polelon, polelat)
kef_x = kef_x + 360
print(kef_y, kef_x)
## get MASIN grid coordinates
MASIN_x, MASIN_y = rotate_pole(np.array(MASIN_lon), np.array(MASIN_lat), polelon, polelat)
MASIN_x = MASIN_x + 360

## set sample points (kefl first, then MASIN at start of ascent)
sample_points = [('grid_longitude', [kef_x[0],MASIN_x[0]]), ('grid_latitude', [kef_y[0], MASIN_y[0]])]
## nearest neighbour interpolate model cubes to get coloumns
theta_coloumns = trajectory.interpolate(theta_cube,sample_points, method =  'nearest')
wind_coloumns = trajectory.interpolate(magwind_cube,sample_points, method =  'nearest')
q_coloumns = trajectory.interpolate(q_cube,sample_points, method =  'nearest')
print(theta_coloumns)

#%% plot figure
fig, (ax0,ax1,ax2) = plt.subplots(1,3, figsize = (14,10), sharey = True)
time_idx = 35
top_alt = 3000
top_idx = 30
## assign U to ax0
ax0.plot(kefl_speed, kefl_h, linestyle = 'dashed', color = 'tab:blue') # keflavik profile
ax0.plot(MASIN_U, MASIN_alt, linestyle = 'dashed', color = 'tab:orange') # MASIN profile
ax0.plot(wind_coloumns[time_idx,:top_idx,0].data, wind_coloumns.coords('altitude')[0].points[:top_idx,0], color = 'tab:blue', label = 'Keflavik')
ax0.plot(wind_coloumns[time_idx,:top_idx,1].data, wind_coloumns.coords('altitude')[0].points[:top_idx,1], color = 'tab:orange', label ='downwind')
ax0.set_ylim(bottom = 0, top = top_alt)
ax0.legend(fontsize = 20)
ax0.set_xticks([0,5,10,15,20])
ax0.set_xlabel(r'$U$, ms$^{-1}$')
ax0.grid()
ax0.set_yticks([500,1000,1500,2000,2500,3000])
ax0.set_yticklabels(['0.5','1.0','1.5','2.0','2.5','3.0'])
ax0.set_ylabel('Altitude, km')
## assign theta to ax1
ax1.plot(kefl_theta, kefl_h, linestyle = 'dashed', color = 'tab:blue') # keflavik profile
ax1.plot(MASIN_theta, MASIN_alt, linestyle = 'dashed', color = 'tab:orange') # MASIN profile
ax1.plot(theta_coloumns[time_idx,:top_idx,0].data, theta_coloumns.coords('altitude')[0].points[:top_idx,0], color = 'tab:blue')
ax1.plot(theta_coloumns[time_idx,:top_idx,1].data, theta_coloumns.coords('altitude')[0].points[:top_idx,1], color = 'tab:orange')
ax1.set_ylim(bottom = 0, top = top_alt)
ax1.set_xlim(right = 295)
ax1.set_xticks([276,280,284,288,292])
ax1.set_xlabel(r'$\theta$, K')
ax1.grid()

## aasign q to ax2
ax2.plot(kefl_q, kefl_h, linestyle = 'dashed', color = 'tab:blue') # keflavik profile
#ax0.plot(MASIN_U, MASIN_alt, linestyle = 'dashed', color = 'tab:orange') # MASIN profile
ax2.plot(q_coloumns[time_idx,:top_idx,0].data*1000, q_coloumns.coords('altitude')[0].points[:top_idx,0], color = 'tab:blue')
ax2.plot(q_coloumns[time_idx,:top_idx,1].data*1000, q_coloumns.coords('altitude')[0].points[:top_idx,1], color = 'tab:orange')
ax2.set_ylim(bottom = 0, top = top_alt)
ax2.set_xlabel(r'$q$, gkg$^{-1}$')
ax2.set_xticks([0,1,2,3,4])
ax2.grid()

fig.suptitle('Case 2 Profiles at 1800 UTC')

plt.tight_layout()

plt.savefig('D:/Project/Figures/PNG/306/CONTROL/u-cc134/profiles/RA1M_obs_upwind_downwind_profiles_1800UTC.png')