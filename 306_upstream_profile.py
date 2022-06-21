# -*- coding: utf-8 -*-
"""
Created on Thu Oct  7 14:12:01 2021

@author: kse18nru
"""
#%%
import xarray as xr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import iris
from iris.analysis.cartography import rotate_pole
import equations
from iris.analysis import trajectory
import pandas as pd
#%% load model data

um_path = 'D:/Project/Model_Data/u-cc134/'
land_cube = iris.load_cube(f'{um_path}RA1M_1p5km_um_land_binary_mask_24hrs_pi_306.nc')
land_data = land_cube.data
land_x = land_cube.coord('grid_longitude').points
land_y = land_cube.coord('grid_latitude').points

theta_cube = iris.load_cube(f'{um_path}RA1M_1p5km_um_air_potential_temperature_24hrs_pi_306.nc', 'air_potential_temperature')
theta4p4_cube = iris.load_cube(f'{um_path}RA1M_4p4km_um_air_potential_temperature_24hrs_pi_306.nc', 'air_potential_temperature')
#%% open MASIN data
obs_dataset = xr.open_dataset('../Obvs_Data/UEA_qc_MASIN_data/MASIN_flightdata_306.nc')
print(obs_dataset)

#%%
idx_0 = 13820; idx_1 = 14740
obs_lat = np.array(obs_dataset.lat_50hz)[::50][idx_0:idx_1]
obs_lon = np.array(obs_dataset.lon_50hz)[::50][idx_0:idx_1]
obs_alt = np.array(obs_dataset.gps_alt_50hz)[::50][idx_0:idx_1]
obs_temp = np.array(obs_dataset.air_temp_1hz)[idx_0:idx_1]
obs_pres = np.array(obs_dataset.static_pressure_1hz)[idx_0:idx_1]
obs_time = np.array(obs_dataset.time_1hz)[idx_0:idx_1]
print(len(obs_alt))
print(obs_time[0], obs_time[-1])
obs_theta = equations.potential_temperature(obs_temp, obs_pres*100)

#%%
print(f'lon: {obs_lon[-1]-obs_lon[0]}')
print(f'lat: {np.mean(obs_lat)}')
#%% rotate obs coordinates onto model grid

polelat = land_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
polelon = land_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
obs_x, obs_y = rotate_pole(np.array(obs_lon), np.array(obs_lat), polelon, polelat)
obs_x = obs_x + 360

#%% interpolate onto two coloumns
sample_points = [('grid_longitude', [obs_x[0], obs_x[-1]]), ('grid_latitude', [obs_y[0], obs_y[-1]])]
theta_coloumns = trajectory.interpolate(theta_cube,sample_points, method =  'nearest')
theta_coldata = theta_coloumns.data
coloumns_alt = theta_coloumns.coord('altitude').points
theta4p4_cols = trajectory.interpolate(theta4p4_cube,sample_points, method = 'nearest')
cols4p4_alt = theta4p4_cols.coord('altitude').points
#%%
print(theta_coloumns)
print(np.shape(coloumns_alt[:,0]))
#%%

fig, ax = plt.subplots(1,1, figsize = (10,10))
ax.contour(land_x, land_y, land_data)
ax.scatter(obs_x, obs_y, c = obs_alt)
#%% save profile data to new csv
## make dataframe
save_obs = False
if save_obs:
    df_data = {'gridx':obs_x,'gridy':obs_y,'lon':obs_lon,'lat':obs_lat,'alt':obs_alt,'airtemp_1hz':obs_temp,'pressure_1hz':obs_pres, 'theta_1hz':obs_theta, 'time':obs_time}
    obs_df = pd.DataFrame(df_data)
    print(obs_df)
    obs_df.to_csv('../Obvs_Data/profiles/MASIN_1hz_downstream_profile_ft306.csv')
#%% plot profiles
fig0, axp1 = plt.subplots(1,1, figsize = (10,14))

#axp0.plot(obs_x, obs_alt)
axp1.plot(obs_theta, obs_alt, label = 'MASIN 1hz obs', linewidth = 2, color = 'k')
axp1.plot(theta_coldata[33,:29,0], coloumns_alt[:29,0], label = 'RA1M 1p5km at segment start')
axp1.plot(theta_coldata[33,:29,1], coloumns_alt[:29,0], label = 'RA1M 1p5km at segment end')
axp1.plot(theta4p4_cols.data[33,:29,0], cols4p4_alt[:29,0], label = 'RA1M 4p4km at segment start')
axp1.plot([270,300], [1500,1500], color = 'k', linestyle = 'dashed')
axp1.legend(fontsize = 20)
axp1.set_xlim(274,290)
axp1.set_ylim(bottom = 0, top = 2700)
axp1.set_ylabel('Altitude, m')
axp1.set_xlabel(r'$\theta$, K')
plt.tight_layout()
#plt.savefig('../Figures/PNG/306/u-cc134/profiles/End_of_flight_profile_theta_1hz_306.png')
#plt.savefig('../Figures/PDF/306/u-cc134/profiles/End_of_flight_profile_theta_1hz_306.pdf')