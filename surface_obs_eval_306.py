# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 13:17:13 2021

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
import cf_units as unit
import datetime 

#%%
obs_path = 'D:/Project/Obvs_Data/groundstations/weather_obs_201803.csv'
stations_path = 'D:/Project/Obvs_Data/groundstations/weather_stations_201803.csv'

obs_df = pd.read_csv(obs_path)
obs_df['TIMESTAMP'] = pd.to_datetime((obs_df['TIMESTAMP']- obs_df.TIMESTAMP[0])*24, unit = 'h', origin =pd.Timestamp('2018-03-01'))
obs_df.rename(columns={'T': 'TEMP'}, inplace=True)
gsta_df = pd.read_csv(stations_path) # gsta: groundstations
gsta_lon = np.array(gsta_df.LON)
gsta_lat = np.array(gsta_df.LAT)
gsta_sname = np.array(gsta_df.SHORTNAME)
print(gsta_df)
print(obs_df[obs_df.STATION == 2738])

#%%
lower_time = obs_df[obs_df.TIMESTAMP <= datetime.datetime(2018,3,20,0,0)]
obs_19df = lower_time[lower_time.TIMESTAMP > datetime.datetime(2018,3,19,0,0)]
print(obs_19df.TIMESTAMP)
#obs_time = unit.date2num(np.array(obs_19df.TIMESTAMP)[0], 'hours since 1970-01-01 00:00:00',unit.CALENDAR_STANDARD)
#%%

station_codes= np.array([2862,2655,2692,2738])

## gsta dataframe first
gsta_condition = np.array(gsta_df.STATION) == 0
for scode in station_codes:
    temp_cond =  np.array(gsta_df.STATION) == scode
    gsta_condition = gsta_condition + temp_cond
four_gsta = gsta_df[gsta_condition]   
gsta_lon = gsta_lon[gsta_condition]
gsta_lat = gsta_lat[gsta_condition]
gsta_sname = gsta_sname[gsta_condition]
## now obs dataframe
obs_condition = np.array(obs_df.STATION) == 0
for scode in station_codes:
    temp_cond = np.array(obs_df.STATION) == scode
    obs_condition = obs_condition + temp_cond
obs_df = obs_df[obs_condition]

#%% subset obs data into stations 
HBV_cond = np.array(obs_19df.STATION) == 2862
HBV_obs = obs_19df[HBV_cond]
AEDEY_cond = np.array(obs_19df.STATION) == 2655
AEDEY_obs = obs_19df[AEDEY_cond]
GJOGR_cond = np.array(obs_19df.STATION) == 2692
GJOGR_obs = obs_19df[GJOGR_cond]
BOLUN_cond = np.array(obs_19df.STATION) == 2738
BOLUN_obs = obs_19df[BOLUN_cond]

#%% model data
um_path = 'D:/Project/Model_Data/u-cc134/'
rh_cube = iris.load_cube(f'{um_path}RA1M_0p5km_um_relative_humidity_24hrs_pg_306.nc')
land_cube = iris.load_cube(f'{um_path}RA1M_0p5km_um_land_binary_mask_24hrs_pi_306.nc')
temp_cube = iris.load_cube(f'{um_path}RA1M_0p5km_um_air_temperature_24hrs_pg_306.nc')
xwps_cube = iris.load_cube(f'{um_path}RA1M_0p5km_um_x_wind_24hrs_pg_306.nc')
ywps_cube = iris.load_cube(f'{um_path}RA1M_0p5km_um_y_wind_24hrs_pg_306.nc')
wsp_cube = iris.load_cube(f'{um_path}RA1M_0p5km_um_x_wind_24hrs_pg_306.nc')
wps_data = (xwps_cube.data**2 + ywps_cube.data[:,:-1]**2)**0.5
wsp_cube.data = wps_data

#%% convert time to datetimes in cubes
print(unit.num2date(temp_cube.coord('time').points, 'hours since 1970-01-01 00:00:00',unit.CALENDAR_STANDARD))
print(temp_cube.coord('time').points)
#%% rotate obs coords onto model grid 
polelat = rh_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
polelon = rh_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
gsta_x, gsta_y = rotate_pole(np.array(gsta_lon), np.array(gsta_lat), polelon, polelat)
gsta_x = gsta_x + 360
#%% map
mapofstations =True
if mapofstations:
    fig, ax = plt.subplots(1,1, figsize = (6.5,10))
    ax.contour(land_cube.coords('grid_longitude')[0].points,land_cube.coords('grid_latitude')[0].points,land_cube.data, colors = 'k')
    ax.scatter(gsta_x,gsta_y)
    for i,txt in enumerate(gsta_sname):
        ax.annotate(txt, (gsta_x[i],gsta_y[i]))
        
#%% Model subsetting
## for temperature
model_x = []; model_y = []
model_x.append(gsta_x[gsta_sname == 'HBV'][0]); model_y.append(gsta_y[gsta_sname == 'HBV'][0])
model_x.append(gsta_x[gsta_sname == 'AEDEY'][0]); model_y.append(gsta_y[gsta_sname == 'AEDEY'][0])
model_x.append(gsta_x[gsta_sname == 'GJOGR'][0]); model_y.append(gsta_y[gsta_sname == 'GJOGR'][0])
model_x.append(gsta_x[gsta_sname == 'BOLUN'][0]); model_y.append(gsta_y[gsta_sname == 'BOLUN'][0])
sample_points =  [('grid_longitude',np.array(model_x)), ('grid_latitude', np.array(model_y))]
print(model_x, model_y)
temp_subcube = trajectory.interpolate(temp_cube,sample_points, method =  'nearest') # temp meaning temperature, not temporary
wsp_subcube = trajectory.interpolate(wsp_cube,sample_points, method =  'nearest')
rh_subcube = trajectory.interpolate(rh_cube,sample_points, method =  'nearest')

#%% calculate T corrections
calc_Tcorrections = True
if calc_Tcorrections:
    print(four_gsta)
    #%%
    Tht3D_cube = iris.load_cube(f'{um_path}RA1M_0p5km_um_air_potential_temperature_24hrs_pi_306.nc', 'air_potential_temperature')
    Tht3D_subcube  = trajectory.interpolate(Tht3D_cube,sample_points, method =  'nearest') 
    model_z = Tht3D_subcube.coord('surface_altitude').points
    print(Tht3D_subcube.coord('surface_altitude'))
    
    alt_HBV = 22.0; alt_AEDEY = 21.0; alt_GJOGR = 31.0; alt_BOLUN = 27.0
    ## recall index 0->HBV, 1->AEDEY, 2->GJOGR, 3->BOLUN
    ## dz = dz(um) - dz(obs)
    dz_HBV = model_z[0] - alt_HBV;
    dz_AEDEY = model_z[1] - alt_AEDEY; 
    dz_GJOGR = model_z[2] - alt_GJOGR; 
    dz_BOLUN = model_z[3] - alt_BOLUN; 
    padia_rate = 6/1000 # pseudo-adiabatic lapse rate ~ 6K/km
    dT_HBV = dz_HBV*padia_rate; 
    dT_AEDEY = dz_AEDEY*padia_rate
    dT_GJOGR = dz_GJOGR*padia_rate
    dT_BOLUN = dz_BOLUN*padia_rate
    print(model_z[0], dz_HBV, dT_HBV)
    print(model_z[1], dz_AEDEY, dT_AEDEY)
    print(model_z[2], dz_GJOGR, dT_GJOGR)
    print(model_z[3], dz_BOLUN, dT_BOLUN)
#%%
print(wsp_subcube[:,0].data)
print(HBV_obs.WSP)

#%% make time (x-axis) list
time_x = np.arange(1,49)*0.5
print(time_x[1::2])

#%% plot temperature
plot_temp = True
if plot_temp:
    lowT = 271.5; highT = 282
    fig, (ax0,ax1,ax2,ax3) = plt.subplots(4,1, figsize = (10,10))
    CtK = 273.15 # Celsius to Kelvin conversion
    ax0.plot(time_x[1::2], np.array(HBV_obs.TEMP) + CtK, label = 'obs')
    ax0.plot(time_x, temp_subcube.data[:,0], label = 'UM 0p5km')
    ax0.set_xlim(left = 0, right = 24)
    ax0.set_xticks(np.arange(0,25))
    ax0.set_xticklabels([])
    ax0.set_ylim(bottom = lowT, top = highT)
    ax0.set_yticks([272,274,276,278,280])
    ax0.set_title('HBV - Most northern, downstream coast')
    ax0.legend(fontsize = 12)
    ax0.set_ylabel('K')
    
    ax1.plot(time_x[1::2], np.array(AEDEY_obs.TEMP) + CtK)
    ax1.plot(time_x, temp_subcube.data[:,1])
    ax1.set_xlim(left = 0, right = 24)
    ax1.set_xticks(np.arange(0,25))
    ax1.set_xticklabels([])
    ax1.set_ylim(bottom = lowT, top = highT)
    ax1.set_yticks([272,274,276,278,280])
    ax1.set_title('AEDEY - Upstream coast on major fjord')
    ax1.set_ylabel('K')
    
    ax2.plot(time_x[1::2], np.array(GJOGR_obs.TEMP) + CtK)
    ax2.plot(time_x, temp_subcube.data[:,2])
    ax2.set_xlim(left = 0, right = 24)
    ax2.set_xticks(np.arange(0,25))
    ax2.set_xticklabels([])
    ax2.set_ylim(bottom = lowT, top = highT)
    ax2.set_yticks([272,274,276,278,280])
    ax2.set_title('GJOGR - Most southern, downstream coast')
    ax2.set_ylabel('K')
    
    ax3.plot(time_x[1::2], np.array(BOLUN_obs.TEMP) + CtK)
    ax3.plot(time_x, temp_subcube.data[:,3])
    ax3.set_xlim(left = 0, right = 24)
    ax3.set_xticks(np.arange(0,25))
    ax3.set_xticklabels(['00:00','','','','','','06:00','','','','','','12:00','','','','','','18:00','','','','','','00:00'])
    ax3.set_ylim(bottom = lowT, top = highT)
    ax3.set_yticks([272,274,276,278,280])
    ax3.set_title('BOLUN - Upstream coast on northern fjord')
    ax3.set_ylabel('K')
    
    fig.suptitle('2m Temperature - 19th March 2018 (case 306)')
    plt.tight_layout()
    plt.savefig('D:/Project/Figures/PNG/306/u-cc134/groundstations/RA1M_2m_temperature_4stations_306.png')

#%% plot wind
plot_wind = True
if plot_wind:
    fig, (ax0,ax1,ax2,ax3) = plt.subplots(4,1, figsize = (10,10))
    CtK = 273.15 # Celsius to Kelvin conversion
    ax0.plot(time_x[1::2], np.array(HBV_obs.WSP), label = 'obs')
    ax0.plot(time_x, wsp_subcube.data[:,0], label = 'UM 0p5km')
    ax0.set_xlim(left = 0, right = 24)
    ax0.set_xticks(np.arange(0,25))
    ax0.set_xticklabels([])
    ax0.set_ylim(bottom = 0, top = 15)
    ax0.set_yticks([0,4,8,12])
    ax0.set_title('HBV - Most northern, downstream coast')
    ax0.legend(fontsize = 12)
    ax0.set_ylabel(r'ms$^{-1}$')
    
    ax1.plot(time_x[1::2], np.array(AEDEY_obs.WSP))
    ax1.plot(time_x, wsp_subcube.data[:,1])
    ax1.set_xlim(left = 0, right = 24)
    ax1.set_xticks(np.arange(0,25))
    ax1.set_xticklabels([])
    ax1.set_ylim(bottom = 0, top = 15)
    ax1.set_yticks([0,4,8,12])
    ax1.set_title('AEDEY - Upstream coast on major fjord')
    ax1.set_ylabel(r'ms$^{-1}$')
    
    ax2.plot(time_x[1::2], np.array(GJOGR_obs.WSP))
    ax2.plot(time_x, wsp_subcube.data[:,2])
    ax2.set_xlim(left = 0, right = 24)
    ax2.set_xticks(np.arange(0,25))
    ax2.set_xticklabels([])
    ax2.set_ylim(bottom = 0, top = 15)
    ax2.set_yticks([0,4,8,12])
    ax2.set_title('GJOGR - Most southern, downstream coast')
    ax2.set_ylabel(r'ms$^{-1}$')
    
    ax3.plot(time_x[1::2], np.array(BOLUN_obs.WSP))
    ax3.plot(time_x, wsp_subcube.data[:,3])
    ax3.set_xlim(left = 0, right = 24)
    ax3.set_xticks(np.arange(0,25))
    ax3.set_xticklabels(['00:00','','','','','','06:00','','','','','','12:00','','','','','','18:00','','','','','','00:00'])
    ax3.set_ylim(bottom = 0, top = 15)
    ax3.set_yticks([0,4,8,12])
    ax3.set_title('BOLUN - Upstream coast on northern fjord')
    ax3.set_ylabel(r'ms$^{-1}$')
    
    fig.suptitle('10m windspeed - 19th March 2018 (case 306)')
    plt.tight_layout()
    plt.savefig('D:/Project/Figures/PNG/306/u-cc134/groundstations/RA1M_10m_windspeed_4stations_306.png')

#%% plot wind
plot_rh = True
if plot_rh:
    fig, (ax0,ax1,ax2,ax3) = plt.subplots(4,1, figsize = (10,10))
    CtK = 273.15 # Celsius to Kelvin conversion
    ax0.plot(time_x[1::2], np.array(HBV_obs.RH), label = 'obs')
    ax0.plot(time_x, rh_subcube.data[:,0], label = 'UM 0p5km')
    ax0.set_xlim(left = 0, right = 24)
    ax0.set_xticks(np.arange(0,25))
    ax0.set_xticklabels([])
    ax0.set_ylim(bottom = 0, top = 105)
    ax0.set_yticks([0,20,40,60,80,100])
    ax0.set_title('HBV - Most northern, downstream coast')
    ax0.legend(fontsize = 12)
    ax0.set_ylabel('%')
    
    ax1.plot(time_x[1::2], np.array(AEDEY_obs.RH))
    ax1.plot(time_x, rh_subcube.data[:,1])
    ax1.set_xlim(left = 0, right = 24)
    ax1.set_xticks(np.arange(0,25))
    ax1.set_xticklabels([])
    ax1.set_ylim(bottom = 0, top = 105)
    ax1.set_yticks([0,20,40,60,80,100])
    ax1.set_title('AEDEY - Upstream coast on major fjord')
    ax1.set_ylabel('%')
    
    ax2.plot(time_x[1::2], np.array(GJOGR_obs.RH))
    ax2.plot(time_x, rh_subcube.data[:,2])
    ax2.set_xlim(left = 0, right = 24)
    ax2.set_xticks(np.arange(0,25))
    ax2.set_xticklabels([])
    ax2.set_ylim(bottom = 0, top = 105)
    ax2.set_yticks([0,20,40,60,80,100])
    ax2.set_title('GJOGR - Most southern, downstream coast')
    ax2.set_ylabel('%')
    
    ax3.plot(time_x[1::2], np.array(BOLUN_obs.RH))
    ax3.plot(time_x, rh_subcube.data[:,3])
    ax3.set_xlim(left = 0, right = 24)
    ax3.set_xticks(np.arange(0,25))
    ax3.set_xticklabels(['00:00','','','','','','06:00','','','','','','12:00','','','','','','18:00','','','','','','00:00'])
    ax3.set_ylim(bottom = 0, top = 105)
    ax3.set_yticks([0,20,40,60,80,100])
    ax3.set_title('BOLUN - Upstream coast on northern fjord')
    ax3.set_ylabel('%')
    
    fig.suptitle('2m Relative Humidity - 19th March 2018 (case 306)')
    plt.tight_layout()
    plt.savefig('D:/Project/Figures/PNG/306/u-cc134/groundstations/RA1M_2m_relative_humidity_4stations_306.png')

