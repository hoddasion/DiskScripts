# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 13:48:24 2021

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
obs_path = 'D:/Project/Obvs_Data/groundstations/snaesfellsnes_obs_201803.csv'
stations_path = 'D:/Project/Obvs_Data/groundstations/snaesfellsnes_stations_info.csv'

obs_df = pd.read_csv(obs_path)
obs_df['TIMESTAMP'] = pd.to_datetime((obs_df['TIMESTAMP']- obs_df.TIMESTAMP[0])*24, unit = 'h', origin =pd.Timestamp('2018-03-01'))

gsta_df = pd.read_csv(stations_path) # gsta: groundstations
gsta_lon = np.array(gsta_df.LON)
gsta_lat = np.array(gsta_df.LAT)
gsta_sname = np.array(gsta_df.SHORTNAME)
print(gsta_df)
print(obs_df[obs_df.STATION == 2050])

#%%
lower_time = obs_df[obs_df.TIMESTAMP <= datetime.datetime(2018,3,13,0,0)]
obs_12df = lower_time[lower_time.TIMESTAMP > datetime.datetime(2018,3,12,0,0)]
print(obs_12df.TIMESTAMP)

#%% subset obs data into stations 
BLAFE_cond = np.array(obs_12df.STATION) == 1936
BLAFE_obs = obs_12df[BLAFE_cond]
STH_cond = np.array(obs_12df.STATION) == 2050
STH_obs = obs_12df[STH_cond]
BULAH_cond = np.array(obs_12df.STATION) == 31932
BULAH_obs = obs_12df[BULAH_cond]
GUFUS_cond = np.array(obs_12df.STATION) == 1919
GUFUS_obs = obs_12df[GUFUS_cond]
print(obs_12df[::12])
#%% testing of time series
test = obs_12df.TEMP[obs_12df.STATION == 31932] + 273.15
test.plot()
#%% model data
um_path = 'D:/Project/Model_Data/u-bu807/nc/Control/'
p_cube = iris.load_cube(f'{um_path}RA1M_0p5km_um_surface_air_pressure_24hrs_g_301.nc')
rh_cube = iris.load_cube(f'{um_path}RA1M_0p5km_um_relative_humidity_24hrs_g_301.nc')
land_cube = iris.load_cube(f'{um_path}RA1M_0p5km_umg1_flt301.nc', 'land_binary_mask')
temp_cube = iris.load_cube(f'{um_path}RA1M_0p5km_um_air_temperature_24hrs_g_301.nc')
xwsp_cube = iris.load_cube(f'{um_path}RA1M_0p5km_um_x_wind_24hrs_g_301.nc')
ywsp_cube = iris.load_cube(f'{um_path}RA1M_0p5km_um_y_wind_24hrs_g_301.nc')
wsp_cube = iris.load_cube(f'{um_path}RA1M_0p5km_um_x_wind_24hrs_g_301.nc')
wsp_data = (xwsp_cube.data**2 + ywsp_cube.data[:,:-1]**2)**0.5
wsp_cube.data = wsp_data

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
model_x.append(gsta_x[gsta_sname == 'GUFUS'][0]); model_y.append(gsta_y[gsta_sname == 'GUFUS'][0])
model_x.append(gsta_x[gsta_sname == 'BLAFE'][0]); model_y.append(gsta_y[gsta_sname == 'BLAFE'][0])
model_x.append(gsta_x[gsta_sname == 'STH'][0]); model_y.append(gsta_y[gsta_sname == 'STH'][0])
model_x.append(gsta_x[gsta_sname == 'BULAH'][0]); model_y.append(gsta_y[gsta_sname == 'BULAH'][0])
sample_points =  [('grid_longitude',np.array(model_x)), ('grid_latitude', np.array(model_y))]
print(model_x, model_y)
temp_subcube = trajectory.interpolate(temp_cube,sample_points, method =  'nearest') # temp meaning temperature, not temporary
wsp_subcube = trajectory.interpolate(wsp_cube,sample_points, method =  'nearest')
rh_subcube = trajectory.interpolate(rh_cube,sample_points, method =  'nearest')
p_subcube = trajectory.interpolate(p_cube, sample_points, method = 'nearest')
#%% measure lapse rate in coloumns
Tht3D_cube = iris.load_cube(f'{um_path}RA1M_0p5km_um_air_potential_temperature_24hrs_i_301.nc', 'air_potential_temperature')
#Pres3D_cube = iris.load_cube(f'{um_path}RA1M_0p5km_um_air_pressure_24hrs_i_301.nc', 'air_pressure')
#%%
import surface_evaluation_functions as sef
#%%
#thtrate_df = sef.measure_theta_rate(Tht3D_cube, sample_points,1000)
#print(thtrate_df[::12])
#presrate_df = sef.measure_pressure_rate(Pres3D_cube, sample_points, 1000)
#(presrate_df)

#%% calculate lapse rate
temp_subdata = temp_subcube.data ## get data for array reformatting
T_1Darray = np.concatenate((temp_subdata[:,0],temp_subdata[:,1],temp_subdata[:,2],temp_subdata[:,3]))

print(temp_subcube)
#lapse_df = sef.lapse_rate(thtrate_df, presrate_df, T = T_1Darray)
#print(lapse_df[48:72])

#%% find difference in altitudes and calculate teperature corrections
calc_tcorrection = False
if calc_tcorrection:
    Tht3D_subcube  = trajectory.interpolate(Tht3D_cube,sample_points, method =  'nearest') 
    model_z = Tht3D_subcube.coord('surface_altitude').points
    print(Tht3D_subcube.coord('surface_altitude'))
    print(gsta_df)
    alt_BLAFE = 17.7; alt_STH = 12.4; alt_BULAH = 25.0; alt_GUFUS = 7.0
    ## recall index 0->GUFUS, 1->BLAFE, 2->STH, 3->BULAH
    ## dz = dz(um) - dz(obs)
    dz_GUFUS = model_z[0] - alt_GUFUS;
    dz_BLAFE = model_z[1] - alt_BLAFE; 
    dz_STH = model_z[2] - alt_STH; 
    dz_BULAH = model_z[3] - alt_BULAH; 
    padia_rate = 6/1000 # pseudo-adiabatic lapse rate ~ 6K/km
    dT_GUFUS = dz_GUFUS*padia_rate; 
    dT_BLAFE = dz_BLAFE*padia_rate
    dT_STH = dz_STH*padia_rate
    dT_BULAH = dz_BULAH*padia_rate
    print(dz_GUFUS, dT_GUFUS)
    print(dz_BLAFE, dT_BLAFE)
    print(dz_STH, dT_STH)
    print(dz_BULAH, dT_BULAH)

#
#%% make time (x-axis) list
time_x = np.arange(1,49)*0.5
print(time_x[1::2])

#%% temp mean analysis
Tmean_analysis = False
if Tmean_analysis:
    fig0, ax0 = plt.subplots(1,1, figsize = (10,4))
    ax0.plot(time_x, temp_subdata[:,0], label = 'GUFUS')
    ax0.plot(time_x, temp_subdata[:,1], label = 'BLAFE')
    ax0.plot(time_x, temp_subdata[:,2], label = 'STH')
    ax0.plot(time_x, temp_subdata[:,3], label = 'BULAH')
    ax0.plot(time_x, np.nanmean(temp_subdata,axis = 1), label = 'mean')
    ax0.legend(fontsize = 12)

#%% plot temperature
plot_temp = False
if plot_temp:
    lowT = 267.5; highT = 276.5; Tticks = [268,270,272,274,276]
    fig, (ax0,ax1,ax2,ax3) = plt.subplots(4,1, figsize = (10,10))
    CtK = 273.15 # Celsius to Kelvin conversion
    ax0.plot(time_x[1::2], np.array(GUFUS_obs.TEMP) + CtK, label = 'obs')
    ax0.plot(time_x, temp_subcube.data[:,0] -0.006, label = 'UM 0p5km')
    ax0.set_xlim(left = 0, right = 24)
    ax0.set_xticks(np.arange(0,25))
    ax0.set_xticklabels([])
    ax0.set_ylim(bottom = lowT, top = highT)
    ax0.set_yticks(Tticks)
    ax0.set_title('GUFUS (1919)')
    ax0.legend(fontsize = 12)
    ax0.set_ylabel('K')
    
    ax1.plot(time_x[1::2], np.array(BLAFE_obs.TEMP) + CtK)
    ax1.plot(time_x, temp_subcube.data[:,1] - 0.063)
    ax1.set_xlim(left = 0, right = 24)
    ax1.set_xticks(np.arange(0,25))
    ax1.set_xticklabels([])
    ax1.set_ylim(bottom = lowT, top = highT)
    ax1.set_yticks(Tticks)
    ax1.set_title('BLAFE (1936)')
    ax1.set_ylabel('K')
    
    ax2.plot(time_x[1::2], np.array(STH_obs.TEMP) + CtK)
    ax2.plot(time_x, temp_subcube.data[:,2] - 0.074)
    ax2.set_xlim(left = 0, right = 24)
    ax2.set_xticks(np.arange(0,25))
    ax2.set_xticklabels([])
    ax2.set_ylim(bottom = lowT, top = highT)
    ax2.set_yticks(Tticks)
    ax2.set_title('STH (2050)')
    ax2.set_ylabel('K')
    
    ax3.plot(time_x[1::2], np.array(BULAH_obs.TEMP) + CtK)
    ax3.plot(time_x, temp_subcube.data[:,3] - 0.157)
    ax3.set_xlim(left = 0, right = 24)
    ax3.set_xticks(np.arange(0,25))
    ax3.set_xticklabels(['00:00','','','','','','06:00','','','','','','12:00','','','','','','18:00','','','','','','00:00'])
    ax3.set_ylim(bottom = lowT, top = highT)
    ax3.set_yticks(Tticks)
    ax3.set_title('BULAH (31932)')
    ax3.set_ylabel('K')
    
    fig.suptitle('2m Temperature - 12th March 2018 (case 301)')
    plt.tight_layout()
    plt.savefig('D:/Project/Figures/PNG/301/u-bu807/groundstations/RA1M_2m_temperature_snaes4stations_301.png')

#%%
plot_wsp = True
if plot_wsp:
    lowT = -0.5; highT = 15.5; Tticks = [0,5,10,15]
    fig, (ax0,ax1,ax2,ax3) = plt.subplots(4,1, figsize = (10,10))
    CtK = 273.15 # Celsius to Kelvin conversion
    ax0.plot(time_x[1::2], np.array(GUFUS_obs.F) , label = 'obs')
    ax0.plot(time_x, wsp_subcube.data[:,0], label = 'UM 0p5km')
    ax0.set_xlim(left = 0, right = 24)
    ax0.set_xticks(np.arange(0,25))
    ax0.set_xticklabels([])
    ax0.set_ylim(bottom = lowT, top = highT)
    ax0.set_yticks(Tticks)
    ax0.set_title('GUFUS (1919)')
    ax0.legend(fontsize = 12)
    ax0.set_ylabel(r'ms$^{-1}$')
    
    ax1.plot(time_x[1::2], np.array(BLAFE_obs.F) )
    ax1.plot(time_x, wsp_subcube.data[:,1])
    ax1.set_xlim(left = 0, right = 24)
    ax1.set_xticks(np.arange(0,25))
    ax1.set_xticklabels([])
    ax1.set_ylim(bottom = lowT, top = highT)
    ax1.set_yticks(Tticks)
    ax1.set_title('BLAFE (1936)')
    ax1.set_ylabel(r'ms$^{-1}$')
    
    ax2.plot(time_x[1::2], np.array(STH_obs.F) )
    ax2.plot(time_x, wsp_subcube.data[:,2])
    ax2.set_xlim(left = 0, right = 24)
    ax2.set_xticks(np.arange(0,25))
    ax2.set_xticklabels([])
    ax2.set_ylim(bottom = lowT, top = highT)
    ax2.set_yticks(Tticks)
    ax2.set_title('STH (2050)')
    ax2.set_ylabel(r'ms$^{-1}$')
    
    ax3.plot(time_x[1::2], np.array(BULAH_obs.F))
    ax3.plot(time_x, wsp_subcube.data[:,3])
    ax3.set_xlim(left = 0, right = 24)
    ax3.set_xticks(np.arange(0,25))
    ax3.set_xticklabels(['00:00','','','','','','06:00','','','','','','12:00','','','','','','18:00','','','','','','00:00'])
    ax3.set_ylim(bottom = lowT, top = highT)
    ax3.set_yticks(Tticks)
    ax3.set_title('BULAH (31932)')
    ax3.set_ylabel(r'ms$^{-1}$')
    
    fig.suptitle('10m Windspeed - 12th March 2018 (case 301)')
    plt.tight_layout()
    plt.savefig('D:/Project/Figures/PNG/301/u-bu807/groundstations/RA1M_10m_windspeed_snaes4stations_301.png')
    
#%%
plot_rh = False
if plot_rh:
    lowT = -0.5; highT = 101; Tticks = [0,20,40,60,80,100]
    fig, (ax0,ax1,ax2,ax3) = plt.subplots(4,1, figsize = (10,10))
    CtK = 273.15 # Celsius to Kelvin conversion
    ax0.plot(time_x[1::2], np.array(GUFUS_obs.RH) , label = 'obs')
    ax0.plot(time_x, rh_subcube.data[:,0], label = 'UM 0p5km')
    ax0.set_xlim(left = 0, right = 24)
    ax0.set_xticks(np.arange(0,25))
    ax0.set_xticklabels([])
    ax0.set_ylim(bottom = lowT, top = highT)
    ax0.set_yticks(Tticks)
    ax0.set_title('GUFUS (1919)')
    ax0.legend(fontsize = 12)
    ax0.set_ylabel('%')
    
    ax1.plot(time_x[1::2], np.array(BLAFE_obs.RH) )
    ax1.plot(time_x, rh_subcube.data[:,1])
    ax1.set_xlim(left = 0, right = 24)
    ax1.set_xticks(np.arange(0,25))
    ax1.set_xticklabels([])
    ax1.set_ylim(bottom = lowT, top = highT)
    ax1.set_yticks(Tticks)
    ax1.set_title('BLAFE (1936)')
    ax1.set_ylabel('%')
    
    ax2.plot(time_x[1::2], np.array(STH_obs.RH) )
    ax2.plot(time_x, rh_subcube.data[:,2])
    ax2.set_xlim(left = 0, right = 24)
    ax2.set_xticks(np.arange(0,25))
    ax2.set_xticklabels([])
    ax2.set_ylim(bottom = lowT, top = highT)
    ax2.set_yticks(Tticks)
    ax2.set_title('STH (2050)')
    ax2.set_ylabel('%')
    
    ax3.plot(time_x[1::2], np.array(BULAH_obs.RH))
    ax3.plot(time_x, rh_subcube.data[:,3])
    ax3.set_xlim(left = 0, right = 24)
    ax3.set_xticks(np.arange(0,25))
    ax3.set_xticklabels(['00:00','','','','','','06:00','','','','','','12:00','','','','','','18:00','','','','','','00:00'])
    ax3.set_ylim(bottom = lowT, top = highT)
    ax3.set_yticks(Tticks)
    ax3.set_title('BULAH (31932)')
    ax3.set_ylabel('%')
    
    fig.suptitle('2m Relative Humidity - 12th March 2018 (case 301)')
    plt.tight_layout()
    plt.savefig('D:/Project/Figures/PNG/301/u-bu807/groundstations/RA1M_2m_relative_humidity_snaes4stations_301.png')
    
#%%
plot_p = False
if plot_p:
    
    print(np.array(GUFUS_obs.PRES).astype(np.float))
    print(type(np.array(GUFUS_obs.PRES)[12]))
    #%%
    
    lowT = 1004; highT = 1013; Tticks = [1005,1010]
    fig, (ax0,ax3) = plt.subplots(2,1, figsize = (10,5))
    CtK = 273.15 # Celsius to Kelvin conversion
    ax0.plot(time_x[1::2], np.array(GUFUS_obs.PRES).astype(np.float), label = 'obs')
    ax0.plot(time_x, p_subcube.data[:,0]/100 -0.126, label = 'UM 0p5km')
    ax0.set_xlim(left = 0, right = 24)
    ax0.set_xticks(np.arange(0,25))
    ax0.set_xticklabels([])
    ax0.set_ylim(bottom = lowT, top = highT)
    ax0.set_yticks(Tticks)
    ax0.set_title('GUFUS (1919)')
    ax0.legend(fontsize = 12)
    ax0.set_ylabel('hPa')
    
    
    
    ax3.plot(time_x[1::2], np.array(STH_obs.PRES).astype(np.float))
    ax3.plot(time_x, p_subcube.data[:,2]/100 - 1.56)
    ax3.set_xlim(left = 0, right = 24)
    ax3.set_xticks(np.arange(0,25))
    ax3.set_xticklabels(['00:00','','','','','','06:00','','','','','','12:00','','','','','','18:00','','','','','','00:00'])
    ax3.set_ylim(bottom = lowT, top = highT)
    ax3.set_yticks(Tticks)
    ax3.set_title('STH (2050)')
    ax3.set_ylabel('hPa')
    
    fig.suptitle('2m Pressure (corr) - 12th March 2018 (case 301)')
    plt.tight_layout()
    plt.savefig('D:/Project/Figures/PNG/301/u-bu807/groundstations/RA1M_2m_pressure_corrected_snaes2stations_301.png')
