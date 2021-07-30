#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 10:05:32 2020

@author: kse18nru
"""

#%% module imports
import iris
from iris.analysis.cartography import rotate_pole
import data_management as dmna
import pandas as pd
import numpy as np
import xarray as xr
from datetime import datetime
import persinterpolate as pert
import matplotlib
import matplotlib.pyplot as plt
import xsection_workshop as workshop
import iris.coord_categorisation
import iris.plot as iplt
import matplotlib.gridspec as gridspec
from scipy.fft import fft, fftfreq, fftshift
from model_foundry import toolkit
#%%
plt.rcParams['font.size'] = 26
resolutions = ['0p5km', '1p5km', '4p4km']
varname_u = 'm01s15i002'
varname_v = 'm01s15i003'
varname_wsp = 'windspeed'
varname_w = 'upward_air_velocity'
varname_theta = 'air_potential_temperature'
varname_q = 'specific_humidity'

lat_lon_ord = ('grid_latitude', 'grid_longitude')
lon_lat_ord = ('grid_longitude', 'grid_latitude')
leg_meta_obs = ((1,26,1, lon_lat_ord),(2,27,2,lon_lat_ord),(3,28,3,lon_lat_ord),(4,29,3,lon_lat_ord),
            (5,30,3, lon_lat_ord),(6,30,3, lon_lat_ord),(8,31,8, lat_lon_ord),(9,32,8, lat_lon_ord))
leg_meta_time = ((1,26,1, lon_lat_ord),(2,27,2,lon_lat_ord),(3,28,3,lon_lat_ord),(8,31,8, lat_lon_ord))
flight = 301
config = 'RA1M'
exp = 'Control'
suite = 'u-bu807'
res = resolutions[2]
#%%
fig = plt.figure(figsize = (16,24))
gs = gridspec.GridSpec(nrows = 8, ncols = 1)
fig.suptitle(f'{res} {varname_w} spectra on flight {flight}')
#fig.suptitle(f'{config}_{res}_{varname_w}_fft_flt{flight}.png')
for j in range(1):
    fileleg =leg_meta_obs[j]
    
    
    #%%
    
    w_cube = dmna.load_model_w(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
    
    
    #%% load and process databse shortrun obs data
    # using pandas
    database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_sr_301.txt', delimiter = ' ')
    
    db_legno = database['legno']
    #legcond = np.array(db_legno == 3) + np.array(db_legno == 4) + np.array(db_legno == 5) + np.array(db_legno == 6)
    legcond = db_legno == fileleg[0]
    db_wsp = np.array(database['wsp'][legcond]) 
    db_w = np.array(database['w'][legcond])
    db_q = np.array(database['q_BUCK'][legcond])
    db_theta = np.array(database['theta'][legcond])
    
    db_lon = np.array(database['lon'][legcond])
    db_lat = np.array(database['lat'][legcond])
    db_alt = np.array(database['altgps'][legcond])
    #%%
    db_meantime = np.array(database['meantime'][legcond])
    print('start-time:',db_meantime[0])
    time_format = '%Y-%m-%dT%H:%M:%S'
    midnight = datetime.timestamp(datetime(year = 2018, month = 3, day = 12, hour = 0, minute = 0, second = 0))
    start_time = datetime.fromtimestamp(db_meantime[0]).strftime(time_format)
    print(start_time)
    #stop_time = db_meantime[-1] + midnight
    
    #%% load qc MASIN 1Hz data
    masin_file = xr.open_dataset('../../Obvs_Data/UEA_qc_MASIN_data/MASIN_flightdata_301.nc')
    print(masin_file)
    #%%
    #obs_time_1hz = np.array(masin_file[''])
    temp_time = np.array(masin_file['posixtime_50hz'])[::50]
    print(temp_time)
    #%%
    time50hz = np.datetime64(temp_time) - midnight
    alt50hz =np.array(masin_file['gps_alt_50hz'])[::50]
    w50hz = np.array(masin_file['w_50hz'])[::50]
    # for t in range(len(temp_time)):
    #     time50hz.append(datetime.timestamp(temp_time[t])-midnight)
    # time50hz = np.array(time50hz)
    print(time50hz)
    #%% rotate db coords
    compcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/RA1M_0p5km_umh1_flt{flight}.nc', varname_u)[0,0]
    polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
    rot_db_lon = rot_db_lon + 360
    
    #%% interpolate model data onto obs altitudes
    time_index = fileleg[1]
    if fileleg[3][0] == 'grid_longitude':
        
        
        matched_w, matched_coors = pert.interpolate_matched_series(w_cube, db_alt, rot_db_lon, time_index)
        simple_w, simple_coors_w = pert.interpolate_simple_series(w_cube, np.nanmean(db_alt), time_index)
        
        
        db_coor = rot_db_lon
    elif fileleg[3][0] == 'grid_latitude':
        
        
        matched_w, matched_coors = pert.interpolate_matched_series(w_cube, db_alt, rot_db_lat, time_index, grid_axis = 'grid_latitude')
        simple_w, simple_coors_w = pert.interpolate_simple_series(w_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
        
        
        db_coor = rot_db_lat
    
    #%% get rough spatial spacing of data points
    ptp_dists = []
    for i in range(len(db_coor)-1):
        step = toolkit.haversine(rot_db_lat[i], rot_db_lon[i], rot_db_lat[i+1], rot_db_lon[i+1]) # returns distance on globe in km
        ptp_dists.append(step)
    ptp_dists = np.array(ptp_dists)
    ptp_mean = np.mean(ptp_dists)
    ptp_std   = np.std(ptp_dists)
    print(f'PTP distance = {ptp_mean} +/- {ptp_std}')
    d = ptp_mean
    #%% dynamically adjust horizontal bounds
    xmin = np.nanmin(simple_coors_w)
    xmax = np.nanmax(simple_coors_w)
    
    #%% perform FFT on matched model dataset
    model_ft = fftshift(fft(matched_w))
    model_freq = fftshift(fftfreq(len(matched_w), d))
    wavelengths = 1/model_freq
    #%% perform FFT on long model dataset
    long_ft = fftshift(fft(simple_w))
    long_freq = fftshift(fftfreq(len(simple_w), d))
    
    #%% perform FFT on obs dataset
    obs_ft = fftshift(fft(db_w))
    obs_freq = fftshift(fftfreq(len(db_w), d))
    #%%
    N = 400
    T = 1/3600
    line_sample = np.linspace(0, N*T, N)
    sine_sample = np.sin(line_sample) 
    sine_fft = fftshift(fft(sine_sample))
    fft_freq = fftshift(fftfreq(N,T))
    #%% plot spectra
    
    ax0 = fig.add_subplot(gs[j,0])
    #ax.plot(fft_freq, np.abs(sine_fft))
    plt.plot(1/model_freq, np.abs(model_ft), label = f'UM leg {fileleg[0]}')#/np.sum(np.abs(model_ft)))
    plt.plot(1/obs_freq, np.abs(obs_ft), label = 'obs')#/np.sum(np.abs(obs_ft)))
    plt.legend()
    #ax.plot(long_freq, np.abs(long_ft))#/np.sum(np.abs(long_ft)))
    plt.xlim(left = 0, right = 80)
    plt.ylim(top = 14)
    #ax.set_xlabel('wavelength [km]')
    #fig.suptitle(f'Spectral analysis of vertical windspeed over leg {fileleg[0]}')
    
    
#%%
ax0.set_xlabel('Wavelength [km]')
fig.tight_layout()
savefile = False
if savefile:
    fig.savefig(f'../../Figures/PNG/{flight}/{suite}/fft/{exp}/{varname_w}/{res}/{config}_{res}_{varname_w}_fft_flt{flight}.png')
