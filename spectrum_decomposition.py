# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 13:19:24 2022

@author: kse18nru
"""

#%% module imports


import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, fftshift
from model_foundry import toolkit
import filtering
import iris
import matplotlib.gridspec as gridspec
import iris.coord_categorisation
import iris.plot as iplt

#%% globals
suite = 'u-cc134'
experiment = 'CONTROL'
config = 'RA1M'
res = '0p5km'
flight = 306
matplotlib.rcParams.update({'font.size': 24})
#%% do obs analysis
obs_analysis= False
if obs_analysis:
    if flight == 306:
        obs_filename = f'IGP_flights_database_obs_sr_{flight}.txt'
        path = 'D:/Project/Obvs_Data/Databases/'
        
        df_obs = pd.read_csv(f'{path}{obs_filename}', delimiter = ' ')
        print(df_obs.columns)
        
        df_model = pd.read_csv(f'D:/Project/Model_Data/{suite}/Matched_interpolated_values_0p5km_sr_306.txt')
        
        print(df_model)
        
        #%% extract variables from obs and model
        
        ## set condition for point selection base don leg
        leg = 5
        condition = np.array(df_obs.legno) == leg
        
        ## extract coordinates
        obs_lon = np.array(df_obs.lon)[condition]
        obs_lat = np.array(df_obs.lat)[condition]
        
        ## extract obs[condition]
        obs_wsp = np.array(df_obs.wsp)[condition]
        obs_theta = np.array(df_obs.theta)[condition]
        obs_q = np.array(df_obs.q_BUCK)[condition]
        obs_w = np.array(df_obs.w)[condition]
        
        ## extract model vars
        um_wsp = np.array(df_model.wsp)[condition]
        um_theta = np.array(df_model.theta)[condition]
        um_q = np.array(df_model.q)[condition]
        um_w = np.array(df_model.w)[condition]
        
        #%% calcualte point to point distance
        ptp_dists = []
        cum_dists = [0]
        for i in range(len(obs_lon)-1):
            step = toolkit.haversine(obs_lat[i], obs_lon[i], obs_lat[i+1], obs_lon[i+1]) # returns distance on globe in km
            ptp_dists.append(step)
            cum_dists.append(cum_dists[i] + ptp_dists[i])
        ptp_dists = np.array(ptp_dists)
        #%% plot time series
        
        fig, (ax0,ax1,ax2,ax3) = plt.subplots(4,1, figsize = (16,16))
        
        ax0.plot(cum_dists, obs_w)
        ax0.plot(cum_dists, um_w)
        
        ax1.plot(cum_dists, um_theta)

#%% do model only analysis  

w_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/upward_air_velocity/{config}_vertslice_{res}_upward_air_velocity_flt{flight}_leg7.nc', 'upward_air_velocity')[29,:25]
hidx = 14
print(w_cube[hidx])
sample_alt = w_cube[hidx].coords('altitude')[0].points[0]
leftidx = 25; rightidx = 70
gridlat = w_cube.coords('grid_latitude')[0].points[leftidx:rightidx]
gridlon = w_cube.coords('grid_longitude')[0].points[leftidx:rightidx]
sample_xbounds = [gridlat[0],gridlat[-1]]
sample_ybounds = [sample_alt, sample_alt]

#%% calcualte point to point distance
ptp_dists = []
cum_dists = [0]
for i in range(len(gridlon)-1):
    step = toolkit.haversine(gridlat[i], gridlon[i], gridlat[i+1], gridlon[i+1]) # returns distance on globe in km
    ptp_dists.append(step)
    cum_dists.append(cum_dists[i] + ptp_dists[i])
ptp_dists = np.array(ptp_dists)
cum_dists = np.array(cum_dists)
print(np.median(ptp_dists))

 #%%
    
length_ticks = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
                1,2,3,4,5,6,7,8,9,
                10,20,30,40,50,60,70,80,90,
                100,200,300,400,500,600,700,800,900,
                1000]
log_ticks = np.log10(length_ticks)
length_ticks_labels = [0.1,'','','','','','','','',
                       1,'','','','','','','','',
                       10,'','','','','','','','',
                       100,'','','','','','','','',
                       1000]

#%% fourier transform 1
w_data = np.array(w_cube[hidx].data)[leftidx:rightidx]
k_num = fftfreq(len(w_data), np.median(ptp_dists)) # due to spatial series, this is the wavenumber
w_fft = fftshift(fft(w_data))

#%% smooth and reperform fft
smooth_window = 5
w_data_smooth = filtering.moving_average(w_data, smooth_window)
w_fft_smooth = fftshift(fft(w_data_smooth))
k_num_smooth = fftfreq(len(w_data_smooth), np.median(ptp_dists))
#%% normalise colorscales
w_norm = matplotlib.colors.Normalize(vmin = -5, vmax = 5)
w_map = matplotlib.cm.ScalarMappable(norm = w_norm, cmap = 'seismic')
#%% plot
geo_axis = 'grid_latitude'
common_kwargs = {'coords' : [geo_axis,'altitude']}
com_cbarkeys = {'pad' : 0.005}
fig = plt.figure(figsize = (16,16))
gs = gridspec.GridSpec(nrows = 4, ncols = 1)
top_alt = 2000

## cross-section panel
ax0 = fig.add_subplot(gs[0,0])
ax0.set_ylim(top = top_alt)
ax0.set_xlim(left = 0.2, right = 1.0)
p_xsec = iplt.contourf(w_cube, **common_kwargs, norm = w_norm, cmap = 'seismic', levels = 40)
ax0.plot(sample_xbounds, sample_ybounds, color = 'k', linestyle = 'dashed')

## spatial series
ax1 = fig.add_subplot(gs[1,0])
ax1.plot(w_cube.coords(geo_axis)[0].points, w_cube[hidx].data)
ax1.plot(gridlat[int((smooth_window - 1)/2): - int((smooth_window - 1)/2)], w_data_smooth)
ax1.set_xlim(left = 0.2, right = 1.0)

## raw fft
ax2 = fig.add_subplot(gs[2,0])
ax2.plot(k_num[k_num >= 0], w_fft[k_num >= 0])
ax2.plot(k_num_smooth[k_num_smooth >= 0], w_fft_smooth[k_num_smooth >= 0])