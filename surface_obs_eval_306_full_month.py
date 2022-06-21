# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 15:24:07 2021

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
#%%
plt.rcParams.update({'font.size': 16})
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
#%%
obs_df.TIMESTAMP = obs_df.TIMESTAMP.round('h')
print(obs_df)

#%% subset obs data into stations 
HBV_cond = np.array(obs_df.STATION) == 2862
HBV_obs = obs_df[HBV_cond]
AEDEY_cond = np.array(obs_df.STATION) == 2655
AEDEY_obs = obs_df[AEDEY_cond]
GJOGR_cond = np.array(obs_df.STATION) == 2692
GJOGR_obs = obs_df[GJOGR_cond]
BOLUN_cond = np.array(obs_df.STATION) == 2738
BOLUN_obs = obs_df[BOLUN_cond]

#%%
plot_temp_tseries=False
if plot_temp_tseries:
    fig, (ax0,ax1,ax2,ax3) = plt.subplots(4,1, figsize = (16,14), sharex = True)
    CtK = 273.15
    window = 2
    ax0.plot(np.array(HBV_obs.TIMESTAMP), np.array(HBV_obs.TEMP) + CtK, color = 'tab:orange')
    ax0.set_ylim(bottom = 267, top = 283)
    ax0.set_title('HBV')
    ax0.set_ylabel('K')
    ax0b = ax0.twinx()
    ax0b.plot(np.array(HBV_obs.TIMESTAMP), np.array(HBV_obs.RH.rolling(window).mean()), color = 'tab:blue')
    ax0b.set_ylim(bottom = 38, top = 102)
    ax0b.set_ylabel('%')
    
    ax1.plot(np.array(AEDEY_obs.TIMESTAMP), np.array(AEDEY_obs.TEMP) + CtK, color = 'tab:orange')
    ax1.set_ylim(bottom = 267, top = 283)
    ax1.set_title('AEDEY')
    ax1.set_ylabel('K')
    ax1b = ax1.twinx()
    ax1b.plot(np.array(AEDEY_obs.TIMESTAMP), np.array(AEDEY_obs.RH.rolling(window).mean()), color = 'tab:blue')
    ax1b.set_ylim(bottom = 38, top = 102)
    ax1b.set_ylabel('%')
    
    ax2.plot(np.array(GJOGR_obs.TIMESTAMP), np.array(GJOGR_obs.TEMP) + CtK, color = 'tab:orange')
    ax2.set_ylim(bottom = 267, top = 283)
    ax2.set_title('GJOGR')
    ax2.set_ylabel('K')
    ax2b = ax2.twinx()
    ax2b.plot(np.array(GJOGR_obs.TIMESTAMP), np.array(GJOGR_obs.RH.rolling(window).mean()), color = 'tab:blue')
    ax2b.set_ylim(bottom = 38, top = 102)
    ax2b.set_ylabel('%')
    
    ax3.plot(np.array(BOLUN_obs.TIMESTAMP), np.array(BOLUN_obs.TEMP) + CtK, color = 'tab:orange')
    ax3.set_ylim(bottom = 267, top = 283)
    ax3.set_title('BOLUN')
    ax3.set_ylabel('K')
    
    ax3b = ax3.twinx()
    ax3b.plot(np.array(BOLUN_obs.TIMESTAMP), np.array(BOLUN_obs.RH.rolling(window).mean()), color = 'tab:blue')
    ax3b.set_ylim(bottom = 38, top = 102)
    ax3b.set_ylabel('%')
    
    ax3.set_xticks(np.array(BOLUN_obs.TIMESTAMP)[0::120])
    
    fig.suptitle('March Air Temperature (orange) and Relative Humidity (blue) Observed Ground Measurements')
    plt.tight_layout()
    
    plt.savefig('D:/Project/Figures/PNG/306/gsta_fullmarch_temp_rh.png')
    
    #%%
plot_wind_tseries=False
if plot_wind_tseries:
    fig, (ax1,ax2,ax3) = plt.subplots(3,1, figsize = (16,13), sharex = True)
    
    window = 6
    
    ax1.plot(np.array(AEDEY_obs.TIMESTAMP), np.array(AEDEY_obs.WSP.rolling(window).mean()), color = 'blue')
    ax1.set_ylim(bottom = 0, top = 24)
    ax1.set_title('AEDEY')
    ax1.set_ylabel(r'ms$^{-1}$')
    ax1b = ax1.twinx() # instantiate a second axes that shares the same x-axis
    ax1b.plot(np.array(AEDEY_obs.TIMESTAMP), np.array(AEDEY_obs.WDIR.rolling(window).mean()), color = 'k', linestyle = 'dashed')
    ax1b.set_ylim(bottom = 0, top = 360)
    ax1b.set_yticks([0,90,180,270,360])
    ax1b.set_yticklabels([r'$0^{\circ}$',r'$90^{\circ}$',r'$180^{\circ}$',r'$270^{\circ}$',r'$360^{\circ}$'])
    
    
    ax2.plot(np.array(GJOGR_obs.TIMESTAMP), np.array(GJOGR_obs.WSP.rolling(window).mean()), color = 'blue')
    ax2.set_ylim(bottom = 0, top = 24)
    ax2.set_title('GJOGR')
    ax2.set_ylabel(r'ms$^{-1}$')
    ax2b = ax2.twinx() # instantiate a second axes that shares the same x-axis
    ax2b.plot(np.array(GJOGR_obs.TIMESTAMP), np.array(GJOGR_obs.WDIR.rolling(window).mean()), color = 'k', linestyle = 'dashed')
    ax2b.set_ylim(bottom = 0, top = 360)
    ax2b.set_yticks([0,90,180,270,360])
    ax2b.set_yticklabels([r'$0^{\circ}$',r'$90^{\circ}$',r'$180^{\circ}$',r'$270^{\circ}$',r'$360^{\circ}$'])
    
    
    ax3.plot(np.array(BOLUN_obs.TIMESTAMP), np.array(BOLUN_obs.WSP.rolling(window).mean()),  color = 'blue')
    ax3.set_ylim(bottom = 0, top = 24)
    ax3.set_title('BOLUN')
    ax3.set_ylabel(r'ms$^{-1}$')
    
    ax3b = ax3.twinx() # instantiate a second axes that shares the same x-axis
    ax3b.plot(np.array(BOLUN_obs.TIMESTAMP), np.array(BOLUN_obs.WDIR.rolling(window).mean()), color = 'k', linestyle = 'dashed')
    ax3b.set_ylim(bottom = 0, top = 360)
    ax3b.set_yticks([0,90,180,270,360])
    ax3b.set_yticklabels([r'$0^{\circ}$',r'$90^{\circ}$',r'$180^{\circ}$',r'$270^{\circ}$',r'$360^{\circ}$'])
    
    ax3.set_xticks(np.array(BOLUN_obs.TIMESTAMP)[0::120])
    
    fig.suptitle('March Windspeed (blue) and Direction (dashed) Observed Ground Measurements')
    plt.tight_layout()
    
    plt.savefig('D:/Project/Figures/PNG/306/gsta_fullmarch_wsp_wdir.png')
#%%   
plt.rcParams.update({'font.size': 20})
plot_scatters_raw=True
if plot_scatters_raw:
    fig,((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2, figsize = (15,15))
    
    ax0.set_ylabel(r'wsp, ms$^{-1}$')
    ax0.set_xlabel('wdir')
    ax0.set_xticks([0,90,180,270,360])
    ax0.scatter(np.array(AEDEY_obs.WDIR),np.array(AEDEY_obs.WSP), label = 'AEDEY')
    ax0.scatter(np.array(GJOGR_obs.WDIR),np.array(GJOGR_obs.WSP), label = 'GJOGR')
    ax0.scatter(np.array(BOLUN_obs.WDIR),np.array(BOLUN_obs.WSP), label = 'BOLUN')
    ax0.legend()
    
    ax1.set_ylabel(r'air temp., $^{\circ}$C')
    ax1.set_xlabel('wdir')
    ax1.set_xticks([0,90,180,270,360])
    ax1.scatter(np.array(AEDEY_obs.WDIR),np.array(AEDEY_obs.TEMP))
    ax1.scatter(np.array(GJOGR_obs.WDIR),np.array(GJOGR_obs.TEMP))
    ax1.scatter(np.array(BOLUN_obs.WDIR),np.array(BOLUN_obs.TEMP))
    
    ax2.set_ylabel(r'air temp., $^{\circ}$C')
    ax2.set_xlabel('wsp, ms$^{-1}$')
    #ax2.set_xticks([0,90,180,270,360])
    ax2.scatter(np.array(AEDEY_obs.WSP),np.array(AEDEY_obs.TEMP))
    ax2.scatter(np.array(GJOGR_obs.WSP),np.array(GJOGR_obs.TEMP))
    ax2.scatter(np.array(BOLUN_obs.WSP),np.array(BOLUN_obs.TEMP))
    
    ax3.set_ylabel(r'air temp., $^{\circ}$C')
    ax3.set_xlabel('rh, %')
    #ax2.set_xticks([0,90,180,270,360])
    ax3.scatter(np.array(AEDEY_obs.RH),np.array(AEDEY_obs.TEMP))
    ax3.scatter(np.array(GJOGR_obs.RH),np.array(GJOGR_obs.TEMP))
    ax3.scatter(np.array(BOLUN_obs.RH),np.array(BOLUN_obs.TEMP))
    
    plt.tight_layout()
    plt.savefig('D:/Project/Figures/PNG/306/obs_only/gsta_raw_scatter_atmostate_full_month.png')

#%%
plot_scatters_6h=True
if plot_scatters_6h:
    ## take 6 hour averages
    window = 6
    AEDEY_wsp_means = np.mean(np.array(AEDEY_obs.WSP).reshape(-1, window), axis=1)
    GJOGR_wsp_means = np.mean(np.array(GJOGR_obs.WSP).reshape(-1, window), axis=1)
    BOLUN_wsp_means = np.mean(np.array(BOLUN_obs.WSP).reshape(-1, window), axis=1)
    
    AEDEY_wdir_means = np.mean(np.array(AEDEY_obs.WDIR).reshape(-1, window), axis=1)
    GJOGR_wdir_means = np.mean(np.array(GJOGR_obs.WDIR).reshape(-1, window), axis=1)
    BOLUN_wdir_means = np.mean(np.array(BOLUN_obs.WDIR).reshape(-1, window), axis=1)
    
    AEDEY_temp_means = np.mean(np.array(AEDEY_obs.TEMP).reshape(-1, window), axis=1)
    GJOGR_temp_means = np.mean(np.array(GJOGR_obs.TEMP).reshape(-1, window), axis=1)
    BOLUN_temp_means = np.mean(np.array(BOLUN_obs.TEMP).reshape(-1, window), axis=1)
    
    AEDEY_rh_means = np.mean(np.array(AEDEY_obs.RH).reshape(-1, window), axis=1)
    GJOGR_rh_means = np.mean(np.array(GJOGR_obs.RH).reshape(-1, window), axis=1)
    BOLUN_rh_means = np.mean(np.array(BOLUN_obs.RH).reshape(-1, window), axis=1)
    
    fig,((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2, figsize = (15,15))
    
    ax0.set_ylabel(r'wsp, ms$^{-1}$')
    ax0.set_xlabel('wdir')
    ax0.set_xticks([0,90,180,270,360])
    ax0.scatter(AEDEY_wdir_means, AEDEY_wsp_means, label = 'AEDEY')
    ax0.scatter(GJOGR_wdir_means, GJOGR_wsp_means, label = 'GJOGR')
    ax0.scatter(BOLUN_wdir_means, BOLUN_wsp_means, label = 'BOLUN')
    ax0.legend()
    
    ax1.set_ylabel(r'air temp., $^{\circ}$C')
    ax1.set_xlabel('wdir')
    ax1.set_xticks([0,90,180,270,360])
    ax1.scatter(AEDEY_wdir_means, AEDEY_temp_means, label = 'AEDEY')
    ax1.scatter(GJOGR_wdir_means, GJOGR_temp_means, label = 'GJOGR')
    ax1.scatter(BOLUN_wdir_means, BOLUN_temp_means, label = 'BOLUN')
    
    ax2.set_ylabel(r'air temp., $^{\circ}$C')
    ax2.set_xlabel('wsp, ms$^{-1}$')
    #ax2.set_xticks([0,90,180,270,360])
    ax2.scatter(AEDEY_wsp_means, AEDEY_temp_means, label = 'AEDEY')
    ax2.scatter(GJOGR_wsp_means, GJOGR_temp_means, label = 'GJOGR')
    ax2.scatter(BOLUN_wsp_means, BOLUN_temp_means, label = 'BOLUN')
    
    ax3.set_ylabel(r'air temp., $^{\circ}$C')
    ax3.set_xlabel('rh, %')
    #ax2.set_xticks([0,90,180,270,360])
    ax3.scatter(AEDEY_rh_means, AEDEY_temp_means, label = 'AEDEY')
    ax3.scatter(GJOGR_rh_means, GJOGR_temp_means, label = 'GJOGR')
    ax3.scatter(BOLUN_rh_means, BOLUN_temp_means, label = 'BOLUN')
    
    plt.tight_layout()
    plt.savefig('D:/Project/Figures/PNG/306/obs_only/gsta_6houravr_scatter_atmostate_full_month.png')

#%%  
plot_raw_heatmap=True
if plot_raw_heatmap:
    temp = np.concatenate((np.array(AEDEY_obs.TEMP),np.array(GJOGR_obs.TEMP),np.array(BOLUN_obs.TEMP)))
    wsp = np.concatenate((np.array(AEDEY_obs.WSP),np.array(GJOGR_obs.WSP),np.array(BOLUN_obs.WSP)))
    wdir = np.concatenate((np.array(AEDEY_obs.WDIR),np.array(GJOGR_obs.WDIR),np.array(BOLUN_obs.WDIR)))
    relh = np.concatenate((np.array(AEDEY_obs.RH),np.array(GJOGR_obs.RH),np.array(BOLUN_obs.RH)))
    
    fig,((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2, figsize = (15,15))
    
    kwargs = {'bins':40, 'cmap':'Blues_r'}
    
    ax0.set_ylabel(r'wsp, ms$^{-1}$')
    ax0.set_xlabel('wdir')
    ax0.set_xticks([0,90,180,270,360])
    ax0.hist2d(wdir,wsp,**kwargs)
    
    ax1.set_ylabel(r'air temp., $^{\circ}$C')
    ax1.set_xlabel('wdir')
    ax1.set_xticks([0,90,180,270,360])
    ax1.hist2d(wdir,temp,**kwargs)
    
    ax2.set_ylabel(r'air temp., $^{\circ}$C')
    ax2.set_xlabel('wsp, ms$^{-1}$')
    ax2.hist2d(wsp,temp,**kwargs)
    
    ax3.set_ylabel(r'air temp., $^{\circ}$C')
    ax3.set_xlabel('rh, %')
    ax3.hist2d(relh[~np.isnan(relh)],temp[~np.isnan(relh)],**kwargs)
    
    plt.tight_layout()
    plt.savefig('D:/Project/Figures/PNG/306/obs_only/gsta_raw_2dhist_atmostate_full_month.png')