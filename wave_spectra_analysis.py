# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 13:07:37 2022

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

#%% globals
suite = 'u-cc134'
experiment = 'CONTROL'
flight = 306
matplotlib.rcParams.update({'font.size': 24})
#%% load obs data
if flight == 306:
    obs_filename = f'IGP_flights_database_obs_sr_{flight}.txt'
    path = 'D:/Project/Obvs_Data/Databases/'
    
    df_obs = pd.read_csv(f'{path}{obs_filename}', delimiter = ' ')
    print(df_obs.columns)
    
    df_model = pd.read_csv(f'D:/Project/Model_Data/{suite}/Matched_interpolated_values_0p5km_sr_306.txt')
    
    print(df_model)
    
    #%%
    
    
    fig,ax = plt.subplots(1,1, figsize = (10,10))
    
    for i in range(10):
        leg = i + 1
        df_lat = df_obs.lat[df_obs.legno == leg]
        df_lon = df_obs.lon[df_obs.legno == leg]
        ax.plot(df_lon,df_lat, label = leg)
    plt.legend()
    plt.show()
    
    
    #%%
    fig,ax = plt.subplots(1,1, figsize = (10,10))
    leg = 5
    df_lat = df_obs.lat[df_obs.legno == leg]
    df_lon = df_obs.lon[df_obs.legno == leg]
    #ax.scatter(df_lon[df_obs.altgps >280],df_lat[df_obs.altgps >280], label = '9')
    #ax.scatter(df_lon[df_obs.altgps <200],df_lat[df_obs.altgps <200])
    ax.scatter(df_lon, df_lat, c = df_obs.altgps[df_obs.legno == leg])
    #%% get rough spatial spacing of data points
    ptp_dists = []
    lat = np.array(df_lat); lon = np.array(df_lon)
    for i in range(len(df_lon)-1):
        step = toolkit.haversine(lat[i], lon[i], lat[i+1], lon[i+1]) # returns distance on globe in km
        ptp_dists.append(step)
    ptp_dists = np.array(ptp_dists)
    ptp_mean = np.mean(ptp_dists)
    ptp_std   = np.std(ptp_dists)
    print(f'PTP distance =( {ptp_mean} +/- {ptp_std} )km')
    d = ptp_mean
    
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
    
    #%% filtering
    #raw_obs_w = np.array(df_obs.w)
    raw_obs_theta = np.array(df_obs.theta)
    SG_poly3_win5_obs_theta = filtering.savitzky_golay(raw_obs_theta, window_size = 3, order = 1)
    
    #raw_um_w = np.array(df_model.w)
    raw_um_theta = np.array(df_model.theta)
    SG_poly3_win5_um_theta = filtering.savitzky_golay(raw_um_theta, window_size = 3, order = 1)
    #%% fft
    fft_plots = False
    if fft_plots:
        freq_spacing = fftfreq(len(raw_obs_theta), d) # due to spatial series, this is the wavenumber
        print(freq_spacing)
        wavelengths = 1/freq_spacing
        print(wavelengths)
        raw_obs_theta_fft = fftshift(fft(np.array(raw_obs_theta)))
        raw_um_theta_fft = fftshift(fft(np.array(raw_um_theta)))
        SG_poly3_win5_obs_theta_fft = fftshift(fft(SG_poly3_win5_obs_theta))
        SG_poly3_win5_um_theta_fft = fftshift(fft(SG_poly3_win5_um_theta))
        fig, (ax0,ax1) = plt.subplots(2,1, figsize = (10,6.5))
        
        ax0.plot(np.log10(wavelengths),np.abs(raw_obs_theta_fft), label = 'obs sr')
        ax0.plot(np.log10(wavelengths),np.abs(raw_um_theta_fft), label = 'UM 0p5km')
        
        
        
        ax0.set_xticks(log_ticks)
        ax0.set_xticklabels(length_ticks_labels)
        ax0.set_xlim(left = np.log10(1))#np.log10(wavelengths[:int((len(wavelengths)/2))])[-1])
        ax0.set_xlabel('Wavelength, km')
        ax0.legend(fontsize = 16)
        
        
        ax1.plot(np.log10(wavelengths), np.abs(SG_poly3_win5_obs_theta_fft))
        ax1.plot(np.log10(wavelengths), np.abs(SG_poly3_win5_um_theta_fft))
        ax1.set_xticks(log_ticks)
        ax1.set_xticklabels(length_ticks_labels)
        ax1.set_xlim(left = np.log10(1))#np.log10(wavelengths[:int((len(wavelengths)/2))])[-1])
        ax1.set_xlabel('Wavelength, km')
        #ax1.set_ylim(top = 100)
        
        fig.suptitle('FFT analysis of Case 2 legs 7-9')
        plt.tight_layout
        #plt.savefig(f'D:/Project/Figures/PNG/{flight}/{experiment}/{suite}/simple_wave_fft_spectra_legs789_sr_cc134.png')
    
    #%%
    
    fig, ax0 = plt.subplots(1,1, figsize = (10,4))
    
    
    ## plot theta timeseries
    