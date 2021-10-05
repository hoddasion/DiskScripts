# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 16:21:32 2021

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

#%% oad obs data
obs_path = 'D:/Project/Obvs_Data/profiles/keflavik_04018BIKF_profiles_upto500mb.csv'
obs_df = pd.read_csv(obs_path)


#%% 
flight  = 306
if flight ==301:
    obs_times = ['20180312T0600Z','20180312T1200Z','20180312T1800Z']
    time_idcs = [11,23,35]
    #%% load model data
    suite = 'u-bu807'
    res = '0p5km'
    um_path = f'D:/Project/Model_Data/{suite}/nc/Control/'
    ## velocity diagnostics
    nwind_cube = iris.load_cube(f'{um_path}RA1M_{res}_um_m01s15i002_24hrs_h_301.nc', 'm01s15i002')
    ewind_cube = iris.load_cube(f'{um_path}RA1M_{res}_um_m01s15i003_24hrs_h_301.nc', 'm01s15i003')
    magwind_cube = (nwind_cube**2 + ewind_cube**2)**0.5 
    ## theta diagnostic
    theta_cube = iris.load_cube(f'{um_path}RA1M_{res}_um_air_potential_temperature_24hrs_i_301.nc', 'air_potential_temperature')
    
    #%% isolate coloumns in model data
    kef_lat = 63.96
    kef_lon = -22.60
    polelat = theta_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = theta_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    kef_x, kef_y = rotate_pole(np.array([kef_lon]), np.array([kef_lat]), polelon, polelat)
    kef_x = kef_x + 360
    print(kef_y, kef_x)
    sample_points = [('grid_longitude', kef_x), ('grid_latitude', kef_y)]
    theta_coloumns = trajectory.interpolate(theta_cube,sample_points, method =  'nearest')
    wind_coloumns = trajectory.interpolate(magwind_cube,sample_points, method =  'nearest')
    #%%
    i = 0
    for obs_time in obs_times:
        print(obs_time)
        idx = time_idcs[i]
        obs_profiles = obs_df[obs_df.Time == obs_time]
        
        obs_h = np.array(obs_profiles.HGHT).astype(np.float)
        
        obs_p = np.array(obs_profiles.PRES).astype(np.float)
        obs_theta = np.array(obs_profiles.THTA).astype(np.float)
        obs_speed =np.array(obs_profiles.SKNT).astype(np.float)*0.5144 # convert from knots to m/s
        
        fig, (ax0,ax1) = plt.subplots(1,2, figsize = (9,8))
        height_ticks = [1000,2000,3000,4000,5000]
        ## ax0: theta
        
        ax0.plot(obs_theta, obs_h, label = 'Kefl. Sounding') # plot obs profile
        ax0.plot(theta_coloumns.data[idx, :39], theta_coloumns.coords('altitude')[0].points[:39], label = 'UM')
        ax0.set_ylim(bottom = 0,top = 5100)
        ax0.set_yticks(height_ticks)
        ax0.set_yticklabels([1,2,3,4,5])
        ax0.set_ylabel('km')
        ax0.set_xlabel(r'$\theta$, $K$')
        ax0.legend()
        #ax0.set_title('Potential temperature')
        ## ax1: wind speed
        ax1.plot(obs_speed, obs_h, label = 'Kefl. Sounding')
        ax1.plot(wind_coloumns.data[idx, :39], wind_coloumns.coords('altitude')[0].points[:39], label = 'UM')
        ax1.set_ylim(bottom = 0, top = 5100)
        ax1.set_yticks(height_ticks)
        ax1.set_yticklabels([])
        ax1.set_xlabel(r'wsp, $ms^{-1}$')
        ax0.legend()
        #ax1.set_title('Windspeed')
        
        plt.tight_layout()
        plt.savefig(f'D:/Project/Figures/PDF/301/{suite}/Keflavik_profiles_{res}_{obs_time}_theta_wsp_301.pdf')
        fig.suptitle(f'Keflavig-UM {res} profiles - {obs_time}')
        plt.tight_layout()
        plt.savefig(f'D:/Project/Figures/PNG/301/{suite}/Keflavik_profiles_{res}_{obs_time}_theta_wsp_301.png')
        i =+ 1
#%%
if flight == 306:
    obs_times = ['20180319T0600Z','20180319T1200Z','20180319T1800Z']
    time_idcs = [11,23,35]
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
    
    #%% isolate coloumns in model data
    kef_lat = 63.96
    kef_lon = -22.60
    polelat = theta_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = theta_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    kef_x, kef_y = rotate_pole(np.array([kef_lon]), np.array([kef_lat]), polelon, polelat)
    kef_x = kef_x + 360
    print(kef_y, kef_x)
    sample_points = [('grid_longitude', kef_x), ('grid_latitude', kef_y)]
    theta_coloumns = trajectory.interpolate(theta_cube,sample_points, method =  'nearest')
    wind_coloumns = trajectory.interpolate(magwind_cube,sample_points, method =  'nearest')
    #%%
    i = 0
    for obs_time in obs_times:
        print(obs_time)
        idx = time_idcs[i]
        obs_profiles = obs_df[obs_df.Time == obs_time]
        
        obs_h = np.array(obs_profiles.HGHT).astype(np.float)
        
        obs_p = np.array(obs_profiles.PRES).astype(np.float)
        obs_theta = np.array(obs_profiles.THTA).astype(np.float)
        obs_speed =np.array(obs_profiles.SKNT).astype(np.float)*0.5144 # convert from knots to m/s
        
        fig, (ax0,ax1) = plt.subplots(1,2, figsize = (9,8))
        height_ticks = [1000,2000,3000,4000,5000]
        ## ax0: theta
        
        ax0.plot(obs_theta, obs_h, label = 'Kefl. Sounding') # plot obs profile
        ax0.plot(theta_coloumns.data[idx, :39], theta_coloumns.coords('altitude')[0].points[:39], label = 'UM')
        ax0.set_ylim(bottom = 0,top = 5100)
        ax0.set_yticks(height_ticks)
        ax0.set_yticklabels([1,2,3,4,5])
        ax0.set_ylabel('km')
        ax0.set_xlabel(r'$\theta$, $K$')
        ax0.legend()
        #ax0.set_title('Potential temperature')
        ## ax1: wind speed
        ax1.plot(obs_speed, obs_h, label = 'Kefl. Sounding')
        ax1.plot(wind_coloumns.data[idx, :39], wind_coloumns.coords('altitude')[0].points[:39], label = 'UM')
        ax1.set_ylim(bottom = 0, top = 5100)
        ax1.set_yticks(height_ticks)
        ax1.set_yticklabels([])
        ax1.set_xlabel(r'wsp, $ms^{-1}$')
        ax0.legend()
        #ax1.set_title('Windspeed')
        
        plt.tight_layout()
        plt.savefig(f'D:/Project/Figures/PDF/{flight}/{suite}/Keflavik_profiles_{res}_{obs_time}_theta_wsp_{flight}.pdf')
        fig.suptitle(f'Keflavig-UM {res} profiles - {obs_time}')
        plt.tight_layout()
        plt.savefig(f'D:/Project/Figures/PNG/{flight}/{suite}/Keflavik_profiles_{res}_{obs_time}_theta_wsp_{flight}.png')
        i =+ 1
    