# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 14:59:12 2021

@author: kse18nru
"""
#%% module imports


import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
#import stat_val_functions

#%% globals
suite = 'u-cf117'
flight = 306
matplotlib.rcParams.update({'font.size': 24})
experiment = 'LONGTAIL'
config = 'LONGTAIL'
#%% load obs data
if flight == 306:
    obs_filename = f'IGP_flights_database_obs_60s_{flight}.txt'
    path = 'D:/Project/Obvs_Data/Databases/'
    
    df_obs = pd.read_csv(f'{path}{obs_filename}', delimiter = ' ')
    
    db_time = np.array(df_obs['starttime'])
    db_alt = np.array(df_obs['altgps'])
    db_theta = np.array(df_obs['theta'])
    db_legno = np.array(df_obs['legno'])
    
    flight_suffix = ''
    df_model = pd.read_csv(f'D:/Project/Model_Data/{suite}/{config}_Matched_interpolated_values_0p5km_60s_{flight}.txt')
    um_theta = np.array(df_model['theta'])
    
    delta_theta = um_theta - db_theta
    
    #df_model_4p4 = pd.read_csv(f'D:/Project/Model_Data/{suite}/{config}_Matched_interpolated_4p4km_values_60s_306.txt')
    #um_theta_4p4 = np.array(df_model_4p4['theta'])
    #delta_theta_4p4 = um_theta_4p4 - db_theta
    
    #%%
    lee_start = 51; lee_end = 148
    lee_time = db_time[lee_start:lee_end]/(60*60)
    lee_median = np.median(delta_theta[lee_start:lee_end])*np.ones(len(lee_time))
    print('median',lee_median[0], '+/-',np.std(delta_theta[lee_start:lee_end]))
    #print(np.std(delta_theta[lee_start:lee_end]))
    isochron_1 = np.arange(-4,5)
    lee_start_array = lee_time[0]*np.ones(len(isochron_1))
    
    lee_end_array = lee_time[-1]*np.ones(len(isochron_1))
    
    print(lee_time[0], lee_time[-1])
    
    #%% plotting
    fig, ax = plt.subplots(1,1, figsize = (20,5))
    
    
    ax1 = ax.twinx() # instantiate a second axes that shares the same x-axis
    ax1.fill_between(db_time/(60*60), 0, db_alt, alpha = 0.2, color = 'dimgrey')
    ax1.set_yticks([0,500,1000,1500,2000])
    ax1.set_yticklabels([0, 0.5, 1.0, 1.5, 2.0])
    ax1.set_ylabel('Altitude km')
    
    ax.plot(db_time/(60*60), np.ones(len(db_time))*0, color = 'k')
    ax.plot(db_time/(60*60),delta_theta, color = 'orangered', label = 'RA1M 0p5km')
    #ax.plot(db_time/(60*60),delta_theta_4p4, color = 'olive' , label = 'RA1M 4p4km')
    ax.set_ylabel(r'$\Delta\theta$ K')
    ax.set_xticks([14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0])
    ax.set_xticklabels(['14:00', '14:30', '15:00', '15:30', '16:00', '16:30', '17:00'])
    #ax.set_xlabel('Clocktime')
    ax.set_xlim(left = db_time[0]/(60*60), right = db_time[-1]/(60*60))
    
    ax.plot(lee_time, lee_median, color = 'blue', label = r'0p5km median = $(-0.96\pm0.51)$K')
    #ax.legend(fontsize = 24)
    ax.plot(lee_start_array, isochron_1, color = 'k', linestyle = 'dashed')
    ax.plot(lee_end_array, isochron_1, color = 'k', linestyle = 'dashed')
    
    ax.text(14.7, 3.0, '14:41')
    ax.text(16.4, 3.0, '16:23')
    #ax.text(15.7, -2.8, r'median = $(-1.15\pm0.51)$K')
    theta_cha = r'$\theta$'
    ax.set_title(f'$\Delta${theta_cha} on matched timeseries - {config} Case 2 - 60s obs')
    ax.legend(fontsize = 19)
    plt.tight_layout()
    fig.savefig(f'D:/Project/Figures/PNG/306/{experiment}/{suite}/stat_val/{config}_0p5km_delta_theta_timeseries_with interval_median_flt{flight}.png')
    #%%
if flight == 301:
    obs_filename = f'IGP_flights_database_obs_60s_{flight}.txt'
    path = 'D:/Project/Obvs_Data/Databases/'
    
    df_obs = pd.read_csv(f'{path}{obs_filename}', delimiter = ' ')
    
    db_time = np.array(df_obs['starttime'])
    db_alt = np.array(df_obs['altgps'])
    db_theta = np.array(df_obs['theta'])
    db_legno = np.array(df_obs['legno'])
    
    flight_suffix = ''
    df_model = pd.read_csv(f'D:/Project/Model_Data/{suite}/Matched_interpolated_0p5km_values_60s_301.txt')
    um_theta = np.array(df_model['theta'])
    
    delta_theta = um_theta - db_theta
    
    df_model_4p4 = pd.read_csv(f'D:/Project/Model_Data/{suite}/Matched_interpolated_4p4km_values_60s_301.txt')
    um_theta_4p4 = np.array(df_model_4p4['theta'])
    delta_theta_4p4 = um_theta_4p4 - db_theta
    print(um_theta - um_theta_4p4)
    #%%
    lee_start = 51; lee_end = 148
    lee_time = db_time[lee_start:lee_end]/(60*60)
    lee_median = np.median(delta_theta[lee_start:lee_end])*np.ones(len(lee_time))
    print(lee_median[0])
    print(np.std(delta_theta[lee_start:lee_end]))
    isochron_1 = np.arange(-4,5)
    lee_start_array = lee_time[0]*np.ones(len(isochron_1))
    
    lee_end_array = lee_time[-1]*np.ones(len(isochron_1))
    
    print(lee_time[0], lee_time[-1])
    
    #%% plotting
    fig, ax = plt.subplots(1,1, figsize = (20,5))
    
    
    ax1 = ax.twinx() # instantiate a second axes that shares the same x-axis
    ax1.fill_between(db_time/(60*60), 0, db_alt, alpha = 0.2, color = 'dimgrey')
    ax1.set_yticks([0,500,1000,1500,2000])
    ax1.set_yticklabels([0, 0.5, 1.0, 1.5, 2.0])
    ax1.set_ylabel('Altitude km')
    
    ax.plot(db_time/(60*60), np.ones(len(db_time))*0, color = 'k')
    ax.plot(db_time/(60*60),delta_theta, color = 'orangered', label = 'RA1M 0p5km')
    ax.plot(db_time/(60*60),delta_theta_4p4, color = 'olive', label = 'RA1M 4p4km' )
    ax.set_ylabel(r'$\Delta\theta$ K')
    ax.set_xticks([12.0,12.5,13.0,13.5,14.0, 14.5, 15.0, 15.5, 16.0, 16.5, 17.0])
    ax.set_xticklabels(['12:00','12:30','13:00','13:30','14:00', '14:30', '15:00', '15:30', '16:00', '16:30', '17:00'])
    #ax.set_xlabel('Clocktime')
    ax.set_xlim(left = db_time[0]/(60*60), right = db_time[-1]/(60*60))
    ax.legend(fontsize = 19)
    #ax.plot(lee_time, lee_median, color = 'blue', label = r'median = $(-1.15\pm0.51)$K')
    #ax.legend()
    #ax.plot(lee_start_array, isochron_1, color = 'k', linestyle = 'dashed')
    #ax.plot(lee_end_array, isochron_1, color = 'k', linestyle = 'dashed')
    
    #ax.text(14.7, 3.0, '14:41')
    #ax.text(16.4, 3.0, '16:23')
    #ax.text(15.7, -2.8, r'median = $(-1.15\pm0.51)$K')
    ax.set_title(r'$\Delta\theta$ matched timeseries - RA1M Case 1 - 60s obs')
    plt.tight_layout()
    fig.savefig(f'D:/Project/Figures/PNG/{flight}/{suite}/stat_val/RA1M_0p5kmand4p4km_delta_theta_timeseries_flt{flight}.png')