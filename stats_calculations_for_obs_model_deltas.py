# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 17:27:08 2022

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
import stat_val_functions

#%%
resolutions = ['0p5km']#,'1p5km','4p4km']
suite = 'u-cf117'
config = 'LONGTAIL'
    
for res in resolutions:
    try:
        df_model_all = pd.read_csv(f'D:/Project/Model_Data/{suite}/{config}_Matched_interpolated_{res}_values_60s_306.txt')
    except:
        df_model_all = pd.read_csv(f'D:/Project/Model_Data/{suite}/{config}_Matched_interpolated_values_{res}_60s_306.txt')
    print(df_model_all)
    
    #%% load obs data
    obs_filename = f'IGP_flights_database_obs_60s_306.txt'
    path = 'D:/Project/Obvs_Data/Databases/'
    
    df_obs_all = pd.read_csv(f'{path}{obs_filename}', delimiter = ' ')
    #df_obs = df_obs_all[['wsp','w','theta','q_BUCK','windstress', 'sh','lh']]
    print(df_obs_all.columns)
    #%%
    for leg in [1,4,7,8,9,13,15]:
        
        if leg == 1:
            name = 'mountain_aloft'
            condition = (df_obs_all['legno'] == 1) | (df_obs_all['legno'] == 2)
        elif leg == 4:
            name = 'direct_lee'
            condition = (df_obs_all['legno'] == 4) | (df_obs_all['legno'] == 5) | (df_obs_all['legno'] == 6) | (df_obs_all['legno'] == 10)
        elif leg == 7:
            name = 'with_wind'
            condition = (df_obs_all['legno'] == 7) | (df_obs_all['legno'] == 8) | (df_obs_all['legno'] == 9) 
        elif leg == 13:
            name = 'leg13'
            condition = (df_obs_all['legno'] == 13) 
        elif leg == 15:
            name = 'leg15'
            condition = (df_obs_all['legno'] == 15) 
        
        df_obs = df_obs_all[condition]
        df_model = df_model_all[condition]
        #%% make subtractions
        ## all deltas defined by OBS - UM
        delta_wsp = -(np.array(df_obs['wsp']) - np.array(df_model['wsp']))
        delta_w   = -(np.array(df_obs['w'])   - np.array(df_model['w']))
        delta_theta = -(np.array(df_obs['theta']) - np.array(df_model['theta']))
        delta_q   = -(np.array(df_obs['q_BUCK']) - np.array(df_model['q'])*1000)
        delta_ws = -(np.array(df_obs['windstress']) - np.array(df_model['ws']))
        delta_sh = -(np.array(df_obs['sh']) - np.array(df_model['sh']))
        delta_lh = -(np.array(df_obs['lh']) - np.array(df_model['lh']))
        
        print('----------------------- testing q values\n',delta_q, np.array(df_obs['q_BUCK']), np.array(df_model['q']))
        #%% make new dataframe for delta and save it
        coloumn_names = ['wsp','w','theta','q','ws','sh','lh']
        data = np.transpose([delta_wsp, delta_w, delta_theta, delta_q, delta_ws, delta_sh, delta_lh])
        df_deltas = pd.DataFrame(data, columns = coloumn_names)
        #print(df_deltas)
        df_deltas.to_csv(f'D:/Project/Model_Data/{suite}/stat_val/deltas/{res}/{config}_obs_matched_deltas_{res}_{name}_306.csv')
        
        #%% compute statistics
        ## we won't use DataFrame.describe, but persoalise it a little
        mean = np.nanmean(data, axis = 0)
        std  = np.nanstd(data, axis = 0)
        median = np.nanmedian(data, axis = 0)
        rms = (np.nanmean(data**2, axis = 0))**0.5
        print(mean,'\n', std,'\n', median, '\n',rms)
        
        statistics =  [mean,std,median,rms]
        labeled_index = ['mean','std','median','rms']
        df_stats = pd.DataFrame(statistics, columns = coloumn_names, index = labeled_index)
        print(df_stats)
        df_stats.to_csv(f'D:/Project/Model_Data/{suite}/stat_val/deltas/{res}/{config}_obs_matched_deltas_statistics_{res}_{name}_306.csv')
        