# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 19:15:11 2022

@author: kse18nru
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
#%%
matplotlib.rcParams.update({'font.size': 28})
#%% load ALL the DATA
if False:
    
    df_mountain_4p4 = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/4p4km/RA1M_obs_matched_deltas_statistics_4p4km_mountain_aloft_306.csv',index_col=[0])[['wsp','theta','q']]
    df_mountain_1p5 = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/1p5km/RA1M_obs_matched_deltas_statistics_1p5km_mountain_aloft_306.csv',index_col=[0])[['wsp','theta','q']]
    df_mountain_0p5 = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/0p5km/RA1M_obs_matched_deltas_statistics_0p5km_mountain_aloft_306.csv',index_col=[0])[['wsp','theta','q']]
    
    df_directlee_4p4 = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/4p4km/RA1M_obs_matched_deltas_statistics_4p4km_direct_lee_306.csv',index_col=[0])[['wsp','theta','q']]
    df_directlee_1p5 = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/1p5km/RA1M_obs_matched_deltas_statistics_1p5km_direct_lee_306.csv',index_col=[0])[['wsp','theta','q']]
    df_directlee_0p5 = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/0p5km/RA1M_obs_matched_deltas_statistics_0p5km_direct_lee_306.csv',index_col=[0])[['wsp','theta','q']]
    
    df_withwind_4p4 = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/4p4km/RA1M_obs_matched_deltas_statistics_4p4km_with_wind_306.csv',index_col=[0])[['wsp','theta','q']]
    df_withwind_1p5 = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/1p5km/RA1M_obs_matched_deltas_statistics_1p5km_with_wind_306.csv',index_col=[0])[['wsp','theta','q']]
    df_withwind_0p5 = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/0p5km/RA1M_obs_matched_deltas_statistics_0p5km_with_wind_306.csv',index_col=[0])[['wsp','theta','q']]
    
    df_leg13_4p4 = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/4p4km/RA1M_obs_matched_deltas_statistics_4p4km_leg13_306.csv',index_col=[0])[['wsp','theta','q']]
    df_leg13_1p5 = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/1p5km/RA1M_obs_matched_deltas_statistics_1p5km_leg13_306.csv',index_col=[0])[['wsp','theta','q']]
    df_leg13_0p5 = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/0p5km/RA1M_obs_matched_deltas_statistics_0p5km_leg13_306.csv',index_col=[0])[['wsp','theta','q']]
    
    df_leg15_4p4 = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/4p4km/RA1M_obs_matched_deltas_statistics_4p4km_leg15_306.csv',index_col=[0])[['wsp','theta','q']]
    df_leg15_1p5 = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/1p5km/RA1M_obs_matched_deltas_statistics_1p5km_leg15_306.csv',index_col=[0])[['wsp','theta','q']]
    df_leg15_0p5 = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/0p5km/RA1M_obs_matched_deltas_statistics_0p5km_leg15_306.csv',index_col=[0])[['wsp','theta','q']]
    
    #%%
    print(df_mountain_1p5)
    print(df_mountain_1p5.wsp['mean'])
    
    #%% group data for means
    
    ### wsp
    wsp_mean_0p5 = [df_mountain_0p5.wsp['mean'],df_directlee_0p5.wsp['mean'],df_withwind_0p5.wsp['mean'],df_leg13_0p5.wsp['mean'],df_leg15_0p5.wsp['mean']]
    wsp_mean_1p5 = [df_mountain_1p5.wsp['mean'],df_directlee_1p5.wsp['mean'],df_withwind_1p5.wsp['mean'],df_leg13_1p5.wsp['mean'],df_leg15_1p5.wsp['mean']]
    wsp_mean_4p4 = [df_mountain_4p4.wsp['mean'],df_directlee_4p4.wsp['mean'],df_withwind_4p4.wsp['mean'],df_leg13_4p4.wsp['mean'],df_leg15_4p4.wsp['mean']]
    ### theta
    theta_mean_0p5 = [df_mountain_0p5.theta['mean'],df_directlee_0p5.theta['mean'],df_withwind_0p5.theta['mean'],df_leg13_0p5.theta['mean'],df_leg15_0p5.theta['mean']]
    theta_mean_1p5 = [df_mountain_1p5.theta['mean'],df_directlee_1p5.theta['mean'],df_withwind_1p5.theta['mean'],df_leg13_1p5.theta['mean'],df_leg15_1p5.theta['mean']]
    theta_mean_4p4 = [df_mountain_4p4.theta['mean'],df_directlee_4p4.theta['mean'],df_withwind_4p4.theta['mean'],df_leg13_4p4.theta['mean'],df_leg15_4p4.theta['mean']]
    ### q
    q_mean_0p5 = [df_mountain_0p5.q['mean'],df_directlee_0p5.q['mean'],df_withwind_0p5.q['mean'],df_leg13_0p5.q['mean'],df_leg15_0p5.q['mean']]
    q_mean_1p5 = [df_mountain_1p5.q['mean'],df_directlee_1p5.q['mean'],df_withwind_1p5.q['mean'],df_leg13_1p5.q['mean'],df_leg15_1p5.q['mean']]
    q_mean_4p4 = [df_mountain_4p4.q['mean'],df_directlee_4p4.q['mean'],df_withwind_4p4.q['mean'],df_leg13_4p4.q['mean'],df_leg15_4p4.q['mean']]
    
    
    
    width = 0.2
    x_locations_center = np.array([1,2,3,4,5])
    x_locations_left = x_locations_center - width
    x_locations_right = x_locations_center + width
    #%% plot figure
    
    fig, (ax0,ax1,ax2) = plt.subplots(3,1, figsize = (20,16))
    
    ## ax0 for wsp
    ax0.bar(x_locations_left, wsp_mean_0p5, width, label = '0p5km')
    ax0.bar(x_locations_center, wsp_mean_1p5, width, label = '1p5km')
    ax0.bar(x_locations_right, wsp_mean_4p4, width, label = '4p4km')
    ax0.plot([0.5,5.5], [0,0], color = 'k', linestyle = 'dashed')
    ax0.set_xlim(left = 0.5, right = 5.5)
    ax0.set_xticks([1,2,3,4,5])
    ax0.set_xticklabels([])
    ax0.set_yticks([-3,-2,-1,0,1,2,3])
    ax0.set_ylabel(r'$\delta U$, ms$^{-1}$')
    #ax0.legend()
    ax0.yaxis.grid()
    
    ## ax1 for theta
    ax1.bar(x_locations_left, theta_mean_0p5, width)
    ax1.bar(x_locations_center, theta_mean_1p5, width)
    ax1.bar(x_locations_right, theta_mean_4p4, width)
    ax1.plot([0.5,5.5], [0,0], color = 'k', linestyle = 'dashed')
    ax1.set_xlim(left = 0.5, right = 5.5)
    ax1.set_xticks([1,2,3,4,5])
    ax1.set_xticklabels([])
    ax1.set_ylabel(r'$\delta\theta$, K')
    ax1.set_yticks([-1.5,-1,-0.5,0,0.5,1,1.5])
    ax1.yaxis.grid()
    
     ## ax2 for q
    ax2.bar(x_locations_left, q_mean_0p5, width)
    ax2.bar(x_locations_center, q_mean_1p5, width)
    ax2.bar(x_locations_right, q_mean_4p4, width)
    ax2.plot([0.5,5.5], [0,0], color = 'k', linestyle = 'dashed')
    ax2.set_xlim(left = 0.5, right = 5.5)
    ax2.set_xticks([1,2,3,4,5])
    ax2.set_xticklabels(['A','B','C','D','E'])
    ax2.set_ylabel(r'$\delta q$, gkg$^{-1}$')
    ax2.set_yticks(-np.array([-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5]))
    ax2.yaxis.grid()
    
    fig.suptitle('Case 2 - RA1M - Mean bias')
    plt.tight_layout()
    
    if False:
        plt.savefig('D:/Project/Figures/PNG/306/CONTROL/u-cc134/stat_val/Case2_RA1M_mean_bias_bar_charts_Uthetaq.png')
    
    #%% group data for rms
    
    ### wsp
    wsp_rms_0p5 = [df_mountain_0p5.wsp['rms'],df_directlee_0p5.wsp['rms'],df_withwind_0p5.wsp['rms'],df_leg13_0p5.wsp['rms'],df_leg15_0p5.wsp['rms']]
    wsp_rms_1p5 = [df_mountain_1p5.wsp['rms'],df_directlee_1p5.wsp['rms'],df_withwind_1p5.wsp['rms'],df_leg13_1p5.wsp['rms'],df_leg15_1p5.wsp['rms']]
    wsp_rms_4p4 = [df_mountain_4p4.wsp['rms'],df_directlee_4p4.wsp['rms'],df_withwind_4p4.wsp['rms'],df_leg13_4p4.wsp['rms'],df_leg15_4p4.wsp['rms']]
    ### theta
    theta_rms_0p5 = [df_mountain_0p5.theta['rms'],df_directlee_0p5.theta['rms'],df_withwind_0p5.theta['rms'],df_leg13_0p5.theta['rms'],df_leg15_0p5.theta['rms']]
    theta_rms_1p5 = [df_mountain_1p5.theta['rms'],df_directlee_1p5.theta['rms'],df_withwind_1p5.theta['rms'],df_leg13_1p5.theta['rms'],df_leg15_1p5.theta['rms']]
    theta_rms_4p4 = [df_mountain_4p4.theta['rms'],df_directlee_4p4.theta['rms'],df_withwind_4p4.theta['rms'],df_leg13_4p4.theta['rms'],df_leg15_4p4.theta['rms']]
    ### q
    q_rms_0p5 = [df_mountain_0p5.q['rms'],df_directlee_0p5.q['rms'],df_withwind_0p5.q['rms'],df_leg13_0p5.q['rms'],df_leg15_0p5.q['rms']]
    q_rms_1p5 = [df_mountain_1p5.q['rms'],df_directlee_1p5.q['rms'],df_withwind_1p5.q['rms'],df_leg13_1p5.q['rms'],df_leg15_1p5.q['rms']]
    q_rms_4p4 = [df_mountain_4p4.q['rms'],df_directlee_4p4.q['rms'],df_withwind_4p4.q['rms'],df_leg13_4p4.q['rms'],df_leg15_4p4.q['rms']]
    
    #%% plot figure
    
    fig, (ax0,ax1,ax2) = plt.subplots(3,1, figsize = (20,16))
    
    ## ax0 for wsp
    ax0.bar(x_locations_left, wsp_rms_0p5, width)
    ax0.bar(x_locations_center, wsp_rms_1p5, width)
    ax0.bar(x_locations_right, wsp_rms_4p4, width)
    #ax0.plot([0.5,5.5], [0,0], color = 'k', linestyle = 'dashed')
    ax0.set_xlim(left = 0.5, right = 5.5)
    ax0.set_xticks([1,2,3,4,5])
    ax0.set_xticklabels(['A','B','C','D','E'])
    ax0.set_yticks([0,2,4,6,8])
    ax0.yaxis.grid()
    ax0.set_ylabel(r'$\delta U$, ms$^{-1}$')
    
    ## ax1 for theta
    ax1.bar(x_locations_left, theta_rms_0p5, width)
    ax1.bar(x_locations_center, theta_rms_1p5, width)
    ax1.bar(x_locations_right, theta_rms_4p4, width)
    #ax1.plot([0.5,5.5], [0,0], color = 'k', linestyle = 'dashed')
    ax1.set_xlim(left = 0.5, right = 5.5)
    ax1.set_xticks([1,2,3,4,5])
    ax1.set_xticklabels(['A','B','C','D','E'])
    ax1.set_yticks([0,0.5,1,1.5,2])
    ax1.yaxis.grid()
    ax1.set_ylabel(r'$\delta\theta$, K')
    
     ## ax2 for q
    ax2.bar(x_locations_left, q_rms_0p5, width)
    ax2.bar(x_locations_center, q_rms_1p5, width)
    ax2.bar(x_locations_right, q_rms_4p4, width)
    #ax2.plot([0.5,5.5], [0,0], color = 'k', linestyle = 'dashed')
    ax2.set_xlim(left = 0.5, right = 5.5)
    ax2.set_xticks([1,2,3,4,5])
    ax2.set_xticklabels(['A','B','C','D','E'])
    ax2.set_yticks([0,0.2,0.4,0.6])
    ax2.yaxis.grid()
    ax2.set_ylabel(r'$\delta q$, gkg$^{-1}$')
    
    fig.suptitle('Case 2 - RA1M - Root mean square bias')
    plt.tight_layout()
    
    if False:
        plt.savefig('D:/Project/Figures/PNG/306/CONTROL/u-cc134/stat_val/Case2_RA1M_rms_bias_bar_charts_Uthetaq.png')
        
#%%
if True:
    ## load all 0p5km CONTROL data
    df_mountain_cont = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/0p5km/RA1M_obs_matched_deltas_statistics_0p5km_mountain_aloft_306.csv',index_col=[0])[['wsp','theta','q']]
    df_directlee_cont = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/0p5km/RA1M_obs_matched_deltas_statistics_0p5km_direct_lee_306.csv',index_col=[0])[['wsp','theta','q']]
    df_withwind_cont = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/0p5km/RA1M_obs_matched_deltas_statistics_0p5km_with_wind_306.csv',index_col=[0])[['wsp','theta','q']]
    df_leg13_cont = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/0p5km/RA1M_obs_matched_deltas_statistics_0p5km_leg13_306.csv',index_col=[0])[['wsp','theta','q']]
    df_leg15_cont = pd.read_csv('D:/Project/Model_Data/u-cc134/stat_val/deltas/0p5km/RA1M_obs_matched_deltas_statistics_0p5km_leg15_306.csv',index_col=[0])[['wsp','theta','q']]
    
    ## load all 0p5km LONGTAIL data
    df_mountain_long = pd.read_csv('D:/Project/Model_Data/u-cf117/stat_val/deltas/0p5km/LONGTAIL_obs_matched_deltas_statistics_0p5km_mountain_aloft_306.csv',index_col=[0])[['wsp','theta','q']]
    df_directlee_long = pd.read_csv('D:/Project/Model_Data/u-cf117/stat_val/deltas/0p5km/LONGTAIL_obs_matched_deltas_statistics_0p5km_direct_lee_306.csv',index_col=[0])[['wsp','theta','q']]
    df_withwind_long = pd.read_csv('D:/Project/Model_Data/u-cf117/stat_val/deltas/0p5km/LONGTAIL_obs_matched_deltas_statistics_0p5km_with_wind_306.csv',index_col=[0])[['wsp','theta','q']]
    df_leg13_long = pd.read_csv('D:/Project/Model_Data/u-cf117/stat_val/deltas/0p5km/LONGTAIL_obs_matched_deltas_statistics_0p5km_leg13_306.csv',index_col=[0])[['wsp','theta','q']]
    df_leg15_long = pd.read_csv('D:/Project/Model_Data/u-cf117/stat_val/deltas/0p5km/LONGTAIL_obs_matched_deltas_statistics_0p5km_leg15_306.csv',index_col=[0])[['wsp','theta','q']]
    
    #%% group data for means
    
    ### wsp
    wsp_mean_cont = [df_mountain_cont.wsp['mean'],df_directlee_cont.wsp['mean'],df_withwind_cont.wsp['mean'],df_leg13_cont.wsp['mean'],df_leg15_cont.wsp['mean']]
    wsp_mean_long = [df_mountain_long.wsp['mean'],df_directlee_long.wsp['mean'],df_withwind_long.wsp['mean'],df_leg13_long.wsp['mean'],df_leg15_long.wsp['mean']]
    
    ### theta
    theta_mean_cont = [df_mountain_cont.theta['mean'],df_directlee_cont.theta['mean'],df_withwind_cont.theta['mean'],df_leg13_cont.theta['mean'],df_leg15_cont.theta['mean']]
    theta_mean_long = [df_mountain_long.theta['mean'],df_directlee_long.theta['mean'],df_withwind_long.theta['mean'],df_leg13_long.theta['mean'],df_leg15_long.theta['mean']]
    
    ### q
    q_mean_cont = [df_mountain_cont.q['mean'],df_directlee_cont.q['mean'],df_withwind_cont.q['mean'],df_leg13_cont.q['mean'],df_leg15_cont.q['mean']]
    q_mean_long = [df_mountain_long.q['mean'],df_directlee_long.q['mean'],df_withwind_long.q['mean'],df_leg13_long.q['mean'],df_leg15_long.q['mean']]
   
    width = 0.2
    x_locations_center = np.array([1,2,3,4,5])
    x_locations_left = x_locations_center - width
    x_locations_right = x_locations_center + width
    
    #%% plot figure
    
    fig, (ax0,ax1,ax2) = plt.subplots(3,1, figsize = (20,16))
    
    ## ax0 for wsp
    ax0.bar(x_locations_left, wsp_mean_cont, width, label = 'CONTROL')
    ax0.bar(x_locations_center, wsp_mean_long, width, label = 'LONGTAIL')
    ax0.plot([0.5,5.5], [0,0], color = 'k', linestyle = 'dashed')
    ax0.set_xlim(left = 0.5, right = 5.5)
    ax0.set_xticks([1,2,3,4,5])
    ax0.set_xticklabels([])
    ax0.set_yticks([-2,-1,0,1,2,3,4])
    ax0.set_ylabel(r'$\delta U$, ms$^{-1}$')
    ax0.legend()
    ax0.yaxis.grid()
    
    ## ax1 for theta
    ax1.bar(x_locations_left, theta_mean_cont, width)
    ax1.bar(x_locations_center, theta_mean_long, width)
    ax1.plot([0.5,5.5], [0,0], color = 'k', linestyle = 'dashed')
    ax1.set_xlim(left = 0.5, right = 5.5)
    ax1.set_xticks([1,2,3,4,5])
    ax1.set_xticklabels([])
    ax1.set_ylabel(r'$\delta\theta$, K')
    ax1.set_yticks([-1.5,-1,-0.5,0,0.5])
    ax1.yaxis.grid()
    
     ## ax2 for q
    ax2.bar(x_locations_left, q_mean_cont, width)
    ax2.bar(x_locations_center, q_mean_long, width)
    ax2.plot([0.5,5.5], [0,0], color = 'k', linestyle = 'dashed')
    ax2.set_xlim(left = 0.5, right = 5.5)
    ax2.set_xticks([1,2,3,4,5])
    ax2.set_xticklabels(['A','B','C','D','E'])
    ax2.set_ylabel(r'$\delta q$, gkg$^{-1}$')
    ax2.set_yticks(-np.array([-0.2,-0.1,0,0.1,0.2,0.3,0.4]))
    ax2.yaxis.grid()
    
    fig.suptitle('Case 2 - CONTROL/LONGTAIL - Mean bias')
    plt.tight_layout()
    
    if False:
        plt.savefig('D:/Project/Figures/PNG/306/LONGTAIL/u-cf117/stat_val/LONGTAIL_vs_CONTROL_mean_bias_bar_charts_Uthetaq.png')
        
        
    #%% group data for rms
    
    ### wsp
    wsp_rms_cont = [df_mountain_cont.wsp['rms'],df_directlee_cont.wsp['rms'],df_withwind_cont.wsp['rms'],df_leg13_cont.wsp['rms'],df_leg15_cont.wsp['rms']]
    wsp_rms_long = [df_mountain_long.wsp['rms'],df_directlee_long.wsp['rms'],df_withwind_long.wsp['rms'],df_leg13_long.wsp['rms'],df_leg15_long.wsp['rms']]
    
    ### theta
    theta_rms_cont = [df_mountain_cont.theta['rms'],df_directlee_cont.theta['rms'],df_withwind_cont.theta['rms'],df_leg13_cont.theta['rms'],df_leg15_cont.theta['rms']]
    theta_rms_long = [df_mountain_long.theta['rms'],df_directlee_long.theta['rms'],df_withwind_long.theta['rms'],df_leg13_long.theta['rms'],df_leg15_long.theta['rms']]
    
    ### q
    q_rms_cont = [df_mountain_cont.q['rms'],df_directlee_cont.q['rms'],df_withwind_cont.q['rms'],df_leg13_cont.q['rms'],df_leg15_cont.q['rms']]
    q_rms_long = [df_mountain_long.q['rms'],df_directlee_long.q['rms'],df_withwind_long.q['rms'],df_leg13_long.q['rms'],df_leg15_long.q['rms']]
    
    #%% plot figure
    
    fig, (ax0,ax1,ax2) = plt.subplots(3,1, figsize = (20,16))
    
    ## ax0 for wsp
    ax0.bar(x_locations_left, wsp_rms_cont, width, label = 'CONTROL')
    ax0.bar(x_locations_center, wsp_rms_long, width, label = 'LONGTAIL')
    #ax0.plot([0.5,5.5], [0,0], color = 'k', linestyle = 'dashed')
    ax0.set_xlim(left = 0.5, right = 5.5)
    ax0.set_xticks([1,2,3,4,5])
    ax0.set_xticklabels([])
    ax0.set_yticks([0,2,4,6,8])
    ax0.set_ylabel(r'$\delta U$, ms$^{-1}$')
    ax0.legend()
    ax0.yaxis.grid()
    
    ## ax1 for theta
    ax1.bar(x_locations_left, theta_rms_cont, width)
    ax1.bar(x_locations_center, theta_rms_long, width)
    #ax1.plot([0.5,5.5], [0,0], color = 'k', linestyle = 'dashed')
    ax1.set_xlim(left = 0.5, right = 5.5)
    ax1.set_xticks([1,2,3,4,5])
    ax1.set_xticklabels([])
    ax1.set_ylabel(r'$\delta\theta$, K')
    ax1.set_yticks([0,0.5,1,1.5,2])
    ax1.yaxis.grid()
    
     ## ax2 for q
    ax2.bar(x_locations_left, q_rms_cont, width)
    ax2.bar(x_locations_center, q_rms_long, width)
    #ax2.plot([0.5,5.5], [0,0], color = 'k', linestyle = 'dashed')
    ax2.set_xlim(left = 0.5, right = 5.5)
    ax2.set_xticks([1,2,3,4,5])
    ax2.set_xticklabels(['A','B','C','D','E'])
    ax2.set_ylabel(r'$\delta q$, gkg$^{-1}$')
    ax2.set_yticks([0,0.1,0.2,0.3,0.4,0.5,0.6])
    ax2.yaxis.grid()
    
    fig.suptitle('Case 2 - CONTROL/LONGTAIL - Root mean square bias')
    plt.tight_layout()
    
    if True:
        plt.savefig('D:/Project/Figures/PNG/306/LONGTAIL/u-cf117/stat_val/LONGTAIL_vs_CONTROL_rms_bias_bar_charts_Uthetaq.png')