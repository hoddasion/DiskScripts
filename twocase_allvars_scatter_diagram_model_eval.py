# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 14:49:39 2021

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

matplotlib.rcParams.update({'font.size': 24})

case1_case2_comp = False
if case1_case2_comp:
    #%% load obs data for case 2
    obs2_filename = 'IGP_flights_database_obs_60s_306.txt'
    path = 'D:/Project/Obvs_Data/Databases/'
    
    df2_obs = pd.read_csv(f'{path}{obs2_filename}', delimiter = ' ')
    
    db2_legs = np.array(df2_obs['legno'])
    cond2_mount = db2_legs <= 2
    cond2_lee   = db2_legs >  2
    
    df2_model = pd.read_csv('D:/Project/Model_Data/u-cc134/Matched_interpolated_values_60s.txt')
    #%% load obs data for case 1
    obs1_filename = 'IGP_flights_database_obs_60s_301.txt'
    path = 'D:/Project/Obvs_Data/Databases/'
    
    df1_obs = pd.read_csv(f'{path}{obs1_filename}', delimiter = ' ')
    
    db1_legs = np.array(df1_obs['legno'])
    cond1_mount = db1_legs <= 1
    cond1_lee   = (db1_legs >  1)*(db1_legs < 11)
    
    df1_model = pd.read_csv('D:/Project/Model_Data/u-bu807/Matched_interpolated_0p5km_values_60s_301.txt')
    print(df1_obs.legno)
    
    fig, ax = plt.subplots(1,1, figsize = (5,5))
    ax.plot(df1_obs.lon[cond1_mount], df1_obs.lat[cond1_mount])
    ax.plot(df1_obs.lon[cond1_lee], df1_obs.lat[cond1_lee])
    #%% make up an array for y=x line
    yx_points = [-400,400]
    #%%
    fig, ((ax8,ax9,ax10,ax11),(ax12,ax13,ax14,ax15),(ax0,ax1,ax2,ax3),(ax4,ax5,ax6,ax7)) = plt.subplots(4,4, figsize = (22,22))
    
    #### Case 2
    ### top row
    ## ax0: windspeed, wsp, case 2
    ax0.scatter(np.array(df2_obs['wsp'])[cond2_lee],np.array(df2_model['wsp'])[cond2_lee])
    ax0.scatter( np.array(df2_obs['wsp'])[cond2_mount],np.array(df2_model['wsp'])[cond2_mount])
    ax0.plot(yx_points, yx_points, color = 'k')
    ax0.set_ylim( bottom = 0,top = 25); ax0.set_xlim(left = 0,right = 25)
    ax0.set_ylabel(r'UM ms$^{-1}$')
    ax0.set_xlabel(r'obs ms$^{-1}$')
    ax0.set_title('Case 2 - wsp')
    
    ## ax1: vert windspeed, w, case 2
    ax1.scatter( np.array(df2_obs['w'])[cond2_lee],np.array(df2_model['w'])[cond2_lee])
    ax1.scatter( np.array(df2_obs['w'])[cond2_mount],np.array(df2_model['w'])[cond2_mount])
    ax1.plot(yx_points, yx_points, color = 'k')
    ax1.set_ylim( bottom = -3,top = 3); ax1.set_xlim(left = -3,right = 3)
    ax1.set_ylabel(r'UM ms$^{-1}$')
    ax1.set_xlabel(r'obs ms$^{-1}$')
    ax1.set_title('Case 2 - w')
    
    ## ax2: potential temperature, theta, case 2
    ax2.scatter( np.array(df2_obs['theta'])[cond2_lee],np.array(df2_model['theta'])[cond2_lee],)
    ax2.scatter( np.array(df2_obs['theta'])[cond2_mount],np.array(df2_model['theta'])[cond2_mount],)
    ax2.plot(yx_points, yx_points, color = 'k')
    ax2.set_ylim( bottom = 275,top = 290); ax2.set_xlim(left = 275,right = 290)
    ax2.set_ylabel('UM K')
    ax2.set_xlabel('obs K')
    ax2.set_title(r'Case 2 - $\theta$')
    
    ## ax3: specific humdity, q, case 2
    ax3.scatter( np.array(df2_obs['q_BUCK'])[cond2_lee],np.array(df2_model['q'])[cond2_lee]*1000,)
    ax3.scatter( np.array(df2_obs['q_BUCK'])[cond2_mount],np.array(df2_model['q'])[cond2_mount]*1000,)
    ax3.plot(yx_points, yx_points, color = 'k')
    ax3.set_ylim( bottom = 0,top = 4.5); ax3.set_xlim(left = 0,right = 4.5)
    ax3.set_ylabel(r'UM gkg$^{-1}$')
    ax3.set_xlabel(r'obs gkg$^{-1}$')
    ax3.set_title('Case 2 - q')
    
    ### second row
    ## ax4: windstress, ws, case 2
    ax4.scatter( np.array(df2_obs['windstress'])[cond2_lee],np.array(df2_model['ws'])[cond2_lee],)
    ax4.scatter( np.array(df2_obs['windstress'])[cond2_mount],np.array(df2_model['ws'])[cond2_mount],)
    ax4.plot(yx_points, yx_points, color = 'k')
    ax4.set_ylim( bottom = 0,top = 4.2); ax4.set_xlim(left = 0,right = 4.2)
    ax4.set_ylabel(r'UM Nm$^{-2}$')
    ax4.set_xlabel(r'obs Nm$^{-2}$')
    ax4.set_title('Case 2 - ws')
    
    ## ax5: sensible heat, sh, case 2
    ax5.scatter( np.array(df2_obs['sh'])[cond2_lee],np.array(df2_model['sh'])[cond2_lee],)
    ax5.scatter( np.array(df2_obs['sh'])[cond2_mount],np.array(df2_model['sh'])[cond2_mount],)
    ax5.plot(yx_points, yx_points, color = 'k')
    ax5.set_ylim( bottom = -160,top = 160); ax5.set_xlim(left = -160,right = 160)
    ax5.set_ylabel(r'UM Wm$^{-2}$')
    ax5.set_xlabel(r'obs Wm$^{-2}$')
    ax5.set_title('Case 2 - sh')
    
    ## ax6: latent heat, lh, case 2
    ax6.scatter( np.array(df2_obs['lh'])[cond2_lee],np.array(df2_model['lh'])[cond2_lee],)
    ax6.scatter( np.array(df2_obs['lh'])[cond2_mount],np.array(df2_model['lh'])[cond2_mount],)
    ax6.plot(yx_points, yx_points, color = 'k')
    ax6.set_ylim( bottom = -160,top = 160); ax6.set_xlim(left = -160,right = 160)
    ax6.set_ylabel(r'UM Wm$^{-2}$')
    ax6.set_xlabel(r'obs Wm$^{-2}$')
    ax6.set_title('Case 2 - lh')
    
    ## ax5: sh + lh, case 2
    ax7.scatter( np.array(df2_obs['sh'])[cond2_lee]+np.array(df2_obs['lh'])[cond2_lee],np.array(df2_model['sh'])[cond2_lee]+np.array(df2_model['lh'])[cond2_lee],)
    ax7.scatter( np.array(df2_obs['sh'])[cond2_mount]+np.array(df2_obs['lh'])[cond2_mount],np.array(df2_model['sh'])[cond2_mount]+np.array(df2_model['lh'])[cond2_mount],)
    ax7.plot(yx_points, yx_points, color = 'k')
    ax7.set_ylim( bottom = -160,top = 160); ax7.set_xlim(left = -160,right = 160)
    ax7.set_ylabel(r'UM Wm$^{-2}$')
    ax7.set_xlabel(r'obs Wm$^{-2}$')
    ax7.set_title('Case 2 - sh + lh')
    
    #### Case 1
    ### top row
    ## ax8: windspeed, wsp, case 1
    ax8.scatter( np.array(df1_obs['wsp'])[cond1_lee],np.array(df1_model['wsp'])[cond1_lee])
    ax8.scatter( np.array(df1_obs['wsp'])[cond1_mount],np.array(df1_model['wsp'])[cond1_mount])
    ax8.plot(yx_points, yx_points, color = 'k')
    ax8.set_ylim( bottom = 0,top = 25); ax8.set_xlim(left = 0,right = 25)
    ax8.set_ylabel(r'UM ms$^{-1}$')
    ax8.set_xlabel(r'obs ms$^{-1}$')
    ax8.set_title('Case 1 - wsp')
    
    ## ax9: vert windspeed, w, case 1
    ax9.scatter( np.array(df1_obs['w'])[cond1_lee],np.array(df1_model['w'])[cond1_lee])
    ax9.scatter( np.array(df1_obs['w'])[cond1_mount],np.array(df1_model['w'])[cond1_mount])
    ax9.plot(yx_points, yx_points, color = 'k')
    ax9.set_ylim( bottom = -3,top = 3); ax9.set_xlim(left = -3,right = 3)
    ax9.set_ylabel(r'UM ms$^{-1}$')
    ax9.set_xlabel(r'obs ms$^{-1}$')
    ax9.set_title('Case 1 - w')
    
    ## ax10: potential temperature, theta, case 1
    ax10.scatter( np.array(df1_obs['theta'])[cond1_lee],np.array(df1_model['theta'])[cond1_lee])
    ax10.scatter( np.array(df1_obs['theta'])[cond1_mount],np.array(df1_model['theta'])[cond1_mount],)
    ax10.plot(yx_points, yx_points, color = 'k')
    ax10.set_ylim( bottom = 270,top = 285); ax10.set_xlim(left = 270,right = 285)
    ax10.set_ylabel('UM K')
    ax10.set_xlabel('obs K')
    ax10.set_title(r'Case 1 - $\theta$')
    
    ## ax11: specific humdity, q, case 1
    ax11.scatter( np.array(df1_obs['q_BUCK'])[cond1_lee],np.array(df1_model['q'])[cond1_lee]*1000)
    ax11.scatter( np.array(df1_obs['q_BUCK'])[cond1_mount],np.array(df1_model['q'])[cond1_mount]*1000)
    ax11.plot(yx_points, yx_points, color = 'k')
    ax11.set_ylim( bottom = 0,top = 4.5); ax11.set_xlim(left = 0,right = 4.5)
    ax11.set_ylabel(r'UM gkg$^{-1}$')
    ax11.set_xlabel(r'obs gkg$^{-1}$')
    ax11.set_title('Case 1 - q')
    
    ### second row
    ## ax12: windstress, ws, case 1
    ax12.scatter( np.array(df1_obs['windstress'])[cond1_lee],np.array(df1_model['ws'])[cond1_lee])
    ax12.scatter( np.array(df1_obs['windstress'])[cond1_mount],np.array(df1_model['ws'])[cond1_mount])
    ax12.plot(yx_points, yx_points, color = 'k')
    ax12.set_ylim( bottom = 0,top = 2.5); ax12.set_xlim(left = 0,right = 2.5)
    ax12.set_ylabel(r'UM Nm$^{-2}$')
    ax12.set_xlabel(r'obs Nm$^{-2}$')
    ax12.set_title('Case 1 - ws')
    
    ## ax13: sensible heat, sh, case 1
    ax13.scatter( np.array(df1_obs['sh'])[cond1_lee],np.array(df1_model['sh'])[cond1_lee])
    ax13.scatter( np.array(df1_obs['sh'])[cond1_mount],np.array(df1_model['sh'])[cond1_mount])
    ax13.plot(yx_points, yx_points, color = 'k')
    ax13.set_ylim( bottom = -250,top = 250); ax13.set_xlim(left = -250,right = 250)
    ax13.set_ylabel(r'UM Wm$^{-2}$')
    ax13.set_xlabel(r'obs Wm$^{-2}$')
    ax13.set_title('Case 1 - sh')
    
    ## ax14: latent heat, lh, case 1
    ax14.scatter( np.array(df1_obs['lh'])[cond1_lee],np.array(df1_model['lh'])[cond1_lee])
    ax14.scatter( np.array(df1_obs['lh'])[cond1_mount],np.array(df1_model['lh'])[cond1_mount])
    ax14.plot(yx_points, yx_points, color = 'k')
    ax14.set_ylim( bottom = -250,top = 250); ax14.set_xlim(left = -250,right = 250)
    ax14.set_ylabel(r'UM Wm$^{-2}$')
    ax14.set_xlabel(r'obs Wm$^{-2}$')
    ax14.set_title('Case 1 - lh')
    
    ## ax15: sh + lh, case 1
    ax15.scatter( np.array(df1_obs['sh'])[cond1_lee]+np.array(df1_obs['lh'])[cond1_lee],np.array(df1_model['sh'])[cond1_lee]+np.array(df1_model['lh'])[cond1_lee])
    ax15.scatter( np.array(df1_obs['sh'])[cond1_mount]+np.array(df1_obs['lh'])[cond1_mount],np.array(df1_model['sh'])[cond1_mount]+np.array(df1_model['lh'])[cond1_mount])
    ax15.plot(yx_points, yx_points, color = 'k')
    ax15.set_ylim( bottom = -250,top = 250); ax15.set_xlim(left = -250,right = 250)
    ax15.set_ylabel(r'UM Wm$^{-2}$')
    ax15.set_xlabel(r'obs Wm$^{-2}$')
    ax15.set_title('Case 1 - sh + lh')
    
    fig.tight_layout()
    plt.savefig('../Figures/Thesis/all_vars_case1bu807_case2cc134_scatter_UMobs.png')

#%%
inter_experiment = True
if inter_experiment:
    config1 = 'LONGTAIL'
    experiment1 = 'LONGTAIL'
    suite1 = 'u-cf117'
    config2 = 'RA1M'
    experiment2 = 'CONTROL'
    flight = 306
    #%% load obs data for case 2
    obs_filename = f'IGP_flights_database_obs_60s_{flight}.txt'
    path = 'D:/Project/Obvs_Data/Databases/'
    
    df_obs = pd.read_csv(f'{path}{obs_filename}', delimiter = ' ')
    
    db_legs = np.array(df_obs['legno'])
    cond_mount = db_legs <= 2
    cond_lee   = db_legs >  2
    cond1_mount = cond_mount
    cond1_lee = cond_lee
    df1_obs = df_obs
    df2_model = pd.read_csv('D:/Project/Model_Data/u-cc134/Matched_interpolated_values_60s.txt')
    #%% load obs data for case 1
    
    
    df1_model = pd.read_csv(f'D:/Project/Model_Data/{suite1}/{config1}_Matched_interpolated_values_0p5km_60s_{flight}.txt')
    print(df_obs.legno)
    
    
    #%% make up an array for y=x line
    yx_points = [-400,400]
    
    #%%
    fig, ((ax8,ax9,ax10,ax11),(ax12,ax13,ax14,ax15),(ax0,ax1,ax2,ax3),(ax4,ax5,ax6,ax7)) = plt.subplots(4,4, figsize = (22,22))
    
    #### Case 2
    ### top row
    ## ax0: windspeed, wsp, case 2
    ax0.scatter(np.array(df_obs['wsp'])[cond_lee],np.array(df2_model['wsp'])[cond_lee])
    ax0.scatter( np.array(df_obs['wsp'])[cond_mount],np.array(df2_model['wsp'])[cond_mount])
    ax0.plot(yx_points, yx_points, color = 'k')
    ax0.set_ylim( bottom = 0,top = 25); ax0.set_xlim(left = 0,right = 25)
    ax0.set_ylabel(r'UM ms$^{-1}$')
    ax0.set_xlabel(r'obs ms$^{-1}$')
    ax0.set_title('Control - wsp')
    
    ## ax1: vert windspeed, w, case 2
    ax1.scatter( np.array(df_obs['w'])[cond_lee],np.array(df2_model['w'])[cond_lee])
    ax1.scatter( np.array(df_obs['w'])[cond_mount],np.array(df2_model['w'])[cond_mount])
    ax1.plot(yx_points, yx_points, color = 'k')
    ax1.set_ylim( bottom = -3,top = 3); ax1.set_xlim(left = -3,right = 3)
    ax1.set_ylabel(r'UM ms$^{-1}$')
    ax1.set_xlabel(r'obs ms$^{-1}$')
    ax1.set_title('Control - w')
    
    ## ax2: potential temperature, theta, case 2
    ax2.scatter( np.array(df_obs['theta'])[cond_lee],np.array(df2_model['theta'])[cond_lee],)
    ax2.scatter( np.array(df_obs['theta'])[cond_mount],np.array(df2_model['theta'])[cond_mount],)
    ax2.plot(yx_points, yx_points, color = 'k')
    ax2.set_ylim( bottom = 275,top = 290); ax2.set_xlim(left = 275,right = 290)
    ax2.set_ylabel('UM K')
    ax2.set_xlabel('obs K')
    ax2.set_title(r'Control - $\theta$')
    
    ## ax3: specific humdity, q, case 2
    ax3.scatter( np.array(df_obs['q_BUCK'])[cond_lee],np.array(df2_model['q'])[cond_lee]*1000,)
    ax3.scatter( np.array(df_obs['q_BUCK'])[cond_mount],np.array(df2_model['q'])[cond_mount]*1000,)
    ax3.plot(yx_points, yx_points, color = 'k')
    ax3.set_ylim( bottom = 0,top = 4.5); ax3.set_xlim(left = 0,right = 4.5)
    ax3.set_ylabel(r'UM gkg$^{-1}$')
    ax3.set_xlabel(r'obs gkg$^{-1}$')
    ax3.set_title('Control - q')
    
    ### second row
    ## ax4: windstress, ws, case 2
    ax4.scatter( np.array(df_obs['windstress'])[cond_lee],np.array(df2_model['ws'])[cond_lee],)
    ax4.scatter( np.array(df_obs['windstress'])[cond_mount],np.array(df2_model['ws'])[cond_mount],)
    ax4.plot(yx_points, yx_points, color = 'k')
    ax4.set_ylim( bottom = 0,top = 4.2); ax4.set_xlim(left = 0,right = 4.2)
    ax4.set_ylabel(r'UM Nm$^{-2}$')
    ax4.set_xlabel(r'obs Nm$^{-2}$')
    ax4.set_title('Control - ws')
    
    ## ax5: sensible heat, sh, case 2
    ax5.scatter( np.array(df_obs['sh'])[cond_lee],np.array(df2_model['sh'])[cond_lee],)
    ax5.scatter( np.array(df_obs['sh'])[cond_mount],np.array(df2_model['sh'])[cond_mount],)
    ax5.plot(yx_points, yx_points, color = 'k')
    ax5.set_ylim( bottom = -250,top = 250); ax5.set_xlim(left = -250,right = 250)
    ax5.set_ylabel(r'UM Wm$^{-2}$')
    ax5.set_xlabel(r'obs Wm$^{-2}$')
    ax5.set_title('Control - sh')
    
    ## ax6: latent heat, lh, case 2
    ax6.scatter( np.array(df_obs['lh'])[cond_lee],np.array(df2_model['lh'])[cond_lee],)
    ax6.scatter( np.array(df_obs['lh'])[cond_mount],np.array(df2_model['lh'])[cond_mount],)
    ax6.plot(yx_points, yx_points, color = 'k')
    ax6.set_ylim( bottom = -250,top = 250); ax6.set_xlim(left = -250,right = 250)
    ax6.set_ylabel(r'UM Wm$^{-2}$')
    ax6.set_xlabel(r'obs Wm$^{-2}$')
    ax6.set_title('Control - lh')
    
    ## ax5: sh + lh, case 2
    ax7.scatter( np.array(df_obs['sh'])[cond_lee]+np.array(df_obs['lh'])[cond_lee],np.array(df2_model['sh'])[cond_lee]+np.array(df2_model['lh'])[cond_lee],)
    ax7.scatter( np.array(df_obs['sh'])[cond_mount]+np.array(df_obs['lh'])[cond_mount],np.array(df2_model['sh'])[cond_mount]+np.array(df2_model['lh'])[cond_mount],)
    ax7.plot(yx_points, yx_points, color = 'k')
    ax7.set_ylim( bottom = -250,top = 250); ax7.set_xlim(left = -250,right = 250)
    ax7.set_ylabel(r'UM Wm$^{-2}$')
    ax7.set_xlabel(r'obs Wm$^{-2}$')
    ax7.set_title('Control - sh + lh')
    
    #### Case 1
    ### top row
    ## ax8: windspeed, wsp, case 1
    ax8.scatter( np.array(df1_obs['wsp'])[cond1_lee],np.array(df1_model['wsp'])[cond1_lee])
    ax8.scatter( np.array(df1_obs['wsp'])[cond1_mount],np.array(df1_model['wsp'])[cond1_mount])
    ax8.plot(yx_points, yx_points, color = 'k')
    ax8.set_ylim( bottom = 0,top = 25); ax8.set_xlim(left = 0,right = 25)
    ax8.set_ylabel(r'UM ms$^{-1}$')
    ax8.set_xlabel(r'obs ms$^{-1}$')
    ax8.set_title('Longtail - wsp')
    
    ## ax9: vert windspeed, w, case 1
    ax9.scatter( np.array(df1_obs['w'])[cond1_lee],np.array(df1_model['w'])[cond1_lee])
    ax9.scatter( np.array(df1_obs['w'])[cond1_mount],np.array(df1_model['w'])[cond1_mount])
    ax9.plot(yx_points, yx_points, color = 'k')
    ax9.set_ylim( bottom = -3,top = 3); ax9.set_xlim(left = -3,right = 3)
    ax9.set_ylabel(r'UM ms$^{-1}$')
    ax9.set_xlabel(r'obs ms$^{-1}$')
    ax9.set_title('Longtail - w')
    
    ## ax10: potential temperature, theta, case 1
    ax10.scatter( np.array(df1_obs['theta'])[cond1_lee],np.array(df1_model['theta'])[cond1_lee])
    ax10.scatter( np.array(df1_obs['theta'])[cond1_mount],np.array(df1_model['theta'])[cond1_mount],)
    ax10.plot(yx_points, yx_points, color = 'k')
    ax10.set_ylim( bottom = 275,top = 290); ax10.set_xlim(left = 275,right = 290)
    ax10.set_ylabel('UM K')
    ax10.set_xlabel('obs K')
    ax10.set_title(r'Longtail- $\theta$')
    
    ## ax11: specific humdity, q, case 1
    ax11.scatter( np.array(df1_obs['q_BUCK'])[cond1_lee],np.array(df1_model['q'])[cond1_lee]*1000)
    ax11.scatter( np.array(df1_obs['q_BUCK'])[cond1_mount],np.array(df1_model['q'])[cond1_mount]*1000)
    ax11.plot(yx_points, yx_points, color = 'k')
    ax11.set_ylim( bottom = 0,top = 4.5); ax11.set_xlim(left = 0,right = 4.5)
    ax11.set_ylabel(r'UM gkg$^{-1}$')
    ax11.set_xlabel(r'obs gkg$^{-1}$')
    ax11.set_title('Longtail - q')
    
    ### second row
    ## ax12: windstress, ws, case 1
    ax12.scatter( np.array(df1_obs['windstress'])[cond1_lee],np.array(df1_model['ws'])[cond1_lee])
    ax12.scatter( np.array(df1_obs['windstress'])[cond1_mount],np.array(df1_model['ws'])[cond1_mount])
    ax12.plot(yx_points, yx_points, color = 'k')
    ax12.set_ylim( bottom = 0,top = 4.2); ax12.set_xlim(left = 0,right =4.2)
    ax12.set_ylabel(r'UM Nm$^{-2}$')
    ax12.set_xlabel(r'obs Nm$^{-2}$')
    ax12.set_title('Longtail - ws')
    
    ## ax13: sensible heat, sh, case 1
    ax13.scatter( np.array(df1_obs['sh'])[cond1_lee],np.array(df1_model['sh'])[cond1_lee])
    ax13.scatter( np.array(df1_obs['sh'])[cond1_mount],np.array(df1_model['sh'])[cond1_mount])
    ax13.plot(yx_points, yx_points, color = 'k')
    ax13.set_ylim( bottom = -250,top = 250); ax13.set_xlim(left = -250,right = 250)
    ax13.set_ylabel(r'UM Wm$^{-2}$')
    ax13.set_xlabel(r'obs Wm$^{-2}$')
    ax13.set_title('Longtail - sh')
    
    ## ax14: latent heat, lh, case 1
    ax14.scatter( np.array(df1_obs['lh'])[cond1_lee],np.array(df1_model['lh'])[cond1_lee])
    ax14.scatter( np.array(df1_obs['lh'])[cond1_mount],np.array(df1_model['lh'])[cond1_mount])
    ax14.plot(yx_points, yx_points, color = 'k')
    ax14.set_ylim( bottom = -250,top = 250); ax14.set_xlim(left = -250,right = 250)
    ax14.set_ylabel(r'UM Wm$^{-2}$')
    ax14.set_xlabel(r'obs Wm$^{-2}$')
    ax14.set_title('Longtail - lh')
    
    ## ax15: sh + lh, case 1
    ax15.scatter( np.array(df1_obs['sh'])[cond1_lee]+np.array(df1_obs['lh'])[cond1_lee],np.array(df1_model['sh'])[cond1_lee]+np.array(df1_model['lh'])[cond1_lee])
    ax15.scatter( np.array(df1_obs['sh'])[cond1_mount]+np.array(df1_obs['lh'])[cond1_mount],np.array(df1_model['sh'])[cond1_mount]+np.array(df1_model['lh'])[cond1_mount])
    ax15.plot(yx_points, yx_points, color = 'k')
    ax15.set_ylim( bottom = -250,top = 250); ax15.set_xlim(left = -250,right = 250)
    ax15.set_ylabel(r'UM Wm$^{-2}$')
    ax15.set_xlabel(r'obs Wm$^{-2}$')
    ax15.set_title('Longtail - sh + lh')
    
    fig.tight_layout()
    plt.savefig('../Figures/Thesis/all_vars_case2_LONGTAIL_CONTROL_scatter_UMobs.png')
#%%  
eight_panel_interexperiment = False
if eight_panel_interexperiment:
    config1 = 'LONGTAIL'
    experiment1 = 'LONGTAIL'
    suite1 = 'u-cf117'
    config2 = 'RA1M'
    experiment2 = 'CONTROL'
    flight = 306
    #%% load obs data for case 2
    obs_filename = f'IGP_flights_database_obs_60s_{flight}.txt'
    path = 'D:/Project/Obvs_Data/Databases/'
    
    df_obs = pd.read_csv(f'{path}{obs_filename}', delimiter = ' ')
    
    db_legs = np.array(df_obs['legno'])
    cond_mount = db_legs <= 2
    cond_lee   = db_legs >  2
    cond1_mount = cond_mount
    cond1_lee = cond_lee
    df1_obs = df_obs
    df2_model = pd.read_csv('D:/Project/Model_Data/u-cc134/Matched_interpolated_values_60s.txt')
    #%% load obs data for case 1
    
    
    df1_model = pd.read_csv(f'D:/Project/Model_Data/{suite1}/{config1}_Matched_interpolated_values_0p5km_60s_{flight}.txt')
    print(df_obs.legno)
    
    fig, ax = plt.subplots(1,1, figsize = (5,5))
    ax.plot(df_obs.lon[cond_mount], df_obs.lat[cond_mount])
    ax.plot(df_obs.lon[cond_lee], df_obs.lat[cond_lee])
    #%% make up an array for y=x line
    yx_points = [-400,400]
    
    #%%
    fig, ((ax0,ax1,ax2,ax3),(ax4,ax5,ax6,ax7)) = plt.subplots(2,4, figsize = (22,11))
    
    #### Case 2
    ### top row
    ## ax0: windspeed, wsp, case 2
    ax0.scatter(np.array(df2_model['wsp'])[cond_lee],np.array(df1_model['wsp'])[cond_lee])
    ax0.scatter(np.array(df2_model['wsp'])[cond_mount],np.array(df1_model['wsp'])[cond_mount])
    ax0.plot(yx_points, yx_points, color = 'k')
    ax0.set_ylim( bottom = 0,top = 25); ax0.set_xlim(left = 0,right = 25)
    ax0.set_ylabel(r'LONG ms$^{-1}$')
    ax0.set_xlabel(r'CON ms$^{-1}$')
    ax0.set_title('wsp')
    
    ## ax1: vert windspeed, w, case 2
    ax1.scatter(np.array(df2_model['w'])[cond_lee],np.array(df1_model['w'])[cond_lee])
    ax1.scatter(np.array(df2_model['w'])[cond_mount],np.array(df1_model['w'])[cond_mount])
    ax1.plot(yx_points, yx_points, color = 'k')
    ax1.set_ylim( bottom = -3,top = 3); ax1.set_xlim(left = -3,right = 3)
    ax1.set_ylabel(r'LONG ms$^{-1}$')
    ax1.set_xlabel(r'CON ms$^{-1}$')
    ax1.set_title('w')
    
    ## ax2: potential temperature, theta, case 2
    ax2.scatter(np.array(df2_model['theta'])[cond_lee],np.array(df1_model['theta'])[cond_lee])
    ax2.scatter(np.array(df2_model['theta'])[cond_mount],np.array(df1_model['theta'])[cond_mount])
    ax2.plot(yx_points, yx_points, color = 'k')
    ax2.set_ylim( bottom = 275,top = 290); ax2.set_xlim(left = 275,right = 290)
    ax2.set_ylabel('LONG K')
    ax2.set_xlabel('CON K')
    ax2.set_title(r'$\theta$')
    
    ## ax3: specific humdity, q, case 2
    ax3.scatter(np.array(df2_model['q'])[cond_lee]*1000,np.array(df1_model['q'])[cond_lee]*1000)
    ax3.scatter(np.array(df2_model['q'])[cond_mount]*1000,np.array(df1_model['q'])[cond_mount]*1000)
    ax3.plot(yx_points, yx_points, color = 'k')
    ax3.set_ylim( bottom = 0,top = 4.5); ax3.set_xlim(left = 0,right = 4.5)
    ax3.set_ylabel(r'LONG gkg$^{-1}$')
    ax3.set_xlabel(r'CON gkg$^{-1}$')
    ax3.set_title('q')
    
    ### second row
    ## ax4: windstress, ws, case 2
    ax4.scatter(np.array(df2_model['ws'])[cond_lee],np.array(df1_model['ws'])[cond_lee])
    ax4.scatter(np.array(df2_model['ws'])[cond_mount],np.array(df1_model['ws'])[cond_mount])
    ax4.plot(yx_points, yx_points, color = 'k')
    ax4.set_ylim( bottom = 0,top = 4.2); ax4.set_xlim(left = 0,right = 4.2)
    ax4.set_ylabel(r'LONG Nm$^{-2}$')
    ax4.set_xlabel(r'CON Nm$^{-2}$')
    ax4.set_title('ws')
    
    ## ax5: sensible heat, sh, case 2
    ax5.scatter(np.array(df2_model['sh'])[cond_lee],np.array(df1_model['sh'])[cond_lee])
    ax5.scatter(np.array(df2_model['sh'])[cond_mount],np.array(df1_model['sh'])[cond_mount])
    ax5.plot(yx_points, yx_points, color = 'k')
    ax5.set_ylim( bottom = -160,top = 160); ax5.set_xlim(left = -160,right = 160)
    ax5.set_ylabel(r'LONG Wm$^{-2}$')
    ax5.set_xlabel(r'CON Wm$^{-2}$')
    ax5.set_title('sh')
    
    ## ax6: latent heat, lh, case 2
    ax6.scatter(np.array(df2_model['lh'])[cond_lee],np.array(df1_model['lh'])[cond_lee])
    ax6.scatter(np.array(df2_model['lh'])[cond_mount],np.array(df1_model['lh'])[cond_mount])
    ax6.plot(yx_points, yx_points, color = 'k')
    ax6.set_ylim( bottom = -160,top = 160); ax6.set_xlim(left = -160,right = 160)
    ax6.set_ylabel(r'LONG Wm$^{-2}$')
    ax6.set_xlabel(r'CON Wm$^{-2}$')
    ax6.set_title('lh')
    
    ## ax5: sh + lh, case 2
    ax7.scatter(np.array(df2_model['sh'])[cond_lee]+np.array(df2_model['lh'])[cond_lee],np.array(df1_model['sh'])[cond_lee]+np.array(df1_model['lh'])[cond_lee])
    ax7.scatter(np.array(df2_model['sh'])[cond_mount]+np.array(df2_model['lh'])[cond_mount],np.array(df1_model['sh'])[cond_mount]+np.array(df1_model['lh'])[cond_mount])
    ax7.plot(yx_points, yx_points, color = 'k')
    ax7.set_ylim( bottom = -160,top = 160); ax7.set_xlim(left = -160,right = 160)
    ax7.set_ylabel(r'LONG Wm$^{-2}$')
    ax7.set_xlabel(r'CON Wm$^{-2}$')
    ax7.set_title('sh + lh')
    
    fig.tight_layout()
    plt.savefig('../Figures/Thesis/all_vars_case2_LONGTAIL_CONTROL_scatter_8panel.png')