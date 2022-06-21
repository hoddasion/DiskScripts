# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 14:17:22 2021

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
res = '4p4km'
flight = '306'
suite = 'u-cc134'
config = 'RA1M'
alt_idx = 30
tstart_idx = 24; tend_idx = 40
raw_all_points_sr = False
raw_all_points_60s = True
seperate_leeside_points_60s = True
seperate_mountain_points_60s = True
#%%
obs_filename = f'IGP_flights_database_obs_60s_{flight}.txt'
path = 'D:/Project/Obvs_Data/Databases/'

df_obs = pd.read_csv(f'{path}{obs_filename}', delimiter = ' ')
db_coors = df_obs[['lon','lat','altgps','starttime']]
condition = (np.array(df_obs['legno'])  <5) * (np.array(df_obs['legno']) >3)
plt.scatter(np.array(db_coors['lon'])[condition], np.array(db_coors['lat'])[condition],c=  np.array(df_obs['legno'])[condition])

#%% 
selection = 3
domain = 'NML'
obspath = 'D:/Project/Obvs_Data/Databases/IGP_flights_database_obs_60s_306.txt'
umpath = f'D:/Project/Model_Data/{suite}/Matched_interpolated_values_60s.txt'
fig_name = f'{config}_{res}_boxplots_{domain}_60s_obs_flt{flight}.png'
stat_val_functions.make_boxplots(selection,res,obspath, umpath, suite, flight, fig_name, domain, config = 'RA1M', obstype = '60s',
                  ws_top = 4.2, hf_top = 250, hf_bttm = -250, save_boxplots = False)

selection = 13
domain = 'leg13'
obspath = 'D:/Project/Obvs_Data/Databases/IGP_flights_database_obs_60s_306.txt'
umpath = f'D:/Project/Model_Data/{suite}/Matched_interpolated_values_60s.txt'
fig_name = f'{config}_{res}_boxplots_{domain}_60s_obs_flt{flight}.png'
stat_val_functions.make_boxplots(selection,res,obspath, umpath, suite, flight, fig_name, domain, config = 'RA1M', obstype = '60s',
                  ws_top = 0.8, hf_top = 50, hf_bttm = -50, save_boxplots = False)

selection = 15
domain = 'leg15'
obspath = 'D:/Project/Obvs_Data/Databases/IGP_flights_database_obs_60s_306.txt'
umpath = f'D:/Project/Model_Data/{suite}/Matched_interpolated_values_60s.txt'
fig_name = f'{config}_{res}_boxplots_{domain}_60s_obs_flt{flight}.png'
stat_val_functions.make_boxplots(selection,res,obspath, umpath, suite, flight, fig_name, domain, config = 'RA1M', obstype = '60s',
                  ws_top = 0.8, hf_top = 40, hf_bttm = -40, save_boxplots = False)



#%%
if raw_all_points_sr:
    ## validate model with all obs points from short run averaged database (no distinction in legs applied)
    #%% load obs data
    obs_filename = f'IGP_flights_database_obs_sr_{flight}.txt'
    path = 'D:/Project/Obvs_Data/Databases/'
    
    df_obs = pd.read_csv(f'{path}{obs_filename}', delimiter = ' ')
    db_coors = df_obs[['lon','lat','altgps','starttime']]
    db_atmost = df_obs[['wsp','wd','w','theta','q_BUCK']]
    db_fluxes = df_obs[['sh','lh','tke','windstress']]
    
    
    df_model = pd.read_csv(f'D:/Project/Model_Data/{suite}/Matched_interpolated_values.txt')
   
    
    print(df_model)
    
    #%% plotting
    fig, ((ax0,ax1,ax2,ax3),(ax4,ax5,ax6,ax7)) = plt.subplots(2,4, figsize = (18,15))
    common_kwargs = {'widths':0.5, 'showmeans' : True,'medianprops':{'color':'red'}}
    obs_kwargs = {'labels' : ['obs']}
    mod_kwargs = {'labels' : ['UM']}
    
    ## horizontal windspeed
    ax0.boxplot(np.array(db_atmost['wsp']), positions = [0], **common_kwargs, **obs_kwargs)
    ax0.boxplot(np.array(df_model['wsp']), positions = [1], **common_kwargs, **mod_kwargs)
    ax0.set_ylabel(r'$ms^{-1}$')
    ax0.set_title('wsp')
    
    ## vertical windspeed
    ax1.boxplot(np.array(db_atmost['w']), positions = [0], **common_kwargs, **obs_kwargs)
    ax1.boxplot(np.array(df_model['w']), positions = [1], **common_kwargs, **mod_kwargs)
    ax1.set_ylabel(r'$ms^{-1}$')
    ax1.set_title('w')
    
    ## potential temperature
    ax2.boxplot(np.array(db_atmost['theta']), positions = [0], **common_kwargs, **obs_kwargs)
    ax2.boxplot(np.array(df_model['theta']), positions = [1], **common_kwargs, **mod_kwargs)
    ax2.set_ylabel(r'$K$')
    ax2.set_title(r'$\theta$')
    
    ## specific humidity
    ax3.boxplot(np.array(db_atmost['q_BUCK']), positions = [0], **common_kwargs, **obs_kwargs)
    ax3.boxplot(np.array(df_model['q'])*1000, positions = [1], **common_kwargs, **mod_kwargs)
    ax3.set_ylabel(r'$gkg^{-1}$')
    ax3.set_title('q')
    
    ## windstress
    ax4.boxplot(np.array(db_fluxes['windstress']), positions = [0], **common_kwargs, **obs_kwargs)
    ax4.boxplot(np.array(df_model['ws']), positions = [1], **common_kwargs, **mod_kwargs)
    ax4.set_ylim(top = 4.2)
    ax4.set_ylabel(r'$Nm^{-2}$')
    ax4.set_title('windstress')
    
    ## sensible heat fluxes
    ax5.boxplot(np.array(db_fluxes['sh']), positions = [0], **common_kwargs, **obs_kwargs)
    ax5.boxplot(np.array(df_model['sh']), positions = [1], **common_kwargs, **mod_kwargs)
    ax5.set_ylim(top = 160, bottom = -160)
    ax5.set_ylabel(r'$Wm^{-2}$')
    ax5.set_title('sh')
    
    ## latet heat fluxes
    ax6.boxplot(np.array(db_fluxes['lh']), positions = [0], **common_kwargs, **obs_kwargs)
    ax6.boxplot(np.array(df_model['lh']), positions = [1], **common_kwargs, **mod_kwargs)
    ax6.set_ylim(top = 160, bottom = -160)
    ax6.set_ylabel(r'$Wm^{-2}$')
    ax6.set_title('lh')
    
    ## combined heat fluxes
    ax7.boxplot(np.array(db_fluxes['lh'])+np.array(db_fluxes['sh']), positions = [0], **common_kwargs, **obs_kwargs)
    ax7.boxplot(np.array(df_model['lh'])+np.array(df_model['sh']), positions = [1], **common_kwargs, **mod_kwargs)
    ax7.set_ylim(top = 110, bottom = -110)
    ax7.set_ylabel(r'$Wm^{-2}$')
    ax7.set_title('sh + lh')
    
    fig.suptitle(f'sr obs, {res} {config} {suite} flt{flight}')
    fig.tight_layout()
    
    save_boxplots = False
    if save_boxplots:
        plt.savefig(f'D:/Project/Figures/PNG/{flight}/{suite}/stat_val/{config}_{res}_boxplots_all_obs_flt{flight}.png')
        
    #%% compute stat values and save to tables
    
    atmos_stats = db_atmost.describe()
    print(atmos_stats)
    atmos_stats.to_csv(f'D:/Project/Model_Data/{suite}/stat_val/atmostate_{res}_all_obs_stats_summary_{flight}.csv')
    flux_stats = db_fluxes.describe()
    flux_stats.to_csv(f'D:/Project/Model_Data/{suite}/stat_val/fluxes_{res}_all_obs_stats_summary_{flight}.csv')
    
    model_stats = df_model.describe()
    model_stats.to_csv(f'D:/Project/Model_Data/{suite}/stat_val/um_{res}_vars_all_legs_stats_summary_{flight}.csv')

#%%
if raw_all_points_60s:
    ## validate model with all obs points from short run averaged database (no distinction in legs applied)
    #%% load obs data
    obs_filename = f'IGP_flights_database_obs_60s_{flight}.txt'
    path = 'D:/Project/Obvs_Data/Databases/'
    
    df_obs = pd.read_csv(f'{path}{obs_filename}', delimiter = ' ')
    db_coors = df_obs[['lon','lat','altgps','starttime']]
    db_atmost = df_obs[['wsp','wd','w','theta','q_BUCK']]
    db_fluxes = df_obs[['sh','lh','tke','windstress']]
    
    flight_suffix = ''
    df_model = pd.read_csv(f'D:/Project/Model_Data/{suite}/Matched_interpolated_{res}_values_60s.txt')
   
    
    print(df_model)
    
    #%% plotting
    fig, ((ax0,ax1,ax2,ax3),(ax4,ax5,ax6,ax7)) = plt.subplots(2,4, figsize = (18,15))
    common_kwargs = {'widths':0.5, 'showmeans' : True,'medianprops':{'color':'red'}}
    obs_kwargs = {'labels' : ['obs'], 'positions' : [0]}
    mod_kwargs = {'labels' : ['UM'], 'positions' : [1]}
    
    ## horizontal windspeed
    ax0.boxplot(np.array(db_atmost['wsp']),  **common_kwargs, **obs_kwargs)
    ax0.boxplot(np.array(df_model['wsp']),**common_kwargs, **mod_kwargs)
    ax0.set_ylabel(r'$ms^{-1}$')
    ax0.set_title('wsp')
    
    ## vertical windspeed
    ax1.boxplot(np.array(db_atmost['w']), **common_kwargs, **obs_kwargs)
    ax1.boxplot(np.array(df_model['w']),  **common_kwargs, **mod_kwargs)
    
    ax1.set_ylabel(r'$ms^{-1}$')
    ax1.set_title('w')
    
    ## potential temperature
    ax2.boxplot(np.array(db_atmost['theta']),  **common_kwargs, **obs_kwargs)
    ax2.boxplot(np.array(df_model['theta']),  **common_kwargs, **mod_kwargs)
    ax2.set_ylabel(r'$K$')
    ax2.set_title(r'$\theta$')
    
    ## specific humidity
    ax3.boxplot(np.array(db_atmost['q_BUCK']),**common_kwargs, **obs_kwargs)
    ax3.boxplot(np.array(df_model['q'])*1000, **common_kwargs, **mod_kwargs)
    ax3.set_ylabel(r'$gkg^{-1}$')
    ax3.set_title('q')
    
   ## windstress
    ax4.boxplot(np.array(db_fluxes['windstress']),  **common_kwargs, **obs_kwargs)
    ax4.boxplot(np.array(df_model['ws']), **common_kwargs, **mod_kwargs)
    ax4.set_ylim(top = 4.2)
    ax4.set_ylabel(r'$Nm^{-2}$')
    ax4.set_title('windstress')
    
    ## sensible heat fluxes
    ax5.boxplot(np.array(db_fluxes['sh']), **common_kwargs, **obs_kwargs)
    ax5.boxplot(np.array(df_model['sh']), **common_kwargs, **mod_kwargs)
    ax5.set_ylim(top = 160, bottom = -160)
    ax5.set_ylabel(r'$Wm^{-2}$')
    ax5.set_title('sh')
    
    ## latet heat fluxes
    ax6.boxplot(np.array(db_fluxes['lh']),  **common_kwargs, **obs_kwargs)
    ax6.boxplot(np.array(df_model['lh']), **common_kwargs, **mod_kwargs)
    ax6.set_ylim(top = 160, bottom = -160)
    ax6.set_ylabel(r'$Wm^{-2}$')
    ax6.set_title('lh')
    
    ## combined heat fluxes
    ax7.boxplot(np.array(db_fluxes['lh'])+np.array(db_fluxes['sh']),**common_kwargs, **obs_kwargs)
    ax7.boxplot(np.array(df_model['lh'])+np.array(df_model['sh']),  **common_kwargs, **mod_kwargs)
    ax7.set_ylim(top = 110, bottom = -110)
    ax7.set_ylabel(r'$Wm^{-2}$')
    ax7.set_title('sh + lh')
    
    fig.suptitle(f'60s obs, {res} {config} {suite} flt{flight}')
    fig.tight_layout()
    
    save_boxplots = False
    if save_boxplots:
        plt.savefig(f'D:/Project/Figures/PNG/{flight}/{suite}/stat_val/{config}_{res}_boxplots_all_60s_obs_flt{flight}.png')
        
    #%% compute stat values and save to tables
    
    atmos_stats = db_atmost.describe()
    print(atmos_stats)
    atmos_stats.to_csv(f'D:/Project/Model_Data/{suite}/stat_val/atmostate_{res}_all_60s_obs_stats_summary_{flight}.csv')
    flux_stats = db_fluxes.describe()
    flux_stats.to_csv(f'D:/Project/Model_Data/{suite}/stat_val/fluxes_{res} all_60s_obs_stats_summary_{flight}.csv')
    
    model_stats = df_model.describe()
    model_stats.to_csv(f'D:/Project/Model_Data/{suite}/stat_val/um_{res}_vars_all_60s_legs_stats_summary_{flight}.csv')
    
#%%
if seperate_mountain_points_60s:
    ## validate model with seperated leeside obs points from short run averaged database 
    #%% load obs data
    obs_filename = f'IGP_flights_database_obs_60s_{flight}.txt'
    path = 'D:/Project/Obvs_Data/Databases/'
    
    df_obs = pd.read_csv(f'{path}{obs_filename}', delimiter = ' ')
    db_coors = df_obs[['lon','lat','altgps','starttime']]
    db_atmost = df_obs[['wsp','wd','w','theta','q_BUCK']]
    db_fluxes = df_obs[['sh','lh','tke','windstress']]
    db_legs = np.array(df_obs['legno'])
    condition = db_legs <=2
    
    df_model = pd.read_csv(f'D:/Project/Model_Data/{suite}/Matched_interpolated_{res}_values_60s.txt')
   
    
    print(df_model)
    
    #%% plotting
    fig, ((ax0,ax1,ax2,ax3),(ax4,ax5,ax6,ax7)) = plt.subplots(2,4, figsize = (18,15))
    common_kwargs = {'widths':0.5, 'showmeans' : True,'medianprops':{'color':'red'}}
    obs_kwargs = {'labels' : ['obs'], 'positions' : [0]}
    mod_kwargs = {'labels' : ['UM'], 'positions' : [1]}
    
    ## horizontal windspeed
    ax0.boxplot(np.array(db_atmost['wsp'])[condition],  **common_kwargs, **obs_kwargs)
    ax0.boxplot(np.array(df_model['wsp'])[condition],**common_kwargs, **mod_kwargs)
    ax0.set_ylabel(r'$ms^{-1}$')
    ax0.set_title('wsp')
    
    ## vertical windspeed
    ax1.boxplot(np.array(db_atmost['w'])[condition], **common_kwargs, **obs_kwargs)
    ax1.boxplot(np.array(df_model['w'])[condition],  **common_kwargs, **mod_kwargs)
    ax1.set_ylabel(r'$ms^{-1}$')
    ax1.set_title('w')
    
    ## potential temperature
    ax2.boxplot(np.array(db_atmost['theta'])[condition],  **common_kwargs, **obs_kwargs)
    ax2.boxplot(np.array(df_model['theta'])[condition],  **common_kwargs, **mod_kwargs)
    ax2.set_ylabel(r'$K$')
    ax2.set_title(r'$\theta$')
    
    ## specific humidity
    ax3.boxplot(np.array(db_atmost['q_BUCK'])[condition],**common_kwargs, **obs_kwargs)
    ax3.boxplot(np.array(df_model['q'])[condition]*1000, **common_kwargs, **mod_kwargs)
    ax3.set_ylabel(r'$gkg^{-1}$')
    ax3.set_title('q')
    
    ## windstress
    ax4.boxplot(np.array(db_fluxes['windstress'])[condition],  **common_kwargs, **obs_kwargs)
    ax4.boxplot(np.array(df_model['ws'])[condition], **common_kwargs, **mod_kwargs)
    ax4.set_ylim(top = 4.2)
    ax4.set_ylabel(r'$Nm^{-2}$')
    ax4.set_title('windstress')
    
    ## sensible heat fluxes
    ax5.boxplot(np.array(db_fluxes['sh'])[condition], **common_kwargs, **obs_kwargs)
    ax5.boxplot(np.array(df_model['sh'])[condition], **common_kwargs, **mod_kwargs)
    ax5.set_ylim(top = 160, bottom = -160)
    ax5.set_ylabel(r'$Wm^{-2}$')
    ax5.set_title('sh')
    
    ## latet heat fluxes
    ax6.boxplot(np.array(db_fluxes['lh'])[condition],  **common_kwargs, **obs_kwargs)
    ax6.boxplot(np.array(df_model['lh'])[condition], **common_kwargs, **mod_kwargs)
    ax6.set_ylim(top = 160, bottom = -160)
    ax6.set_ylabel(r'$Wm^{-2}$')
    ax6.set_title('lh')
    
    ## combined heat fluxes
    ax7.boxplot(np.array(db_fluxes['lh'])[condition]+np.array(db_fluxes['sh'])[condition],**common_kwargs, **obs_kwargs)
    ax7.boxplot(np.array(df_model['lh'])[condition]+np.array(df_model['sh'])[condition],  **common_kwargs, **mod_kwargs)
    ax7.set_ylim(top = 110, bottom = -110)
    ax7.set_ylabel(r'$Wm^{-2}$')
    ax7.set_title('sh + lh')
    
    fig.suptitle(f'60s mountain obs, {res} {config} {suite} flt{flight}')
    fig.tight_layout()
    
    save_boxplots = False
    if save_boxplots:
        plt.savefig(f'D:/Project/Figures/PNG/{flight}/{suite}/stat_val/{config}_{res}_boxplots_all_60s_mountain_obs_flt{flight}.png')
        
    #%% compute stat values and save to tables
    
    atmos_stats = db_atmost[condition].describe()
    print(atmos_stats)
    atmos_stats.to_csv(f'D:/Project/Model_Data/{suite}/stat_val/atmostate__{res}_all_60s_mountain_obs_stats_summary_{flight}.csv')
    flux_stats = db_fluxes[condition].describe()
    flux_stats.to_csv(f'D:/Project/Model_Data/{suite}/stat_val/fluxes_{res}_all_60s_mountain_obs_stats_summary_{flight}.csv')
    
    model_stats = df_model[condition].describe()
    model_stats.to_csv(f'D:/Project/Model_Data/{suite}/stat_val/um_{res}_vars_all_60s_mountain_legs_stats_summary_{flight}.csv')
    
if seperate_leeside_points_60s:
    ## validate model with seperated leeside obs points from short run averaged database 
    #%% load obs data
    obs_filename = f'IGP_flights_database_obs_60s_{flight}.txt'
    path = 'D:/Project/Obvs_Data/Databases/'
    
    df_obs = pd.read_csv(f'{path}{obs_filename}', delimiter = ' ')
    db_coors = df_obs[['lon','lat','altgps','starttime']]
    db_atmost = df_obs[['wsp','wd','w','theta','q_BUCK']]
    db_fluxes = df_obs[['sh','lh','tke','windstress']]
    db_legs = np.array(df_obs['legno'])
    condition = db_legs >2
    
    df_model = pd.read_csv(f'D:/Project/Model_Data/{suite}/Matched_interpolated_{res}_values_60s.txt')
   
    
    print(df_model)
    
    #%% plotting
    fig, ((ax0,ax1,ax2,ax3),(ax4,ax5,ax6,ax7)) = plt.subplots(2,4, figsize = (18,15))
    common_kwargs = {'widths':0.5, 'showmeans' : True,'medianprops':{'color':'red'}}
    obs_kwargs = {'labels' : ['obs'], 'positions' : [0]}
    mod_kwargs = {'labels' : ['UM'], 'positions' : [1]}
    
    ## horizontal windspeed
    ax0.boxplot(np.array(db_atmost['wsp'])[condition],  **common_kwargs, **obs_kwargs)
    ax0.boxplot(np.array(df_model['wsp'])[condition],**common_kwargs, **mod_kwargs)
    ax0.set_ylabel(r'$ms^{-1}$')
    ax0.set_title('wsp')
    
    ## vertical windspeed
    ax1.boxplot(np.array(db_atmost['w'])[condition], **common_kwargs, **obs_kwargs)
    ax1.boxplot(np.array(df_model['w'])[condition],  **common_kwargs, **mod_kwargs)
    ax1.set_ylabel(r'$ms^{-1}$')
    ax1.set_title('w')
    
    ## potential temperature
    ax2.boxplot(np.array(db_atmost['theta'])[condition],  **common_kwargs, **obs_kwargs)
    ax2.boxplot(np.array(df_model['theta'])[condition],  **common_kwargs, **mod_kwargs)
    ax2.set_ylabel(r'$K$')
    ax2.set_title(r'$\theta$')
    
    ## specific humidity
    ax3.boxplot(np.array(db_atmost['q_BUCK'])[condition],**common_kwargs, **obs_kwargs)
    ax3.boxplot(np.array(df_model['q'])[condition]*1000, **common_kwargs, **mod_kwargs)
    ax3.set_ylabel(r'$gkg^{-1}$')
    ax3.set_title('q')
    
    ## windstress
    ax4.boxplot(np.array(db_fluxes['windstress'])[condition],  **common_kwargs, **obs_kwargs)
    ax4.boxplot(np.array(df_model['ws'])[condition], **common_kwargs, **mod_kwargs)
    ax4.set_ylim(top = 4.2)
    ax4.set_ylabel(r'$Nm^{-2}$')
    ax4.set_title('windstress')
    
    ## sensible heat fluxes
    ax5.boxplot(np.array(db_fluxes['sh'])[condition], **common_kwargs, **obs_kwargs)
    ax5.boxplot(np.array(df_model['sh'])[condition], **common_kwargs, **mod_kwargs)
    ax5.set_ylim(top = 160, bottom = -160)
    ax5.set_ylabel(r'$Wm^{-2}$')
    ax5.set_title('sh')
    
    ## latet heat fluxes
    ax6.boxplot(np.array(db_fluxes['lh'])[condition],  **common_kwargs, **obs_kwargs)
    ax6.boxplot(np.array(df_model['lh'])[condition], **common_kwargs, **mod_kwargs)
    ax6.set_ylim(top = 160, bottom = -160)
    ax6.set_ylabel(r'$Wm^{-2}$')
    ax6.set_title('lh')
    
    ## combined heat fluxes
    ax7.boxplot(np.array(db_fluxes['lh'])[condition]+np.array(db_fluxes['sh'])[condition],**common_kwargs, **obs_kwargs)
    ax7.boxplot(np.array(df_model['lh'])[condition]+np.array(df_model['sh'])[condition],  **common_kwargs, **mod_kwargs)
    ax7.set_ylim(top = 110, bottom = -110)
    ax7.set_ylabel(r'$Wm^{-2}$')
    ax7.set_title('sh + lh')
    
    fig.suptitle(f'60s leeside obs, {res} {config} {suite} flt{flight}')
    fig.tight_layout()
    
    save_boxplots =False
    if save_boxplots:
        plt.savefig(f'D:/Project/Figures/PNG/{flight}/{suite}/stat_val/{config}_{res}_boxplots_all_60s_leeside_obs_flt{flight}.png')
        
    #%% compute stat values and save to tables
    
    atmos_stats = db_atmost[condition].describe()
    print(atmos_stats)
    atmos_stats.to_csv(f'D:/Project/Model_Data/{suite}/stat_val/atmostate_{res}_all_60s_leeside_obs_stats_summary_{flight}.csv')
    flux_stats = db_fluxes[condition].describe()
    flux_stats.to_csv(f'D:/Project/Model_Data/{suite}/stat_val/fluxes_{res}_all_60s_leeside_obs_stats_summary_{flight}.csv')
    
    model_stats = df_model[condition].describe()
    model_stats.to_csv(f'D:/Project/Model_Data/{suite}/stat_val/um_{res}_vars_all_60s_leeside_legs_stats_summary_{flight}.csv')