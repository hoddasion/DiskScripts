# -*- coding: utf-8 -*-
"""
Created on Tue Jul 27 13:57:27 2021

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

#%%



#%%
def make_boxplots(selection,res,obspath, umpath, suite, flight, fig_name, domain, config = 'RA1M', obstype = '60s', save_boxplots = True,
                  ws_top = 4.2, hf_top = 160, hf_bttm = -160):
    
    ## validate model with all obs points from short run averaged database (no distinction in legs applied)
    #%% load obs data
    
    
    
    df_obs = pd.read_csv(f'{obspath}', delimiter = ' ')
    db_coors = df_obs[['lon','lat','altgps','starttime']]
    db_atmost = df_obs[['wsp','wd','w','theta','q_BUCK']]
    db_fluxes = df_obs[['sh','lh','tke','windstress']]
    db_legno = np.array(df_obs['legno'])
    data_size = len(np.array(db_atmost['wsp']))
    df_model = pd.read_csv(f'{umpath}')
   
    condition = db_legno == -1
    for select in selection:
        temp_cond = db_legno == select
        condition = condition + temp_cond
    print(df_model)
    
    #%% plotting
    fig, ((ax0,ax1,ax2,ax3),(ax4,ax5,ax6,ax7)) = plt.subplots(2,4, figsize = (18,15))
    common_kwargs = {'widths':0.5, 'showmeans' : True,'medianprops':{'color':'red'}}
    obs_kwargs = {'labels' : ['obs']}
    mod_kwargs = {'labels' : ['UM']}
    
    ## horizontal windspeed
    ax0.boxplot(np.array(db_atmost['wsp'])[condition], positions = [0], **common_kwargs, **obs_kwargs)
    ax0.boxplot(np.array(df_model['wsp'])[condition], positions = [1], **common_kwargs, **mod_kwargs)
    ax0.set_ylabel(r'$ms^{-1}$')
    ax0.set_title('wsp')
    
    ## vertical windspeed
    ax1.boxplot(np.array(db_atmost['w'])[condition], positions = [0], **common_kwargs, **obs_kwargs)
    ax1.boxplot(np.array(df_model['w'])[condition], positions = [1], **common_kwargs, **mod_kwargs)
    ax1.set_ylabel(r'$ms^{-1}$')
    ax1.set_title('w')
    
    ## potential temperature
    ax2.boxplot(np.array(db_atmost['theta'])[condition], positions = [0], **common_kwargs, **obs_kwargs)
    ax2.boxplot(np.array(df_model['theta'])[condition], positions = [1], **common_kwargs, **mod_kwargs)
    ax2.set_ylabel(r'$K$')
    ax2.set_title(r'$\theta$')
    
    ## specific humidity
    ax3.boxplot(np.array(db_atmost['q_BUCK'])[condition], positions = [0], **common_kwargs, **obs_kwargs)
    ax3.boxplot(np.array(df_model['q'])[condition]*1000, positions = [1], **common_kwargs, **mod_kwargs)
    ax3.set_ylabel(r'$gkg^{-1}$')
    ax3.set_title('q')
    
    ## windstress
    ax4.boxplot(np.array(db_fluxes['windstress'])[condition], positions = [0], **common_kwargs, **obs_kwargs)
    ax4.boxplot(np.array(df_model['ws'])[condition], positions = [1], **common_kwargs, **mod_kwargs)
    ax4.set_ylim(top = ws_top)
    ax4.set_ylabel(r'$Nm^{-2}$')
    ax4.set_title('windstress')
    
    ## sensible heat fluxes
    ax5.boxplot(np.array(db_fluxes['sh'])[condition], positions = [0], **common_kwargs, **obs_kwargs)
    ax5.boxplot(np.array(df_model['sh'])[condition], positions = [1], **common_kwargs, **mod_kwargs)
    ax5.set_ylim(top = hf_top, bottom = hf_bttm)
    ax5.set_ylabel(r'$Wm^{-2}$')
    ax5.set_title('sh')
    
    ## latet heat fluxes
    ax6.boxplot(np.array(db_fluxes['lh'])[condition], positions = [0], **common_kwargs, **obs_kwargs)
    ax6.boxplot(np.array(df_model['lh'])[condition], positions = [1], **common_kwargs, **mod_kwargs)
    ax6.set_ylim(top = hf_top, bottom = hf_bttm)
    ax6.set_ylabel(r'$Wm^{-2}$')
    ax6.set_title('lh')
    
    ## combined heat fluxes
    ax7.boxplot(np.array(db_fluxes['lh'])[condition]+np.array(db_fluxes['sh'])[condition], positions = [0], **common_kwargs, **obs_kwargs)
    ax7.boxplot(np.array(df_model['lh'])[condition]+np.array(df_model['sh'])[condition], positions = [1], **common_kwargs, **mod_kwargs)
    ax7.set_ylim(top = hf_top, bottom = hf_bttm)
    ax7.set_ylabel(r'$Wm^{-2}$')
    ax7.set_title('sh + lh')
    
    fig.suptitle(f'{obstype} obs, {domain}, {res} {config} flt{flight}, size {data_size}')
    fig.tight_layout()
    
    
    if save_boxplots:
        plt.savefig(f'D:/Project/Figures/PNG/{flight}/{suite}/stat_val/{fig_name}')
        
    #%% compute stat values and save to tables
    
    atmos_stats = db_atmost.describe()
    print(atmos_stats)
    atmos_stats.to_csv(f'D:/Project/Model_Data/{suite}/stat_val/atmostate_{res}_all_obs_stats_summary_{flight}.csv')
    flux_stats = db_fluxes.describe()
    flux_stats.to_csv(f'D:/Project/Model_Data/{suite}/stat_val/fluxes_{res}_all_obs_stats_summary_{flight}.csv')
    
    model_stats = df_model.describe()
    model_stats.to_csv(f'D:/Project/Model_Data/{suite}/stat_val/um_{res}_vars_all_legs_stats_summary_{flight}.csv')
    
    return fig