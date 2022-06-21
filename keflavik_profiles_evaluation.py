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
import glm_functions as glmf

#%% oad obs data
obs_path = 'D:/Project/Obvs_Data/profiles/'
knots_conv = 0.5144
matplotlib.rcParams.update({'font.size': 24})

#%%
use_full_month = False
if use_full_month:
    filename = 'keflavik_March2018_profiles.csv'
    variable_list = ['Time','PRES','HGHT','TEMP','RELH','DRCT','SKNT','THTA','THTE']
    obs_df = pd.read_csv(f'{obs_path}{filename}').drop(0) # load data as dataframe and drop row containing units
    obs_df = obs_df[obs_df.columns & variable_list].dropna() # select variables and drop rows containing nans
    obs_df.Time = pd.to_datetime(obs_df.Time).dt.strftime('%m%dT%H') # format string times into datetimes and simplify (no year or minutes needed)
    # cast all numeric elements as floats (currently strings)
    obs_df[['PRES','HGHT','TEMP','RELH','DRCT','SKNT','THTA','THTE']] = obs_df[['PRES','HGHT','TEMP','RELH','DRCT','SKNT','THTA','THTE']].astype(float)
    # take means of soundings below 1000m (i.e. grouped by time step)
    sounding_means = obs_df[(obs_df.HGHT < 1000) & (obs_df.HGHT >500)].groupby('Time', as_index=False).agg({'PRES':'mean','TEMP':'mean','RELH':'mean', 'DRCT':'mean','SKNT':'mean','THTA':'mean','THTE':'mean'})
    # take maxs of soundings below 1000m + median of direction (i.e. grouped by time step)
    sounding_max = obs_df[(obs_df.HGHT < 1000)].groupby('Time', as_index=False).agg({'PRES':'max','TEMP':'max','RELH':'max', 'DRCT':'median','SKNT':'max','THTA':'max','THTE':'max'})
    sounding_times = np.array(obs_df.Time.drop_duplicates()) # extract times of soundings
    #print(sounding_times)
    print(sounding_means)
    
    #%% subset 12th and 19th points
    means_19th = sounding_means[(sounding_means.Time == '0319T06') | (sounding_means.Time == '0319T12') | (sounding_means.Time == '0319T18')]
    means_12th = sounding_means[(sounding_means.Time == '0312T06') | (sounding_means.Time == '0312T12') | (sounding_means.Time == '0312T18')]
    
    print(means_19th)
    #%% plot
    fig, ((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2, figsize = (16,16))
    
    ax0.set_xticks([0,90,180,270,360])
    ax0.scatter(sounding_means.DRCT,sounding_means.SKNT*knots_conv)
    ax0.scatter(means_12th.DRCT,means_12th.SKNT*knots_conv,marker = 'x', s = 400, label = '12th March')
    ax0.scatter(means_19th.DRCT,means_19th.SKNT*knots_conv,marker= 'x', s = 400, label = '19th March')
    ax0.legend()
    
    ax1.set_xticks([0,90,180,270,360])
    ax1.scatter(sounding_means.DRCT,sounding_means.THTA)
    ax1.scatter(means_12th.DRCT,means_12th.THTA,marker = 'x', s = 400, label = '12th March')
    ax1.scatter(means_19th.DRCT,means_19th.THTA,marker= 'x', s = 400, label = '19th March')
    
#%% 
use_flights = True
if use_flights:
    flight  = 306
    filename = 'keflavik_04018BIKF_profiles_upto500mb.csv'
    obs_df = pd.read_csv(f'{obs_path}{filename}')
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
        
        #%% use global 
        use_global = False
        if use_global:
            
    
            cubes = iris.load(f'D:/Project/Model_Data/{suite}/RA1M_glm_1_flt{flight}.nc')
            print(cubes)
            print(cubes[5])
            print(cubes[25])
            cubes = None
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
        ### also load 4p4km data
        nwind4p4_cube = iris.load_cube(f'{um_path}RA1M_4p4km_um_y_wind_24hrs_ph_306.nc', 'y_wind')[:,:40]
        ewind4p4_cube = iris.load_cube(f'{um_path}RA1M_4p4km_um_x_wind_24hrs_ph_306.nc', 'x_wind')[:,:40]
        magwind4p4_cube = iris.load_cube(f'{um_path}RA1M_4p4km_um_x_wind_24hrs_ph_306.nc', 'x_wind')[:,:40]
        mag_data = (nwind4p4_cube.data[:,:,:-1,:]**2 + ewind4p4_cube.data**2)**0.5 
        magwind4p4_cube.data = mag_data
        theta4p4_cube = iris.load_cube(f'{um_path}RA1M_4p4km_um_air_potential_temperature_24hrs_pi_306.nc', 'air_potential_temperature')[:,:40]
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
        theta4p4_cols = trajectory.interpolate(theta4p4_cube,sample_points, method = 'nearest')
        wind4p4_cols = trajectory.interpolate(magwind4p4_cube, sample_points, method = 'nearest')
        #%%
        seperate_panels = False
        if seperate_panels:
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
                ax0.plot(theta_coloumns.data[idx, :39], theta_coloumns.coords('altitude')[0].points[:39], label = 'RA1M 1p5km')
                ax0.plot(theta4p4_cols.data[idx,:39], theta4p4_cols.coords('altitude')[0].points[:39], label = 'RA1M 4p4km')
                ax0.set_ylim(bottom = 0,top = 5100)
                ax0.set_yticks(height_ticks)
                ax0.set_yticklabels([1,2,3,4,5])
                ax0.set_ylabel('km')
                ax0.set_xlabel(r'$\theta$, $K$')
                ax0.legend()
                #ax0.set_title('Potential temperature')
                ## ax1: wind speed
                ax1.plot(obs_speed, obs_h, label = 'Kefl. Sounding')
                ax1.plot(wind_coloumns.data[idx, :39], wind_coloumns.coords('altitude')[0].points[:39], label = 'RA1M 1p5km')
                ax1.plot(wind4p4_cols.data[idx,:39], wind4p4_cols.coords('altitude')[0].points[:39], label = 'RA1M 4p4km')
                ax1.set_ylim(bottom = 0, top = 5100)
                ax1.set_yticks(height_ticks)
                ax1.set_yticklabels([])
                ax1.set_xlabel(r'wsp, $ms^{-1}$')
                ax0.legend(fontsize = 12)
                #ax1.set_title('Windspeed')
                
                plt.tight_layout()
                #plt.savefig(f'D:/Project/Figures/PDF/{flight}/{suite}/Keflavik_profiles_{res}_{obs_time}_theta_wsp_{flight}.pdf')
                fig.suptitle(f'Keflavig-UM {res} profiles - {obs_time}')
                plt.tight_layout()
                #plt.savefig(f'D:/Project/Figures/PNG/{flight}/{suite}/Keflavik_profiles_{res}_{obs_time}_theta_wsp_{flight}.png')
                i =+ 1
        
    #%%    
        combined_panels = True
        if combined_panels:
            obs_profiles06 = obs_df[obs_df.Time == obs_times[0]]
            obs_h06 = np.array(obs_profiles06.HGHT).astype(np.float)
            obs_p06 = np.array(obs_profiles06.PRES).astype(np.float)
            obs_theta06 = np.array(obs_profiles06.THTA).astype(np.float)
            obs_speed06 =np.array(obs_profiles06.SKNT).astype(np.float)*0.5144 # convert from knots to m/s
            
            obs_profiles12 = obs_df[obs_df.Time == obs_times[1]]
            obs_h12 = np.array(obs_profiles12.HGHT).astype(np.float)
            obs_p12 = np.array(obs_profiles12.PRES).astype(np.float)
            obs_theta12 = np.array(obs_profiles12.THTA).astype(np.float)
            obs_speed12 =np.array(obs_profiles12.SKNT).astype(np.float)*0.5144 # convert from knots to m/s
            
            obs_profiles18 = obs_df[obs_df.Time == obs_times[2]]
            obs_h18 = np.array(obs_profiles18.HGHT).astype(np.float)
            obs_p18 = np.array(obs_profiles18.PRES).astype(np.float)
            obs_theta18 = np.array(obs_profiles18.THTA).astype(np.float)
            obs_speed18 =np.array(obs_profiles18.SKNT).astype(np.float)*0.5144 # convert from knots to m/s
            
            height_ticks = [1000,2000,3000,4000]#,5000]
            ### plot
            fig, ((ax0,ax1,ax2),(ax6,ax7,ax8),(ax3,ax4,ax5)) = plt.subplots(3,3, squeeze = True, figsize = (16,20))
            idx = time_idcs
            ### top row is or theta
            Ttop = 4100
            ax0.plot(obs_theta06, obs_h06, label = 'Kefl. Sounding', color = 'k') # plot obs profile
            ax0.plot(theta_coloumns.data[idx[0], :39], theta_coloumns.coords('altitude')[0].points[:39], label = 'RA1M 1p5km')
            ax0.plot(theta4p4_cols.data[idx[0], :39], theta4p4_cols.coords('altitude')[0].points[:39], label = 'RA1M 4p4km')
            ax0.set_ylim(bottom = 0,top = Ttop)
            ax0.set_yticks(height_ticks)
            ax0.set_yticklabels([1,2,3,4])#,5])
            ax0.set_ylabel('km')
            ax0.set_xlabel(r'$\theta$, $K$')
            ax0.set_xlim(left = 275, right = 302)
            ax0.legend(fontsize = 12)
            
            ax1.plot(obs_theta12, obs_h12, label = 'Kefl. Sounding', color = 'k') # plot obs profile
            ax1.plot(theta_coloumns.data[idx[1], :39], theta_coloumns.coords('altitude')[0].points[:39], label = 'RA1M 1p5km')
            ax1.plot(theta4p4_cols.data[idx[1], :39], theta4p4_cols.coords('altitude')[0].points[:39], label = 'RA1M 4p4km')
            ax1.set_ylim(bottom = 0,top = Ttop)
            ax1.set_yticks(height_ticks)
            ax1.set_yticklabels([])
            ax1.set_xlabel(r'$\theta$, $K$')
            ax1.set_xlim(left = 275, right = 302)
            
            ax2.plot(obs_theta18, obs_h18, label = 'Kefl. Sounding', color = 'k') # plot obs profile
            ax2.plot(theta_coloumns.data[idx[2], :39], theta_coloumns.coords('altitude')[0].points[:39], label = 'RA1M 1p5km')
            ax2.plot(theta4p4_cols.data[idx[2], :39], theta4p4_cols.coords('altitude')[0].points[:39], label = 'RA1M 4p4km')
            ax2.set_ylim(bottom = 0,top = Ttop)
            ax2.set_yticks(height_ticks)
            ax2.set_yticklabels([])
            ax2.set_xlabel(r'$\theta$, $K$')
            ax2.set_xlim(left = 275, right = 302)
            
            ## give top row axes sub titles with times on them
            ax0.set_title(obs_times[0])
            ax1.set_title(obs_times[1])
            ax2.set_title(obs_times[2])
            
            
            ### bottom row is wsp
            ax3.plot(obs_speed06, obs_h06, label = 'Kefl. Sounding', color = 'k')
            ax3.plot(wind_coloumns.data[idx[0], :39], wind_coloumns.coords('altitude')[0].points[:39], label = 'RA1M 1p5km')
            ax3.plot(wind4p4_cols.data[idx[0],:39], wind4p4_cols.coords('altitude')[0].points[:39], label = 'RA1M 4p4km')
            ax3.set_ylim(bottom = 0, top = Ttop)
            ax3.set_yticks(height_ticks)
            ax3.set_yticklabels([1,2,3,4])#,5])
            ax3.set_xlabel(r'wsp, $ms^{-1}$')
            ax3.set_ylabel('km')
            ax3.set_xlim(left = 0, right = 12.5)
            
            ax4.plot(obs_speed12, obs_h12, label = 'Kefl. Sounding', color = 'k')
            ax4.plot(wind_coloumns.data[idx[1], :39], wind_coloumns.coords('altitude')[0].points[:39], label = 'RA1M 1p5km')
            ax4.plot(wind4p4_cols.data[idx[1],:39], wind4p4_cols.coords('altitude')[0].points[:39], label = 'RA1M 4p4km')
            ax4.set_ylim(bottom = 0, top = Ttop)
            ax4.set_yticks(height_ticks)
            ax4.set_yticklabels([])
            ax4.set_xlabel(r'wsp, $ms^{-1}$')
            ax4.set_xlim(left = 0, right = 12.5)
            
            ax5.plot(obs_speed18, obs_h18, label = 'Kefl. Sounding', color = 'k')
            ax5.plot(wind_coloumns.data[idx[2], :39], wind_coloumns.coords('altitude')[0].points[:39], label = 'RA1M 1p5km')
            ax5.plot(wind4p4_cols.data[idx[2],:39], wind4p4_cols.coords('altitude')[0].points[:39], label = 'RA1M 4p4km')
            ax5.set_ylim(bottom = 0, top = Ttop)
            ax5.set_yticks(height_ticks)
            ax5.set_yticklabels([])
            ax5.set_xlabel(r'wsp, $ms^{-1}$')
            ax5.set_xlim(left = 0, right = 12.5)
            
            ## height cut theta profiles
            
            small_hticks = [500,1000,1500,2000]
            ax6.plot(obs_theta06, obs_h06, label = 'Kefl. Sounding', color = 'k') # plot obs profile
            ax6.plot(theta_coloumns.data[idx[0], :39], theta_coloumns.coords('altitude')[0].points[:39], label = 'RA1M 1p5km')
            ax6.plot(theta4p4_cols.data[idx[0], :39], theta4p4_cols.coords('altitude')[0].points[:39], label = 'RA1M 4p4km')
            ax6.set_ylim(bottom = 0,top = 2000)
            ax6.set_yticks(small_hticks)
            ax6.set_yticklabels([0.5,1.0,1.5,2.0])#,5])
            ax6.set_ylabel('km')
            ax6.set_xlabel(r'$\theta$, $K$')
            ax6.set_xlim(left = 275, right = 281)
            
            
            ax7.plot(obs_theta12, obs_h12, label = 'Kefl. Sounding', color = 'k') # plot obs profile
            ax7.plot(theta_coloumns.data[idx[1], :39], theta_coloumns.coords('altitude')[0].points[:39], label = 'RA1M 1p5km')
            ax7.plot(theta4p4_cols.data[idx[1], :39], theta4p4_cols.coords('altitude')[0].points[:39], label = 'RA1M 4p4km')
            ax7.set_ylim(bottom = 0,top = 2000)
            ax7.set_yticks(small_hticks)
            ax7.set_yticklabels([0.5,1.0,1.5,2.0])
            ax7.set_xlabel(r'$\theta$, $K$')
            ax7.set_xlim(left = 275, right = 281)
            
            ax8.plot(obs_theta18, obs_h18, label = 'Kefl. Sounding', color = 'k') # plot obs profile
            ax8.plot(theta_coloumns.data[idx[2], :39], theta_coloumns.coords('altitude')[0].points[:39], label = 'RA1M 1p5km')
            ax8.plot(theta4p4_cols.data[idx[2], :39], theta4p4_cols.coords('altitude')[0].points[:39], label = 'RA1M 4p4km')
            ax8.set_ylim(bottom = 0,top = 2000)
            ax8.set_yticks(small_hticks)
            ax8.set_yticklabels([0.5,1.0,1.5,2.0])
            ax8.set_xlabel(r'$\theta$, $K$')
            ax8.set_xlim(left = 275, right = 281)
            
            plt.tight_layout()
            experiment = 'CONTROL'
            #plt.savefig(f'D:/Project/Figures/PDF/{flight}/{experiment}/{suite}/Keflavik_profiles_06to18_theta_wsp_{flight}_3x3.pdf')
            plt.savefig(f'D:/Project/Figures/PNG/{flight}/{experiment}/{suite}/Keflavik_profiles_06to18_theta_wsp_{flight}_3x3.png')
           
