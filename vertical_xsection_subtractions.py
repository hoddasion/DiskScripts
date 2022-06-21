# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 13:25:00 2022

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
import matplotlib.colors as colors
import warnings
from datetime import datetime


warnings.filterwarnings("ignore")
#%%
plt.rcParams['font.size'] = 28
resolutions = ['0p5km']#,'1p5km', '4p4km']


lat_lon_ord = ('grid_latitude', 'grid_longitude')
lon_lat_ord = ('grid_longitude', 'grid_latitude')


alt_idx = 35
long = ''



#%% initiate resolution and flight leg  loops
for res in resolutions:
    for leg in [1,4,7,13,15]:
        ## load control data
        ## load available atmostate vertsclice cubes 
        factor_suffix = '_f2'
        control_xwind_cube = iris.load_cube(f'D:/Project/Model_Data/u-cc134/Vertical/x_wind/RA1M_vertslice_{res}_x_wind_flt306_leg{leg}{factor_suffix}.nc', 'x_wind')[:,:alt_idx]
        control_ywind_cube = iris.load_cube(f'D:/Project/Model_Data/u-cc134/Vertical/y_wind/RA1M_vertslice_{res}_y_wind_flt306_leg{leg}{factor_suffix}.nc', 'y_wind')[:,:alt_idx]
        control_theta_cube = iris.load_cube(f'D:/Project/Model_Data/u-cc134/Vertical/air_potential_temperature/RA1M_vertslice_{res}_air_potential_temperature_flt306_leg{leg}{factor_suffix}.nc', 'air_potential_temperature')[:,:alt_idx]
        control_q_cube = iris.load_cube(f'D:/Project/Model_Data/u-cc134/Vertical/specific_humidity/RA1M_vertslice_{res}_specific_humidity_flt306_leg{leg}{factor_suffix}.nc', 'specific_humidity')[:,:alt_idx]
        control_w_cube = iris.load_cube(f'D:/Project/Model_Data/u-cc134/Vertical/upward_air_velocity/RA1M_vertslice_{res}_upward_air_velocity_flt306_leg{leg}{factor_suffix}.nc', 'upward_air_velocity')[:,:alt_idx]
        ## compute wind magnitude
        wind_data = (control_xwind_cube.data**2 + control_ywind_cube.data**2)**0.5
        control_wind_cube = iris.load_cube(f'D:/Project/Model_Data/u-cc134/Vertical/x_wind/RA1M_vertslice_{res}_x_wind_flt306_leg{leg}{factor_suffix}.nc', 'x_wind')[:,:alt_idx]
        control_wind_cube.data = wind_data
        
        ## load available flux vertslice cubes
        control_sh_cube = iris.load_cube(f'D:/Project/Model_Data/u-cc134/Vertical/upward_heat_flux_in_air/RA1M_vertslice_{res}_upward_heat_flux_in_air_flt306_leg{leg}{factor_suffix}.nc', 'upward_heat_flux_in_air')[:,:alt_idx]
        control_lh_cube = iris.load_cube(f'D:/Project/Model_Data/u-cc134/Vertical/upward_latent_heat_flux_in_air/RA1M_vertslice_{res}_upward_latent_heat_flux_in_air_flt306_leg{leg}{factor_suffix}.nc', 'upward_latent_heat_flux_in_air')[:,:alt_idx]
        control_wse_cube = iris.load_cube(f'D:/Project/Model_Data/u-cc134/Vertical/atmosphere_downward_eastward_stress/RA1M_vertslice_{res}_atmosphere_downward_eastward_stress_flt306_leg{leg}{factor_suffix}.nc', 'atmosphere_downward_eastward_stress')[:,:alt_idx]
        control_wsn_cube = iris.load_cube(f'D:/Project/Model_Data/u-cc134/Vertical/atmosphere_downward_northward_stress/RA1M_vertslice_{res}_atmosphere_downward_northward_stress_flt306_leg{leg}{factor_suffix}.nc', 'atmosphere_downward_northward_stress')[:,:alt_idx]
        ## compute stress magnitude
        control_ws_cube = iris.load_cube(f'D:/Project/Model_Data/u-cc134/Vertical/atmosphere_downward_eastward_stress/RA1M_vertslice_{res}_atmosphere_downward_eastward_stress_flt306_leg{leg}{factor_suffix}.nc', 'atmosphere_downward_eastward_stress')[:,:alt_idx]
        ws_data = (control_wse_cube.data**2 + control_wsn_cube.data**2)**0.5
        control_ws_cube.data = ws_data
        
        print('Control data loaded')
        
        ## load experiment data
        experiment = 'LONGTAIL'
        suite = 'u-cf117'
        config = 'LONGTAIL'
        factor_suffix = ''
        
        test_xwind_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/x_wind/{config}_vertslice_{res}_x_wind_flt306_leg{leg}{factor_suffix}.nc', 'x_wind')[:,:alt_idx]
        test_ywind_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/y_wind/{config}_vertslice_{res}_y_wind_flt306_leg{leg}{factor_suffix}.nc', 'y_wind')[:,:alt_idx]
        test_theta_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/air_potential_temperature/{config}_vertslice_{res}_air_potential_temperature_flt306_leg{leg}{factor_suffix}.nc', 'air_potential_temperature')[:,:alt_idx]
        test_q_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/specific_humidity/{config}_vertslice_{res}_specific_humidity_flt306_leg{leg}{factor_suffix}.nc', 'specific_humidity')[:,:alt_idx]
        test_w_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/upward_air_velocity/{config}_vertslice_{res}_upward_air_velocity_flt306_leg{leg}{factor_suffix}.nc', 'upward_air_velocity')[:,:alt_idx]
        ## compute wind magnitude
        wind_data = (test_xwind_cube.data**2 + test_ywind_cube.data**2)**0.5
        test_wind_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/x_wind/{config}_vertslice_{res}_x_wind_flt306_leg{leg}{factor_suffix}.nc', 'x_wind')[:,:alt_idx]
        test_wind_cube.data = wind_data
        ## load available flux vertslice cubes
        test_sh_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/upward_heat_flux_in_air/{config}_vertslice_{res}_upward_heat_flux_in_air_flt306_leg{leg}{factor_suffix}.nc', 'upward_heat_flux_in_air')[:,:alt_idx]
        test_lh_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/upward_latent_heat_flux_in_air/{config}_vertslice_{res}_upward_latent_heat_flux_in_air_flt306_leg{leg}{factor_suffix}.nc', 'upward_latent_heat_flux_in_air')[:,:alt_idx]
        test_wse_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/atmosphere_downward_eastward_stress/{config}_vertslice_{res}_atmosphere_downward_eastward_stress_flt306_leg{leg}{factor_suffix}.nc', 'atmosphere_downward_eastward_stress')[:,:alt_idx]
        test_wsn_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/atmosphere_downward_northward_stress/{config}_vertslice_{res}_atmosphere_downward_northward_stress_flt306_leg{leg}{factor_suffix}.nc', 'atmosphere_downward_northward_stress')[:,:alt_idx]
        ## compute stress magnitude
        test_ws_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/atmosphere_downward_eastward_stress/{config}_vertslice_{res}_atmosphere_downward_eastward_stress_flt306_leg{leg}{factor_suffix}.nc', 'atmosphere_downward_eastward_stress')[:,:alt_idx]
        ws_data = (test_wse_cube.data**2 + test_wsn_cube.data**2)**0.5
        test_ws_cube.data = ws_data
        print('Test data loaded')
        
        #%% compute subtractions
        
        diff_wind_cube = test_wind_cube - control_wind_cube
        diff_theta_cube = test_theta_cube - control_theta_cube
        diff_q_cube = test_q_cube - control_q_cube
        diff_w_cube = test_w_cube - control_w_cube
        
        diff_ws_cube = test_ws_cube - control_ws_cube
        diff_sh_cube = test_sh_cube - control_sh_cube
        diff_lh_cube = test_lh_cube - control_lh_cube
        
        print('Differences computed')
        #%% format date time for file name and figure title
        time = diff_w_cube.coord('time')
        dates = time.units.num2date(time.points)
        
        #%% get scalar distance ticklist for x axis labelling
        tick_coors, tick_dists = workshop.dist_tick_generator(diff_w_cube, 'grid_latitude', ticknumber = 6, decprecision = 2)
        #%% loop through time steps
        mass_plots = False
        if mass_plots:
            for t in range(48):
                print('t',t)
                #%% plot these cubes
                if t > -1: # for now single out single time
                    
                    #%% set colorbar normalisations
                    uni_cmap = 'bwr'             
                    
                    ## atmostate
                    wspmin = -10; wspmax = -wspmin
                    wmin = -3; wmax = -wmin
                    Tmin = -3; Tmax = -Tmin
                    qmin = -2; qmax = -qmin
                    q_levels = np.arange(0,50)*0.5
                    wsp_norm = matplotlib.colors.Normalize(vmin = wspmin, vmax = wspmax)
                    w_norm = matplotlib.colors.Normalize(vmin = wmin, vmax = wmax)
                    theta_norm = matplotlib.colors.Normalize(vmin = Tmin, vmax = Tmax)
                    q_norm = matplotlib.colors.Normalize(vmin = qmin, vmax = qmax)
                    
                    wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = uni_cmap)
                    w_map = matplotlib.cm.ScalarMappable(norm = w_norm, cmap = uni_cmap)
                    theta_map = matplotlib.cm.ScalarMappable(norm = theta_norm, cmap = uni_cmap)
                    q_map = matplotlib.cm.ScalarMappable(norm = q_norm, cmap = uni_cmap)
                    
                    ## fluxes
                    wsmin = -3; wsmax = -wsmin
                    shmin = -400; shmax = -shmin
                    lhmin = -400; lhmax = -lhmin
                    
                    ws_norm = matplotlib.colors.Normalize(vmin = wsmin, vmax = wsmax)
                    sh_norm = matplotlib.colors.Normalize(vmin = shmin, vmax = shmax)
                    lh_norm = matplotlib.colors.Normalize(vmin = lhmin, vmax = lhmax)
                    
                    ws_map = matplotlib.cm.ScalarMappable(norm = ws_norm, cmap = uni_cmap)
                    sh_map = matplotlib.cm.ScalarMappable(norm = sh_norm, cmap = uni_cmap)
                    lh_map = matplotlib.cm.ScalarMappable(norm = lh_norm, cmap = uni_cmap)
                    #%%
                    model_coors = diff_w_cube.coords('grid_latitude')[0]
                    xmin = np.nanmin(model_coors.points)
                    xmax = np.nanmax(model_coors.points)
                    #%% commence plotting
                        
                    geo_axis = 'grid_latitude'
                    common_keywargs = {'coords' : [geo_axis,'altitude'], 'levels' : 40}#,'norm':colors.CenteredNorm()}
                    com_cbarkeys = {'pad' : 0.005}
                    fig = plt.figure(figsize = (20,28))
                    gs = gridspec.GridSpec(nrows = 7, ncols = 1)
                    top_alt = 4000
                    ## wsp
                    
                    ax0 = fig.add_subplot(gs[0,0])
                    ax0.set_ylim(top = top_alt)
                    ax0.set_title('horizontal windspeed anomaly')
                    ax0.set_xticks([])
                    ax0.set_ylabel('m')
                    ax0.set_xlim(left = xmin, right = xmax)
                    pwsp = iplt.contourf(diff_wind_cube[t], **common_keywargs, cmap = uni_cmap, vmin = wspmin, vmax = wspmax)
                    wsp_cbar = plt.colorbar(wsp_map, ax = ax0, **com_cbarkeys)
                    wsp_cbar.ax.set(ylabel = r'$ms^{-1}$')
                    
                    ax1 = fig.add_subplot(gs[1,0])
                    ax1.set_ylim(top = top_alt)
                    ax1.set_title('vertical windspeed anomaly')
                    ax1.set_xticks([])
                    ax1.set_ylabel('m')
                    ax1.set_xlim(left = xmin, right = xmax)
                    pw = iplt.contourf(diff_w_cube[t], **common_keywargs, cmap = uni_cmap, vmin = wmin, vmax = wmax)
                    w_cbar = plt.colorbar(w_map, ax = ax1, **com_cbarkeys)
                    w_cbar.ax.set(ylabel = r'$ms^{-1}$')
                    
                    ax2 = fig.add_subplot(gs[2,0])
                    ax2.set_ylim(top = top_alt)
                    ax2.set_title('potential temperature anomaly')
                    ax2.set_xticks([])
                    ax2.set_ylabel('m')
                    ax2.set_xlim(left = xmin, right = xmax)
                    pw = iplt.contourf(diff_theta_cube[t], **common_keywargs, cmap = uni_cmap, vmin = Tmin, vmax = Tmax)
                    theta_cbar = plt.colorbar(theta_map, ax = ax2, **com_cbarkeys)
                    theta_cbar.ax.set(ylabel = 'K')
                    
                    ax3 = fig.add_subplot(gs[3,0])
                    ax3.set_ylim(top = top_alt)
                    ax3.set_title('specific humidity anomaly')
                    ax3.set_xticks([])
                    ax3.set_ylabel('m')
                    ax3.set_xlim(left = xmin, right = xmax)
                    pw = iplt.contourf(diff_q_cube[t]*1000, **common_keywargs, cmap = uni_cmap, vmin = qmin, vmax = qmax)
                    q_cbar = plt.colorbar(q_map, ax = ax3, **com_cbarkeys)
                    q_cbar.ax.set(ylabel = '$gkg^{-1}$')
                    
                    ax4 = fig.add_subplot(gs[4,0])
                    ax4.set_ylim(top = top_alt)
                    ax4.set_title('windstress anomaly')
                    ax4.set_xticks([])
                    ax4.set_ylabel('m')
                    ax4.set_xlim(left = xmin, right = xmax)
                    pws = iplt.contourf(diff_ws_cube[t], **common_keywargs, cmap = uni_cmap, vmin = wsmin, vmax = wsmax)
                    ws_cbar = plt.colorbar(ws_map, ax = ax4, **com_cbarkeys)
                    
                    ax5 = fig.add_subplot(gs[5,0])
                    ax5.set_ylim(top = top_alt)
                    ax5.set_title('sensible heat flux anomaly')
                    ax5.set_xticks([])
                    ax5.set_ylabel('m')
                    ax5.set_xlim(left = xmin, right = xmax)
                    psh = iplt.contourf(diff_sh_cube[t], **common_keywargs, cmap = uni_cmap, vmin = shmin, vmax = shmax)
                    sh_cbar = plt.colorbar(sh_map, ax = ax5, **com_cbarkeys)
                    
                    ax6 = fig.add_subplot(gs[6,0])
                    ax6.set_ylim(top = top_alt)
                    ax6.set_title('latent heat flux anomaly')
                    ax6.set_xticks([])
                    ax6.set_ylabel('m')
                    ax6.set_xlim(left = xmin, right = xmax)
                    plh = iplt.contourf(diff_lh_cube[t], **common_keywargs, cmap = uni_cmap, vmin = lhmin, vmax = lhmax)
                    lh_cbar = plt.colorbar(sh_map, ax = ax6, **com_cbarkeys)
                    
                    ax6.set_xticks(tick_coors)
                    ax6.set_xticklabels(tick_dists)
                    ax6.set_xlabel('Scalar distance, km')
                    
                    fig.suptitle(f'CONTROL vs {config} - leg {leg} - {dates[t]}')
                    
                    plt.tight_layout()
                    
                    if True:
                        timepoint = dates[t].strftime('H%HM%M')
                        plt.savefig(f'D:/Project/Figures/PNG/306/{config}/{suite}/vertical/subtractions/{res}/mass/leg{leg}/{config}_CONTROL_subtraction_vert_xsections_{timepoint}_306_leg{leg}.png')
                        
                    plt.close()
        
        
        #%% start aggregate plots section
        aggregate_6hour_plots = True
        if aggregate_6hour_plots:
            ## aggregate cube data in 6 hour means
            for chunk in range(8):
                ss = 6
                mean_wind_cube = diff_wind_cube[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
                mean_theta_cube = diff_theta_cube[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
                mean_q_cube = diff_q_cube[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
                mean_ws_cube = diff_ws_cube[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
                mean_sh_cube = diff_sh_cube[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
                mean_lh_cube = diff_lh_cube[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
                mean_w_cube = diff_w_cube[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
                print(mean_theta_cube)
                
                surface = mean_wind_cube.coords('surface_altitude')[0].points
                x_coors = mean_wind_cube.coords('grid_latitude')[0].points
                
                #%% set colorbar normalisations
                uni_cmap = 'bwr'             
                
                ## atmostate
                wspmin = -10; wspmax = -wspmin
                wmin = -3; wmax = -wmin
                Tmin = -3; Tmax = -Tmin
                qmin = -2; qmax = -qmin
                q_levels = np.arange(0,50)*0.5
                wsp_norm = matplotlib.colors.Normalize(vmin = wspmin, vmax = wspmax)
                w_norm = matplotlib.colors.Normalize(vmin = wmin, vmax = wmax)
                theta_norm = matplotlib.colors.Normalize(vmin = Tmin, vmax = Tmax)
                q_norm = matplotlib.colors.Normalize(vmin = qmin, vmax = qmax)
                
                wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = uni_cmap)
                w_map = matplotlib.cm.ScalarMappable(norm = w_norm, cmap = uni_cmap)
                theta_map = matplotlib.cm.ScalarMappable(norm = theta_norm, cmap = uni_cmap)
                q_map = matplotlib.cm.ScalarMappable(norm = q_norm, cmap = uni_cmap)
                
                ## fluxes
                wsmin = -3; wsmax = -wsmin
                shmin = -400; shmax = -shmin
                lhmin = -400; lhmax = -lhmin
                
                ws_norm = matplotlib.colors.Normalize(vmin = wsmin, vmax = wsmax)
                sh_norm = matplotlib.colors.Normalize(vmin = shmin, vmax = shmax)
                lh_norm = matplotlib.colors.Normalize(vmin = lhmin, vmax = lhmax)
                
                ws_map = matplotlib.cm.ScalarMappable(norm = ws_norm, cmap = uni_cmap)
                sh_map = matplotlib.cm.ScalarMappable(norm = sh_norm, cmap = uni_cmap)
                lh_map = matplotlib.cm.ScalarMappable(norm = lh_norm, cmap = uni_cmap)
                
                #%%
                model_coors = diff_w_cube.coords('grid_latitude')[0]
                xmin = np.nanmin(model_coors.points)
                xmax = np.nanmax(model_coors.points)
                    
                #%% commence plotting
                        
                geo_axis = 'grid_latitude'
                common_keywargs = {'coords' : [geo_axis,'altitude'], 'levels' : 40}#,'norm':colors.CenteredNorm()}
                com_cbarkeys = {'pad' : 0.005}
                fig = plt.figure(figsize = (20,28))
                gs = gridspec.GridSpec(nrows = 7, ncols = 1)
                top_alt = 4000
                ## wsp
                
                ax0 = fig.add_subplot(gs[0,0])
                ax0.set_ylim(top = top_alt)
                ax0.set_title('horizontal windspeed anomaly')
                ax0.set_xticks([])
                ax0.set_ylabel('m')
                ax0.set_xlim(left = xmin, right = xmax)
                pwsp = iplt.contourf(mean_wind_cube, **common_keywargs, cmap = uni_cmap, vmin = wspmin, vmax = wspmax)
                wsp_cbar = plt.colorbar(wsp_map, ax = ax0, **com_cbarkeys)
                wsp_cbar.ax.set(ylabel = r'$ms^{-1}$') 
                
                ax1 = fig.add_subplot(gs[1,0])
                ax1.set_ylim(top = top_alt)
                ax1.set_title('vertical windspeed anomaly')
                ax1.set_xticks([])
                ax1.set_ylabel('m')
                ax1.set_xlim(left = xmin, right = xmax)
                pw = iplt.contourf(mean_w_cube, **common_keywargs, cmap = uni_cmap, vmin = wmin, vmax = wmax)
                w_cbar = plt.colorbar(w_map, ax = ax1, **com_cbarkeys)
                w_cbar.ax.set(ylabel = r'$ms^{-1}$')
                
                ax2 = fig.add_subplot(gs[2,0])
                ax2.set_ylim(top = top_alt)
                ax2.set_title('potential temperature anomaly')
                ax2.set_xticks([])
                ax2.set_ylabel('m')
                ax2.set_xlim(left = xmin, right = xmax)
                pw = iplt.contourf(mean_theta_cube, **common_keywargs, cmap = uni_cmap, vmin = Tmin, vmax = Tmax)
                theta_cbar = plt.colorbar(theta_map, ax = ax2, **com_cbarkeys)
                theta_cbar.ax.set(ylabel = 'K')
                
                ax3 = fig.add_subplot(gs[3,0])
                ax3.set_ylim(top = top_alt)
                ax3.set_title('specific humidity anomaly')
                ax3.set_xticks([])
                ax3.set_ylabel('m')
                ax3.set_xlim(left = xmin, right = xmax)
                pw = iplt.contourf(mean_q_cube*1000, **common_keywargs, cmap = uni_cmap, vmin = qmin, vmax = qmax)
                q_cbar = plt.colorbar(q_map, ax = ax3, **com_cbarkeys)
                q_cbar.ax.set(ylabel = '$gkg^{-1}$')
                
                ax4 = fig.add_subplot(gs[4,0])
                ax4.set_ylim(top = top_alt)
                ax4.set_title('windstress anomaly')
                ax4.set_xticks([])
                ax4.set_ylabel('m')
                ax4.set_xlim(left = xmin, right = xmax)
                pws = iplt.contourf(mean_ws_cube, **common_keywargs, cmap = uni_cmap, vmin = wsmin, vmax = wsmax)
                ws_cbar = plt.colorbar(ws_map, ax = ax4, **com_cbarkeys)
                
                ax5 = fig.add_subplot(gs[5,0])
                ax5.set_ylim(top = top_alt)
                ax5.set_title('sensible heat flux anomaly')
                ax5.set_xticks([])
                ax5.set_ylabel('m')
                ax5.set_xlim(left = xmin, right = xmax)
                psh = iplt.contourf(mean_sh_cube, **common_keywargs, cmap = uni_cmap, vmin = shmin, vmax = shmax)
                sh_cbar = plt.colorbar(sh_map, ax = ax5, **com_cbarkeys)
                
                ax6 = fig.add_subplot(gs[6,0])
                ax6.set_ylim(top = top_alt)
                ax6.set_title('latent heat flux anomaly')
                ax6.set_xticks([])
                ax6.set_ylabel('m')
                ax6.set_xlim(left = xmin, right = xmax)
                plh = iplt.contourf(mean_lh_cube, **common_keywargs, cmap = uni_cmap, vmin = lhmin, vmax = lhmax)
                lh_cbar = plt.colorbar(sh_map, ax = ax6, **com_cbarkeys)
                
                ax6.set_xticks(tick_coors)
                ax6.set_xticklabels(tick_dists)
                ax6.set_xlabel('Scalar distance, km')
                
                timepoint = mean_wind_cube.coord('time').bounds[0]
                print(dates[chunk*ss:ss+ss*chunk])
                boundstart = dates[chunk*ss:ss+ss*chunk][0].strftime('H%HM%M')
                print(boundstart)
                boundend = dates[chunk*ss:ss+ss*chunk][-1].strftime('H%HM%M')
                
                ## add in surface altitude outline
                for k in [ax0,ax1,ax2,ax3,ax4,ax5,ax6]:
                    k.plot(x_coors,surface, color = 'k', linewidth = 2)
                
                fig.suptitle(f'CONTROL vs {config} - leg {leg} - 3 hour mean {boundstart}-{boundend}')
                
                plt.tight_layout()
                
                if True:
                    
                    plt.savefig(f'D:/Project/Figures/PNG/306/{config}/{suite}/vertical/subtractions/{res}/three_hour/leg{leg}/{config}_CONTROL_subtraction_vert_xsections_{boundstart}_306_leg{leg}.png')
                      
                #plt.close()