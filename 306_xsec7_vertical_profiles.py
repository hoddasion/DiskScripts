# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 14:03:52 2021

@author: kse18nru

A script to specifically generate vertical temprature profiles along the leg 7 transsect.

Need geographical location for referencce
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
#%%

resolutions = ['0p5km','1p5km','4p4km']
suite = 'u-cc134'
#expt = 'Control'
config = 'RA1M'

#%%
for res in resolutions:
    print(f'{res}')
    #%% start by getting the locations right first
    ## need 2D orography or landmask field
    orog_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{config}_{res}_um_surface_altitude_0_24hrs_pi_306.nc', 'surface_altitude')
    print(orog_cube)
    
    orog_data = orog_cube.data
    orog_lon = orog_cube.coords('grid_longitude')[0].points
    orog_lat = orog_cube.coords('grid_latitude')[0].points
    
    
    theta_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/air_potential_temperature/{config}_vertslice_{res}_air_potential_temperature_flt306_leg7long.nc', 'air_potential_temperature')[:,:40]
    xwind_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/x_wind/{config}_vertslice_{res}_x_wind_flt306_leg7long.nc', 'x_wind')[:,:40]
    ywind_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/y_wind/{config}_vertslice_{res}_y_wind_flt306_leg7long.nc', 'y_wind')[:,:40]
    wwind_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/upward_air_velocity/{config}_vertslice_{res}_upward_air_velocity_flt306_leg7long.nc', 'upward_air_velocity')[:,:40]
    hoz_wind_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/upward_air_velocity/{config}_vertslice_{res}_upward_air_velocity_flt306_leg7long.nc', 'upward_air_velocity')[:,:40]
    hoz_wind_data = (xwind_cube.data**2+ ywind_cube.data**2)**0.5
    hoz_wind_cube.data = hoz_wind_data
    q_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/specific_humidity/{config}_vertslice_{res}_specific_humidity_flt306_leg7long.nc', 'specific_humidity')[:,:40]
    #print(hoz_wind_data)
    mag_wind = (wwind_cube**2 + hoz_wind_cube**2)**0.5
    
    #%% isolate time coords
    time_coor = theta_cube.coords('time')[0]
    dates = time_coor.units.num2date(time_coor.points)
    
    # hours since midnight; the added 0.5 is to correct for the first output being at 00:30
    hrs_s_midn = time_coor.points - time_coor.points[0] + 0.5
    print(hrs_s_midn)
    #%% get scalar distance ticklist for x axis labelling
    tick_coors, tick_dists = workshop.dist_tick_generator(theta_cube, 'grid_latitude', ticknumber = 6, decprecision = 2)
    #%%
    Tmin = np.nanmin(theta_cube.data); Tmax = np.nanmax(theta_cube.data)
    wspmin = 0; wspmax = np.nanpercentile(hoz_wind_data, 99.99)
    theta_norm = matplotlib.colors.Normalize(vmin = Tmin, vmax = Tmax)
    theta_map = matplotlib.cm.ScalarMappable(norm = theta_norm, cmap = 'plasma')
    wsp_norm = matplotlib.colors.Normalize(vmin = wspmin, vmax = wspmax)
    wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = 'Oranges')
    #%%
    plot_theta_long_transect = False
    if plot_theta_long_transect:
        for tidx in range(48):
            #%% exploratory plotting
            levels = np.arange(0,300,2)
            common_keywargs = {'coords' : ['grid_latitude','altitude'], 'levels' : 20}#levels}
            
            fig, ax = plt.subplots(1,1, figsize = (15,4))
            line_c = iplt.contourf(theta_cube[tidx], **common_keywargs, cmap = 'plasma',vmin = Tmin, vmax = Tmax)#, colors = 'k')
            ax.plot(theta_cube.coords('grid_latitude')[0].points,theta_cube.coords('surface_altitude')[0].points, color = 'k', linestyle = 'dashed')
            ax.set_ylim(top = 5300)
            
            theta_cbar = plt.colorbar(theta_map, ax = ax, pad = 0.005)
            theta_cbar.ax.set(ylabel = 'K')
          
            ax.set_xticks(tick_coors)
            ax.set_xticklabels(tick_dists)
            ax.set_yticks([1000,2000,3000,4000,5000])
            ax.set_yticklabels([1,2,3,4,5])
            ax.set_ylabel('Altitude [km]')
            ax.set_xlabel('Horz. Distance [km]')
            ax.set_title(f'{res} Theta on leg 7 long transect - {dates[tidx]}')
        
            savefig = True
            if savefig:
                print(f'Saving {dates[tidx]}')
                hours = int(hrs_s_midn[tidx] //1)
                minutes = int((hrs_s_midn[tidx] - hours) *60 )
                plt.savefig(f'D:/Project/Figures/PNG/306/{suite}/Vertical/theta/leg7_long_transect/{res}/color/{config}_{res}_air_potential_temperature_{suite[2:]}_TH{hours}M{minutes}_flt306.png')
                plt.close()
    #%%
    plot_theta_long_transect_bw = False
    if plot_theta_long_transect_bw:
        for tidx in range(48):
            #%% exploratory plotting
            levels = np.arange(0,300,2)
            common_keywargs = {'coords' : ['grid_latitude','altitude'], 'levels' : levels}
            
            fig, ax = plt.subplots(1,1, figsize = (15,4))
            line_c = iplt.contour(theta_cube[tidx], **common_keywargs, colors = 'k')
            ax.plot(theta_cube.coords('grid_latitude')[0].points,theta_cube.coords('surface_altitude')[0].points, color = 'k', linestyle = 'dashed')
            ax.set_ylim(top = 5300)
            
            
             # Use the line contours to place contour labels.
       
            ax.clabel(
                line_c,  # Typically best results when labelling line contours.
                colors=['black'],
                manual=False,  # Automatic placement vs manual placement.
                inline=False,  # Cut the line where the label will be placed.
                fmt=' {:.0f}K '.format,  # Labes as integers, with some extra space.
                )
            ax.set_xticks(tick_coors)
            ax.set_xticklabels(tick_dists)
            ax.set_yticks([1000,2000,3000,4000,5000])
            ax.set_yticklabels([1,2,3,4,5])
            ax.set_ylabel('Altitude [km]')
            ax.set_xlabel('Horz. Distance [km]')
            ax.set_title(f'{res} Theta on leg 7 long transect - {dates[tidx]}')
        
            savefig = True
            if savefig:
                print(f'Saving {dates[tidx]}')
                hours = int(hrs_s_midn[tidx] //1)
                minutes = int((hrs_s_midn[tidx] - hours) *60 )
                plt.savefig(f'D:/Project/Figures/PNG/306/{suite}/Vertical/theta/leg7_long_transect/{res}/bw/{config}_{res}_air_potential_temperature_{suite[2:]}_TH{hours}M{minutes}_bw__flt306.png')
                plt.close()

    
    
    #%%
    plot_map = False
    if plot_map:
        res = '0p5km'
        if res == '0p5km':
            fig, ax = plt.subplots(1,1, figsize = (8,14))
            ax.contour(orog_lon, orog_lat, orog_data, color = 'k')
            q0 = ax.plot(theta_cube.coords('grid_longitude')[0].points, theta_cube.coords('grid_latitude')[0].points, color = 'black', linewidth = 5)
            ax.set_title('Leg 7 long transect - 0p5km orography')
            ax.set_xlim(right = orog_lon[-1])
            ax.set_xticks([])
            ax.set_yticks([])
            plt.savefig(f'D:/Project/Figures/PNG/306/{suite}/Vertical/theta/leg7_long_transect/leg7_long_transect_location_simple.png')
            
    plot_theta_profiles = False
    if plot_theta_profiles:
        #%%
        theta_lon = theta_cube.coords('grid_longitude')[0].points
        theta_lat = theta_cube.coords('grid_latitude')[0].points 
        
        for tidx in range(48):
            #%% subset profile coords
            ## define three points: upstream (UP), middle fjord (FJ), and downstream (DN)
            ## need to use nearest neighbour interpolation, or linear to minimise deviation error at coarser resolutions
            ## units in grid latitudes
            UPlat = -0.7
            DNlat = 1.1
            FJlat = 0.170
            theta_coloumns = pert.interpolate_timeseries_coloumns1D(theta_cube,[UPlat,FJlat,DNlat], 'grid_latitude', method = 'linear')
            tcoloumn_data = theta_coloumns.data
            theta_alt = theta_coloumns.coords('altitude')[0].points
            print(theta_coloumns[30,:,1])
            print(np.shape(theta_alt))
            #%% plot a cross-section first and impose profile locations onto it for reference
            
            levels = np.arange(0,300,2)
            common_keywargs = {'coords' : ['grid_latitude','altitude'], 'levels' : 20}#levels}
            
            fig = plt.figure(figsize=(8, 16)) 
            gs = gridspec.GridSpec(2,1, height_ratios=[4, 1])
            
            ## plot profiles
            ax0 = plt.subplot(gs[0])
            ax0.plot(tcoloumn_data[tidx,:,0], theta_alt[:,0], linewidth = 3, label = 'Upwind')
            ax0.plot(tcoloumn_data[tidx,:,1], theta_alt[:,1], linewidth = 3, label = 'Fjord')
            ax0.plot(tcoloumn_data[tidx,:,2], theta_alt[:,2], linewidth = 3, label = 'Downwind')
            ax0.legend()
            ax0.set_yticks([1000,1500,2000,2500,3000,3500,4000,4500,5000])
            ax0.set_yticklabels([1,'',2,'',3,'',4,'',5])
            ax0.set_ylabel('Altitude [km]')
            ax0.set_xlabel(r'$\theta$ [K]')
            ax0.set_ylim(top = 5300, bottom = 0)
            ax0.set_xlim(left = 272, right = 302.5)
            ax0.grid(True)
            
            ## plot reference transect
            ax1 = plt.subplot(gs[1])
            line_c = iplt.contourf(theta_cube[tidx], **common_keywargs, cmap = 'plasma',vmin = Tmin, vmax = Tmax)#, colors = 'k')
            ax1.plot(theta_cube.coords('grid_latitude')[0].points,theta_cube.coords('surface_altitude')[0].points, color = 'k', linestyle = 'dashed')
            ax1.plot(np.ones(len(theta_alt))*UPlat,theta_alt, color = 'k')
            ax1.plot(np.ones(len(theta_alt))*FJlat,theta_alt, color = 'k')
            ax1.plot(np.ones(len(theta_alt))*DNlat,theta_alt, color = 'k')
            
            ax1.set_ylim(top = 5300, bottom = 0)
            ax1.set_xticks([])#tick_coors)
            #ax1.set_xticklabels(tick_dists)
            ax1.set_yticks([1000,2000,3000,4000,5000])
            ax1.set_yticklabels([1,2,3,4,5])
            ax1.set_ylabel('Altitude [km]')
            #ax1.set_xlabel('Horz. Distance [km]')
            theta_cbar = plt.colorbar(theta_map, ax = ax1, orientation = 'horizontal', pad = 0.05)
            theta_cbar.ax.set(xlabel =r'$\theta$ [K]')
            
            ax1.text(-0.68,4500,'UP',fontsize = 25)
            ax1.text(0.19,4500,'FJ',fontsize = 25)
            ax1.text(0.93,4500,'DN',fontsize = 25)
            
            fig.suptitle(f'{res} Theta profiles leg 7 long transect - {dates[tidx]}', fontsize = 17)
            fig.tight_layout()
            
            savefig = True
            if savefig:
                print(f'Saving {dates[tidx]}')
                hours = int(hrs_s_midn[tidx] //1)
                minutes = int((hrs_s_midn[tidx] - hours) *60 )
                
                
                plt.savefig(f'D:/Project/Figures/PNG/306/{suite}/Vertical/theta/leg7_long_transect/{res}/profiles/{config}_{res}_air_potential_temperature_{suite[2:]}_TH{hours}M{minutes}_profiles__flt306.png')
                
                
                plt.close()
            else:
                plt.show()
                
    plot_wsp_profiles = True
    if plot_wsp_profiles:
        #%%
        theta_lon = theta_cube.coords('grid_longitude')[0].points
        theta_lat = theta_cube.coords('grid_latitude')[0].points 
        
        for tidx in range(48):
            #%% subset profile coords
            ## define three points: upstream (UP), middle fjord (FJ), and downstream (DN)
            ## need to use nearest neighbour interpolation, or linear to minimise deviation error at coarser resolutions
            ## units in grid latitudes
            UPlat = -0.7
            DNlat = 1.1
            FJlat = 0.170
            wsp_coloumns = pert.interpolate_timeseries_coloumns1D(hoz_wind_cube,[UPlat,FJlat,DNlat], 'grid_latitude', method = 'linear')
            wspcoloumn_data = wsp_coloumns.data
            wsp_alt = wsp_coloumns.coords('altitude')[0].points
            
            #%% plot a cross-section first and impose profile locations onto it for reference
            
            levels = np.arange(0,300,2)
            common_keywargs = {'coords' : ['grid_latitude','altitude'], 'levels' : 20}#levels}
            
            fig = plt.figure(figsize=(8, 16)) 
            gs = gridspec.GridSpec(2,1, height_ratios=[4, 1])
            
            ## plot profiles
            ax0 = plt.subplot(gs[0])
            ax0.plot(wspcoloumn_data[tidx,:,0], wsp_alt[:,0], linewidth = 3, label = 'Upwind')
            ax0.plot(wspcoloumn_data[tidx,:,1], wsp_alt[:,1], linewidth = 3, label = 'Fjord')
            ax0.plot(wspcoloumn_data[tidx,:,2], wsp_alt[:,2], linewidth = 3, label = 'Downwind')
            ax0.legend()
            ax0.set_yticks([1000,1500,2000,2500,3000,3500,4000,4500,5000])
            ax0.set_yticklabels([1,'',2,'',3,'',4,'',5])
            ax0.set_ylabel('Altitude [km]')
            ax0.set_xlabel(r'$wsp$ [$ms^{-1}$]')
            ax0.set_ylim(top = 5300, bottom = 0)
            ax0.set_xlim(left = 0, right = 30)
            ax0.grid(True)
            
            ## plot reference transect
            ax1 = plt.subplot(gs[1])
            line_c = iplt.contourf(hoz_wind_cube[tidx], **common_keywargs, cmap = 'Oranges',vmin = wspmin, vmax = wspmax)#, colors = 'k')
            ax1.plot(theta_cube.coords('grid_latitude')[0].points,theta_cube.coords('surface_altitude')[0].points, color = 'k', linestyle = 'dashed')
            ax1.plot(np.ones(len(wsp_alt))*UPlat,wsp_alt, color = 'k')
            ax1.plot(np.ones(len(wsp_alt))*FJlat,wsp_alt, color = 'k')
            ax1.plot(np.ones(len(wsp_alt))*DNlat,wsp_alt, color = 'k')
            
            ax1.set_ylim(top = 5300, bottom = 0)
            ax1.set_xticks([])#tick_coors)
            #ax1.set_xticklabels(tick_dists)
            ax1.set_yticks([1000,2000,3000,4000,5000])
            ax1.set_yticklabels([1,2,3,4,5])
            ax1.set_ylabel('Altitude [km]')
            #ax1.set_xlabel('Horz. Distance [km]')
            theta_cbar = plt.colorbar(wsp_map, ax = ax1, orientation = 'horizontal', pad = 0.05)
            theta_cbar.ax.set(xlabel =r'$wsp$ [$ms^{-1}$]')
            
            ax1.text(-0.68,4500,'UP',fontsize = 25)
            ax1.text(0.19,4500,'FJ',fontsize = 25)
            ax1.text(0.93,4500,'DN',fontsize = 25)
            
            fig.suptitle(f'{res} wsp profiles leg 7 long transect - {dates[tidx]}', fontsize = 17)
            fig.tight_layout()
            
            savefig = False
            if savefig:
                print(f'Saving {dates[tidx]}')
                hours = int(hrs_s_midn[tidx] //1)
                minutes = int((hrs_s_midn[tidx] - hours) *60 )
                
                
                plt.savefig(f'D:/Project/Figures/PNG/306/{suite}/Vertical/wsp/leg7_long_transect/{res}/profiles/{config}_{res}_wsp_{suite[2:]}_TH{hours}M{minutes}_profiles__flt306.png')
                
                
                plt.close()
            else:
                plt.show()