# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 13:19:49 2021

@author: kse18nru
"""
#%%
import iris
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import iris.plot as iplt
import xsection_workshop as workshop
import math
import matplotlib
#%% Global settings
plt.rcParams.update({'font.size': 20})

#%% unique function definiton

def dynamic_normalisation(cube, level_idx, cmap,symm = False, fix_mask = False, start_at_zero = False,  percentile = 1):
    
    # extract data from cube for single level by index
    data = cube.data[:,level_idx]
    
    ## correct masks (e.g. zeros in temperature)
    if fix_mask:
        zeroes = np.where(data == 0) 
        data[zeroes] = np.nan
    ## first symmetric use case
    if symm:
        val_min = -np.nanpercentile(-data,100*percentile)
        val_max = np.nanpercentile(data,100*percentile)
        # identify which has greater magnitude
        if -val_min >= val_max:
            val_max = -val_min
        else:
            val_min = - val_max
        # apply saturation adjustments
        
        
        # proceed with normalisation
        norm = matplotlib.colors.Normalize(vmin = val_min, vmax = val_max)
        new_map = matplotlib.cm.ScalarMappable(norm = norm, cmap = cmap)
            
    else:
        ## second for non-symmetric use case
        if start_at_zero:
            val_min = 0 
        else:
            val_min = -np.nanpercentile(-data, 100*percentile) 
        val_max = np.nanpercentile(data, 100*percentile) # simultaneously adjust for desired saturation at each end 
        norm = matplotlib.colors.Normalize(vmin = val_min, vmax = val_max)
        new_map = matplotlib.cm.ScalarMappable(norm = norm, cmap = cmap)
    
    return norm, new_map
#%% load data

for res in ['0p5km','1p5km','4p4km']:
#   
    suite = 'u-cc134'
    flight = 306
    datapath = f'D:/Project/Model_Data/{suite}/'
    filename = f'RA1M_{res}_cc134_24hrs_pressurevars_{flight}.nc'
    
    time_step_idx = 2
    x_wind = iris.load_cube(f'{datapath}{filename}', 'x_wind')[2::time_step_idx]
    y_wind = iris.load_cube(f'{datapath}{filename}', 'y_wind')[2::time_step_idx]
    air_temperature = iris.load_cube(f'{datapath}{filename}', 'air_temperature')[2::time_step_idx]
    specific_humidity = iris.load_cube(f'{datapath}{filename}', 'specific_humidity')[2::time_step_idx]
    w_wind = iris.load_cube(f'{datapath}{filename}', 'upward_air_velocity')[2::time_step_idx]
    geopotential = iris.load_cube(f'{datapath}{filename}', 'geopotential_height')[2::time_step_idx]
    landmask = iris.load_cube(f'{datapath}{filename}', 'land_binary_mask')
    
    #%% calculate wind magnitude
    mag_wind = (x_wind**2 + y_wind**2)**0.5
    
    cube_times = x_wind.coords('time')[0]
    dates = cube_times.units.num2date(cube_times.points)
    timepoints = []
    for t in range(len(dates)):
        timepoint = dates[t].strftime('%H:%M')
        timepoints.append(timepoint)
    times = np.array(timepoints)
    #print(times)
    
    #%% extract coordinates for var datacube
    var_lon = x_wind.coords('grid_longitude')[0].points
    var_lat = x_wind.coords('grid_latitude')[0].points
    
    #%% create scalar distanc axis labels and posotions
        
    tick_lon_coors, tick_lon_dists = workshop.dist_tick_generator(w_wind, 'grid_longitude', ticknumber = 7, decprecision = 0, flat = True)
    tick_lat_coors, tick_lat_dists = workshop.dist_tick_generator(w_wind, 'grid_latitude', ticknumber = 7, decprecision = 0, flat = True)
    #print(tick_lon_coors, tick_lon_dists)
    temp_lon = []
    temp_lat = []
    for k in range(len(tick_lon_dists)):
        
        temp_lon.append(math.trunc(tick_lon_dists[k]))
    for k in range(len(tick_lat_dists)):
        
        temp_lat.append(math.trunc(tick_lat_dists[k]))
    
    tick_lat_dists = temp_lat
    tick_lon_dists = temp_lon
    #%%
    for j in range(11):
        i = j + 5
        pressure =  x_wind.coords('pressure')[0].points[i]
        
        #%% normalisation at each level
        wspnorm, wspmap = dynamic_normalisation(mag_wind, i, 'Oranges', start_at_zero = True)
        wnorm, wmap = dynamic_normalisation(w_wind, i, 'seismic', symm = True, percentile = 0.9995)
        Tnorm, Tmap = dynamic_normalisation(air_temperature, i, 'plasma', percentile = 0.9995, fix_mask = True)
        qnorm, qmap = dynamic_normalisation(specific_humidity*1000, i, 'Blues', start_at_zero = True)
        
        #%% plot at each tme step
        for t in range(11):
            print(f'{res} {pressure} tidx{t}')
            
            if res == '0p5km':
                figsize = (16,20)
                geolevels = np.arange(0,20000,5)
                xstep = 20; ystep = 30
            if res == '1p5km':
                figsize = (20,20)
                geolevels = np.arange(0,20000,10)
                xstep = 20; ystep = 20
            if res == '4p4km':
                figsize = (20,20)
                geolevels = np.arange(0,20000,20)
                xstep = 20; ystep = 20
                #%%
            fig, axes = plt.subplots(2,2,figsize = figsize)
            #print(axes)
            #fig.suptitle(f'{pressure}hPa')
            #%%
            ax0 = axes[0,0]
            q0 = ax0.pcolormesh(var_lon, var_lat, mag_wind[t][i].data, cmap = 'Oranges', shading = 'auto', norm = wspnorm)
            qq0 = ax0.quiver(var_lon[::xstep],var_lat[::ystep],x_wind.data[t,i,::ystep,::xstep], y_wind.data[t,i, ::ystep, ::xstep])
            q1 = ax0.contour(var_lon, var_lat, geopotential[t][i].data, levels = geolevels, colors = 'k')
            #ax0.set_ylim(bottom = var_lat[0], top = var_lat[-1])
            #ax0.set_xlim(left = var_lon[0], right = var_lon[-1])
            
            ax0.contour(landmask.coords('grid_longitude')[0].points, landmask.coords('grid_latitude')[0].points, landmask.data, colors = 'k')
            
            ax1 = axes[0,1]
            q3 = ax1.pcolormesh(var_lon,var_lat, w_wind.data[t][i], cmap = 'seismic', norm = wnorm, shading  = 'auto')
            q4 = ax1.contour(var_lon, var_lat, geopotential[t][i].data, levels = geolevels, colors = 'k')
            ax1.contour(landmask.coords('grid_longitude')[0].points, landmask.coords('grid_latitude')[0].points, landmask.data, colors = 'k')
    
            
            ax2 = axes[1,0]
            q5 = ax2.pcolormesh(var_lon,var_lat, air_temperature.data[t][i], cmap = 'plasma', norm = Tnorm, shading = 'auto')
            q6 = ax2.contour(var_lon, var_lat, geopotential[t][i].data, levels = geolevels, colors = 'k')
            ax2.contour(landmask.coords('grid_longitude')[0].points, landmask.coords('grid_latitude')[0].points, landmask.data, colors = 'k')
        
            ax3 = axes[1,1]
            q7 = ax3.pcolormesh(var_lon, var_lat, specific_humidity.data[t][i]*1000, cmap = 'Blues',norm = qnorm, shading = 'auto')
            q8 = ax3.contour(var_lon, var_lat, geopotential[t][i].data, levels = geolevels, colors = 'k')
            ax3.contour(landmask.coords('grid_longitude')[0].points, landmask.coords('grid_latitude')[0].points, landmask.data, colors = 'k')
            plt.clabel(q8,fmt = '%1.0f')
            ax0.set_xticks([])
            ax1.set_xticks([]); ax1.set_yticks([])
            ax3.set_yticks([])
            ax0.set_yticks(tick_lat_coors); ax0.set_yticklabels(tick_lat_dists)
            ax2.set_yticks(tick_lat_coors); ax2.set_yticklabels(tick_lat_dists)
            ax2.set_xticks(tick_lon_coors); ax2.set_xticklabels(tick_lon_dists)
            ax3.set_xticks(tick_lon_coors); ax3.set_xticklabels(tick_lon_dists)
            
            cbar_keywargs = {'orientation':'vertical'}
            cbar0 = plt.colorbar(q0, ax = ax0, **cbar_keywargs)
            cbar1 = plt.colorbar(q3, ax = ax1, **cbar_keywargs)
            cbar2 = plt.colorbar(q5, ax = ax2, **cbar_keywargs)
            cbar3 = plt.colorbar(q7, ax = ax3, **cbar_keywargs)
            cbar0.ax.set(ylabel = r'Horz. Windsp. $ms^{-1}$')
            cbar1.ax.set(ylabel = r'Vert. Windsp. $ms^{-1}$')
            cbar2.ax.set(ylabel = r'Temp. $K$')
            cbar3.ax.set(ylabel = r'Spec.Hum. $gkg^{-1}$')
            
            title = f'{res} {pressure:.0f}hPa {times[t]} UTC'
            fig.suptitle(title, y = 0.91)
            plt.subplots_adjust(wspace=0.06, hspace=0.05)
            
            pngpath = f'D:/Project/Figures/PNG/{flight}/{suite}/P_levels/Control/atmostate/fulldomain/{res}/{pressure:.0f}hPa/RA1M_{res}_control_atmostate_p_level_fulldomain_{pressure:.0f}hPa_H{times[t][0:2]}M{times[t][3:]}_{flight}.png'
            print(pngpath)
            plt.savefig(pngpath)
            plt.show()
            
            
           