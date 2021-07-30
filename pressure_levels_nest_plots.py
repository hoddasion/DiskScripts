# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 09:56:24 2021

@author: kse18nru
"""

import iris
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import iris.plot as iplt
import xsection_workshop as workshop
import math
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

#%% Global settings
plt.close()
plt.rcParams.update({'font.size': 20})
flight = 306
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


#%% main
suite = 'u-cc134'
datapath = f'D:/Project/Model_Data/{suite}/'
filename0p5 = f'RA1M_0p5km_cc134_24hrs_pressurevars_{flight}.nc'

filename1p5 = f'RA1M_1p5km_cc134_24hrs_pressurevars_{flight}.nc'

filename4p4 = f'RA1M_4p4km_cc134_24hrs_pressurevars_{flight}.nc'

#%%
wsp_plots = False
if wsp_plots:
    #t = 5; h = 7
    timestep_idx = 2
    
    u0p5 = iris.load_cube(f'{datapath}{filename0p5}', 'x_wind')[2::timestep_idx]
    v0p5 = iris.load_cube(f'{datapath}{filename0p5}', 'y_wind')[2::timestep_idx]
    wsp0p5 = (u0p5**2+v0p5**2)**0.5
    u1p5 = iris.load_cube(f'{datapath}{filename1p5}', 'x_wind')[2::timestep_idx]
    v1p5 = iris.load_cube(f'{datapath}{filename1p5}', 'y_wind')[2::timestep_idx]
    wsp1p5 = (u1p5**2+v1p5**2)**0.5
    u4p4 = iris.load_cube(f'{datapath}{filename4p4}', 'x_wind')[2::timestep_idx]
    v4p4 = iris.load_cube(f'{datapath}{filename4p4}', 'y_wind')[2::timestep_idx]
    wsp4p4 = (u4p4**2+v4p4**2)**0.5
    landcube0p5 = iris.load_cube(f'{datapath}{filename0p5}', 'land_binary_mask')
    landcube1p5 = iris.load_cube(f'{datapath}{filename1p5}', 'land_binary_mask')
    landcube4p4 = iris.load_cube(f'{datapath}{filename4p4}', 'land_binary_mask')
    
    geop0p5 = iris.load_cube(f'{datapath}{filename0p5}', 'geopotential_height')[2::timestep_idx]
    geop1p5 = iris.load_cube(f'{datapath}{filename1p5}', 'geopotential_height')[2::timestep_idx]
    geop4p4 = iris.load_cube(f'{datapath}{filename4p4}', 'geopotential_height')[2::timestep_idx]
    #%% extract data and coords
    lat0p5 = wsp0p5.coords('grid_latitude')[0].points
    lon0p5 = wsp0p5.coords('grid_longitude')[0].points
    lat1p5 = wsp1p5.coords('grid_latitude')[0].points
    lon1p5 = wsp1p5.coords('grid_longitude')[0].points
    lat4p4 = wsp4p4.coords('grid_latitude')[0].points
    lon4p4 = wsp4p4.coords('grid_longitude')[0].points
    
    data0p5 = wsp0p5.data
    data1p5 = wsp1p5.data
    data4p4 = wsp4p4.data    
    data0p5[np.where(data0p5 == 0)] = np.nan
    data1p5[np.where(data1p5 == 0)] = np.nan
    data4p4[np.where(data4p4 == 0)] = np.nan
    udata0p5 = u0p5.data
    udata1p5 = u1p5.data
    udata4p4 = u4p4.data
    vdata0p5 = v0p5.data
    vdata1p5 = v1p5.data
    vdata4p4 = v4p4.data
    udata0p5[np.where(udata0p5 == 0)] = np.nan
    udata1p5[np.where(udata1p5 == 0)] = np.nan
    udata4p4[np.where(udata4p4 == 0)] = np.nan
    vdata0p5[np.where(vdata0p5 == 0)] = np.nan
    vdata1p5[np.where(vdata1p5 == 0)] = np.nan
    vdata4p4[np.where(vdata4p4 == 0)] = np.nan
   
    landlon4p4 = landcube4p4.coords('grid_longitude')[0].points
    landlat4p4 = landcube4p4.coords('grid_latitude')[0].points
    
    
    land4p4 = landcube4p4.data
    
    #%%
    for h in range(16):
        pressure = u4p4.coords('pressure')[0].points[h]
        for t in range(11):
            cube_times = u4p4.coords('time')[0]
            dates = cube_times.units.num2date(cube_times.points)
            timepoints = []
            for m in range(len(dates)):
                timepoint = dates[m].strftime('%H:%M')
                timepoints.append(timepoint)
            times = np.array(timepoints)
            print(times)
            #%% normalisations
            norm4p4, map4p4 = dynamic_normalisation(wsp4p4, h, 'Greens', start_at_zero = True)
            norm1p5, map1p5 = dynamic_normalisation(wsp1p5, h, 'Oranges', start_at_zero = True)
            norm0p5, map0p5 = dynamic_normalisation(wsp0p5, h, 'Blues', start_at_zero = True)
            
            #%%
            ystep4 = 15; xstep4 = 15
            ystep1 = 30; xstep1 = 30
            ystep0 = 40; xstep0 = 40
            geolevels = np.arange(0,20000,10)
            #%% plot
            kwargs = {'shading' : 'auto'}
            fig,(ax, cax0, cax1,cax2) = plt.subplots(1,4, figsize = (24,18), gridspec_kw={'width_ratios':[1,0.05,0.05,0.05]})
            q0 = ax.pcolormesh(lon4p4, lat4p4,data4p4[t][h], norm = norm4p4, cmap = 'Greens', **kwargs)
            ax.contour(landlon4p4, landlat4p4, land4p4, colors = 'k')
            q8 = ax.contour(lon4p4,lat4p4, geop4p4.data[t,h], levels = geolevels, colors = 'k')
            ax.quiver(lon4p4[::xstep4], lat4p4[::ystep4], udata4p4[t,h,::xstep4,::ystep4], vdata4p4[t,h,::xstep4,::ystep4], 
                      norm = norm4p4, pivot = 'mid', width = 0.002)
            plt.clabel(q8,fmt = '%1.0f')
            
            q1 = ax.pcolormesh(lon1p5, lat1p5,data1p5[t][h], norm = norm1p5, cmap = 'Oranges', **kwargs)
            ax.quiver(lon1p5[::xstep1], lat1p5[::ystep1], udata1p5[t,h,::xstep1,::ystep1], vdata1p5[t,h,::xstep1,::ystep1], 
                      norm = norm4p4, width = 0.002, pivot = 'mid')
            
            
            
            q2 = ax.pcolormesh(lon0p5, lat0p5,data0p5[t][h], norm = norm0p5, cmap = 'Blues', **kwargs)
            
            ax.set_yticks([])
            ax.set_xticks([])
            
            cbar_keywargs = {'orientation':'vertical', 'pad' : 0.005}
            cbar0 = plt.colorbar(q0, cax =cax0, **cbar_keywargs)
            cbar1 = plt.colorbar(q1, cax =cax1, **cbar_keywargs)
            cbar2 = plt.colorbar(q2, cax =cax2, **cbar_keywargs)
                    
            cbar2.ax.set(ylabel = r'Horz. Windsp. $ms^{-1}$')
            title = f'Horizontal wind {pressure:.0f}hPa {times[t]} UTC'
            fig.suptitle(title, y = 0.91)
            plt.subplots_adjust(wspace=0.13, hspace=0.05)
            plt.savefig(f'D:/Project/Figures/PNG/{flight}/{suite}/P_levels/Control/windvector_pl/nests/{pressure:.0f}hPa/Nests_horz_windspeed_Control_pl{pressure:.0f}hPa_H{times[t][0:2]}M{times[t][3:]}_flt{flight}.png')
            plt.savefig(f'D:/Project/Figures/PDF/{flight}/{suite}/P_levels/Control/windvector_pl/nests/{pressure:.0f}hPa/Nests_horz_windspeed_Control_pl{pressure:.0f}hPa_H{times[t][0:2]}M{times[t][3:]}_flt{flight}.pdf')
            plt.close()
            
w_plots = False
if w_plots:
    
    #t = 5; h = 7
    timestep_idx = 2
    
    
    w0p5 = iris.load_cube(f'{datapath}{filename0p5}', 'upward_air_velocity')[2::timestep_idx]
    
    w1p5 = iris.load_cube(f'{datapath}{filename1p5}', 'upward_air_velocity')[2::timestep_idx]
    
    w4p4 = iris.load_cube(f'{datapath}{filename4p4}', 'upward_air_velocity')[2::timestep_idx]
    
    landcube0p5 = iris.load_cube(f'{datapath}{filename0p5}', 'land_binary_mask')
    landcube1p5 = iris.load_cube(f'{datapath}{filename1p5}', 'land_binary_mask')
    landcube4p4 = iris.load_cube(f'{datapath}{filename4p4}', 'land_binary_mask')
    
    geop0p5 = iris.load_cube(f'{datapath}{filename0p5}', 'geopotential_height')[2::timestep_idx]
    geop1p5 = iris.load_cube(f'{datapath}{filename1p5}', 'geopotential_height')[2::timestep_idx]
    geop4p4 = iris.load_cube(f'{datapath}{filename4p4}', 'geopotential_height')[2::timestep_idx]
    #%% extract data and coords
    lat0p5 = w0p5.coords('grid_latitude')[0].points
    lon0p5 = w0p5.coords('grid_longitude')[0].points
    lat1p5 = w1p5.coords('grid_latitude')[0].points
    lon1p5 = w1p5.coords('grid_longitude')[0].points
    lat4p4 = w4p4.coords('grid_latitude')[0].points
    lon4p4 = w4p4.coords('grid_longitude')[0].points
    
    data0p5 = w0p5.data
    data1p5 = w1p5.data
    data4p4 = w4p4.data    
    # data0p5[np.where(data0p5 == 0)] = np.nan
    # data1p5[np.where(data1p5 == 0)] = np.nan
    # data4p4[np.where(data4p4 == 0)] = np.nan
   
    
   
    landlon4p4 = landcube4p4.coords('grid_longitude')[0].points
    landlat4p4 = landcube4p4.coords('grid_latitude')[0].points
    
    
    land4p4 = landcube4p4.data
    
    #%%
    for j in range(48):
        h = j +10
        pressure = w4p4.coords('pressure')[0].points[h]
        for t in range(11):
            cube_times = w4p4.coords('time')[0]
            dates = cube_times.units.num2date(cube_times.points)
            timepoints = []
            for m in range(len(dates)):
                timepoint = dates[m].strftime('%H:%M')
                timepoints.append(timepoint)
            times = np.array(timepoints)
            print(times)
            #%% normalisations
            norm4p4, map4p4 = dynamic_normalisation(w4p4, h, 'PuOr_r', symm = True, percentile = 0.999)
            norm1p5, map1p5 = dynamic_normalisation(w1p5, h, 'PiYG_r', symm = True, percentile = 0.999)
            norm0p5, map0p5 = dynamic_normalisation(w0p5, h, 'seismic', symm = True, percentile = 0.999)
            
            #%%
            ystep4 = 15; xstep4 = 15
            ystep1 = 30; xstep1 = 30
            ystep0 = 40; xstep0 = 40
            geolevels = np.arange(0,20000,10)
            #%% plot
            kwargs = {'shading' : 'auto'}
            fig,(ax, cax0, cax1,cax2) = plt.subplots(1,4, figsize = (24,18), gridspec_kw={'width_ratios':[1,0.05,0.05,0.05]})
            q0 = ax.pcolormesh(lon4p4, lat4p4,data4p4[t][h], norm = norm4p4, cmap = 'PuOr_r', **kwargs)
            ax.contour(landlon4p4, landlat4p4, land4p4, colors = 'k')
            q8 = ax.contour(lon4p4,lat4p4, geop4p4.data[t,h], levels = geolevels, colors = 'k')
            
            plt.clabel(q8,fmt = '%1.0f')
            
            q1 = ax.pcolormesh(lon1p5, lat1p5,data1p5[t][h], norm = norm1p5, cmap = 'PiYG_r', **kwargs)
            
            
            
            
            q2 = ax.pcolormesh(lon0p5, lat0p5,data0p5[t][h], norm = norm0p5, cmap = 'seismic', **kwargs)
            
            ax.set_yticks([])
            ax.set_xticks([])
            
            cbar_keywargs = {'orientation':'vertical', 'pad' : 0.01}
            cbar0 = plt.colorbar(q0, cax =cax0, **cbar_keywargs)
            cbar1 = plt.colorbar(q1, cax =cax1, **cbar_keywargs)
            cbar2 = plt.colorbar(q2, cax =cax2, **cbar_keywargs)
                    
            cbar2.ax.set(ylabel = r'Vert. Windsp. $ms^{-1}$')
            title = f'Vertical wind {pressure:.0f}hPa {times[t]} UTC'
            fig.suptitle(title, y = 0.91)
            plt.subplots_adjust(wspace=0.2, hspace=0.05)
            plt.savefig(f'D:/Project/Figures/PNG/{flight}/{suite}/P_levels/Control/upward_air_velocity/nests/{pressure:.0f}hPa/Nests_vert_windspeed_Control_pl{pressure:.0f}hPa_H{times[t][0:2]}M{times[t][3:]}_flt{flight}.png')
            plt.savefig(f'D:/Project/Figures/PDF/{flight}/{suite}/P_levels/Control/upward_air_velocity/nests/{pressure:.0f}hPa/Nests_vert_windspeed_Control_pl{pressure:.0f}hPa_H{times[t][0:2]}M{times[t][3:]}_flt{flight}.pdf')
            plt.close()
            print(title)
            
T_plots = True
if T_plots:
    
    #t = 5; h = 7
    timestep_idx = 2
    
    
    T0p5 = iris.load_cube(f'{datapath}{filename0p5}', 'air_temperature')[2::timestep_idx]
    
    T1p5 = iris.load_cube(f'{datapath}{filename1p5}', 'air_temperature')[2::timestep_idx]
    
    T4p4 = iris.load_cube(f'{datapath}{filename4p4}', 'air_temperature')[2::timestep_idx]
    
    landcube0p5 = iris.load_cube(f'{datapath}{filename0p5}', 'land_binary_mask')
    landcube1p5 = iris.load_cube(f'{datapath}{filename1p5}', 'land_binary_mask')
    landcube4p4 = iris.load_cube(f'{datapath}{filename4p4}', 'land_binary_mask')
    
    geop0p5 = iris.load_cube(f'{datapath}{filename0p5}', 'geopotential_height')[2::timestep_idx]
    geop1p5 = iris.load_cube(f'{datapath}{filename1p5}', 'geopotential_height')[2::timestep_idx]
    geop4p4 = iris.load_cube(f'{datapath}{filename4p4}', 'geopotential_height')[2::timestep_idx]
    #%% extract data and coords
    lat0p5 = T0p5.coords('grid_latitude')[0].points
    lon0p5 = T0p5.coords('grid_longitude')[0].points
    lat1p5 = T1p5.coords('grid_latitude')[0].points
    lon1p5 = T1p5.coords('grid_longitude')[0].points
    lat4p4 = T4p4.coords('grid_latitude')[0].points
    lon4p4 = T4p4.coords('grid_longitude')[0].points
    
    data0p5 = T0p5.data
    data1p5 = T1p5.data
    data4p4 = T4p4.data    
    data0p5[np.where(data0p5 == 0)] = np.nan
    data1p5[np.where(data1p5 == 0)] = np.nan
    data4p4[np.where(data4p4 == 0)] = np.nan
   
    
   
    landlon4p4 = landcube4p4.coords('grid_longitude')[0].points
    landlat4p4 = landcube4p4.coords('grid_latitude')[0].points
    
    
    land4p4 = landcube4p4.data
    
    #%%
    for j in range(16):
        h = j + 9
        pressure = T4p4.coords('pressure')[0].points[h]
        for t in range(11):
            cube_times = T4p4.coords('time')[0]
            dates = cube_times.units.num2date(cube_times.points)
            timepoints = []
            for m in range(len(dates)):
                timepoint = dates[m].strftime('%H:%M')
                timepoints.append(timepoint)
            times = np.array(timepoints)
            print(times)
            #%% normalisations
            norm4p4, map4p4 = dynamic_normalisation(T4p4, h, 'plasma', symm = False, percentile = 0.999)
            norm1p5, map1p5 = dynamic_normalisation(T1p5, h, 'cividis', symm = False, percentile = 0.999)
            norm0p5, map0p5 = dynamic_normalisation(T0p5, h, 'viridis', symm = False, percentile = 0.999)
            
            #%%
            ystep4 = 15; xstep4 = 15
            ystep1 = 30; xstep1 = 30
            ystep0 = 40; xstep0 = 40
            geolevels = np.arange(0,20000,10)
            #%% plot
            kwargs = {'shading' : 'auto'}
            fig,(ax, cax0, cax1,cax2) = plt.subplots(1,4, figsize = (24,18), gridspec_kw={'width_ratios':[1,0.05,0.05,0.05]})
            q0 = ax.pcolormesh(lon4p4, lat4p4,data4p4[t][h], norm = norm4p4, cmap = 'plasma', **kwargs)
            ax.contour(landlon4p4, landlat4p4, land4p4, colors = 'k')
            q8 = ax.contour(lon4p4,lat4p4, geop4p4.data[t,h], levels = geolevels, colors = 'k')
            
            plt.clabel(q8,fmt = '%1.0f')
            
            q1 = ax.pcolormesh(lon1p5, lat1p5,data1p5[t][h], norm = norm1p5, cmap = 'cividis', **kwargs)
            
            
            
            
            q2 = ax.pcolormesh(lon0p5, lat0p5,data0p5[t][h], norm = norm0p5, cmap = 'viridis', **kwargs)
            
            ax.set_yticks([])
            ax.set_xticks([])
            
            cbar_keywargs = {'orientation':'vertical', 'pad' : 0.01}
            cbar0 = plt.colorbar(q0, cax =cax0, **cbar_keywargs)
            cbar1 = plt.colorbar(q1, cax =cax1, **cbar_keywargs)
            cbar2 = plt.colorbar(q2, cax =cax2, **cbar_keywargs)
                    
            cbar2.ax.set(ylabel = r'Air Temp. $K$')
            title = f'Air Temperature {pressure:.0f}hPa {times[t]} UTC'
            fig.suptitle(title, y = 0.91)
            plt.subplots_adjust(wspace=0.2, hspace=0.05)
            plt.savefig(f'D:/Project/Figures/PNG/{flight}/{suite}/P_levels/Control/air_temperature/nests/{pressure:.0f}hPa/Nests_air_temperature_Control_pl{pressure:.0f}hPa_H{times[t][0:2]}M{times[t][3:]}_flt{flight}.png')
            plt.savefig(f'D:/Project/Figures/PDF/{flight}/{suite}/P_levels/Control/air_temperature/nests/{pressure:.0f}hPa/Nests_air_temperature_Control_pl{pressure:.0f}hPa_H{times[t][0:2]}M{times[t][3:]}_flt{flight}.pdf')
            plt.close()
            print(title)