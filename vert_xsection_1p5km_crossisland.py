# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 12:07:00 2021

@author: kse18nru
"""

#%% module imports
import sys
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
import model_foundry_toolkit as toolkit
import data_management as dmna
import time
import model_foundry as foundry
import matplotlib.gridspec as gridspec
from matplotlib.ticker import FormatStrFormatter
#%% Global settings
plt.rcParams.update({'font.size': 14})


#%% load data
resolutions = ['1p5km']
for res in resolutions:
    max_level = 50
    datapath = 'D:/Project/Model_Data/u-bu807/nc/Control/'
    thetacube = iris.load_cube(f'{datapath}RA1M_{res}_um_air_potential_temperature_24hrs_i_301.nc', 'air_potential_temperature')[:,:max_level]
    ucube = iris.load_cube(f'{datapath}RA1M_{res}_um_m01s15i002_24hrs_h_301.nc', 'm01s15i002')[:,:max_level]
    vcube = iris.load_cube(f'{datapath}RA1M_{res}_um_m01s15i003_24hrs_h_301.nc', 'm01s15i003')[:,:max_level]
    wspcube = (ucube**2+vcube**2)**0.5
    qcube = iris.load_cube(f'{datapath}RA1M_{res}_um_specific_humidity_24hrs_i_301.nc', 'specific_humidity')[:,:max_level]
    wcube = iris.load_cube(f'{datapath}RA1M_{res}_um_upward_air_velocity_24hrs_h_301.nc', 'upward_air_velocity')[:,:max_level]
    orogcube = iris.load_cube(f'{datapath}RA1M_{res}_um_air_potential_temperature_24hrs_i_301.nc', 'surface_altitude')
    land_cube = iris.load_cube('../../Model_Data/u-bu807/nc/Control/surface_altitude/RA1M_1p5km_surface_altitude_fulldomain_flt301.nc', 'land_binary_mask')
    
    
    #%% seperate out coordinates
    gridlon = orogcube.coords('grid_longitude')[0].points
    gridlat = orogcube.coords('grid_latitude')[0].points
    ## single out parallel for xsection
    print('no. of lat =', len(gridlat))
    latidx = 195
    print(f'gridlat at idx {latidx} =', gridlat[latidx])
    ## make new lat array equal length ot lon array
    monolat = np.ones(len(gridlon))*gridlat[latidx]
    lm_data = np.copy(land_cube.data)
    or_data = np.copy(orogcube.data)
    xorog = orogcube[latidx]
    
    #%% get earth coord field
    polelat = wcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = wcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    lons, lats = foundry.modf.unrotate_coords(gridlon, gridlat, polelon, polelat)
    location_plot = False
    if location_plot:
        #%% create arrays for parallel boundary lines
        lower_idx = latidx - 20
        upper_idx = latidx + 20
        lower_monolat = np.ones(len(gridlon))*gridlat[lower_idx]
        upper_monolat = np.ones(len(gridlon))*gridlat[upper_idx]
        #%% mask out sea points from surface altitude
        #print(lm_data)
        or_data[lm_data == 0] = np.nan
        #print(or_data)
        sea_points = (lm_data -1)*(-1)
        
        
        
        #%% configure orogcube for xsection
        
        print(xorog)
        xorog_max = orogcube[lower_idx:upper_idx].collapsed('grid_latitude', iris.analysis.MAX)
        
        #%% locator plot
        fig = plt.figure(figsize = (12,12))
                
        gs = gridspec.GridSpec(nrows = 2, ncols = 1, height_ratios = [1.0,0.2])
        ax0 = fig.add_subplot(gs[0,0])
        ax0.pcolormesh(gridlon, gridlat, sea_points, cmap = 'GnBu', shading = 'auto')
        ax0.pcolormesh(gridlon, gridlat, or_data,  cmap = 'summer', shading = 'auto')
        ax0.scatter(gridlon, monolat, marker = '.',s = 10,color = 'r')
        ax0.scatter(gridlon, lower_monolat, marker = '.', s = 10, color = 'grey')
        ax0.scatter(gridlon, upper_monolat, marker = '.', s= 10, color = 'grey')
        ax0.set_xticks([])
        ax0.set_yticks([])
        qlons = ax0.contour(gridlon, gridlat, lons, colors = 'k')
        qlats = ax0.contour(gridlon, gridlat, lats, colors = 'k')
        label_keys = {'fontsize' : 10, 'inline' : 1,'fmt' : '%1.1f'}
        ax0.clabel(qlons, **label_keys)
        ax0.clabel(qlats, **label_keys)
        
        ax1 = fig.add_subplot(gs[1,0])
        ax1.plot(xorog.coords('grid_longitude')[0].points, xorog.data, color = 'k', label = 'Cross-section')
        ax1.plot(xorog_max.coords('grid_longitude')[0].points, xorog_max.data, color = 'k', linestyle = 'dotted', label =  'Axis maximum')
        ax1.legend(fontsize = 'small')
        ax1.set_xlim(left = xorog_max.coords('grid_longitude')[0].points[0], right = xorog_max.coords('grid_longitude')[0].points[-1])
        ax1.set_xticks(gridlon[::100]); print(gridlon[::100])
        xticklabels = []
        for o in range(5):
            xticklabels.append(f'{np.abs(lons[latidx, ::100][o]):.1f}W')
        ax1.set_xticklabels(xticklabels); print(xticklabels)
        ax1.set_xlabel('Degrees Longitude')
        ax1.set_ylim(bottom = 0)
        ax1.set_yticks([500,1000,1500,2000])
        ax1.set_ylabel('Altitude [m]')
        ax0.set_title('1p5km domain cross-island x-section location')
        plt.subplots_adjust(wspace=0.06, hspace=0.05)
        plt.savefig('D:/Project/Figures/PDF/301/u-bu807/Vertical/cross_island/1p5km_crossisland_xsection_location.pdf')
        plt.savefig('D:/Project/Figures/PDF/301/u-bu807/Vertical/cross_island/1p5km_crossisland_xsection_location.pdf')
        
    #%%
    plot_xsection = True
    if plot_xsection:
        tidx = 30
        
        #%% normalisations
        
        wspmin = 0; wspmax = 22
        wmin = -3; wmax = -wmin
        Tmin = 268; Tmax = 305
        qmin = 0; qmax = 3
        
        wsp_norm = matplotlib.colors.Normalize(vmin = wspmin, vmax = wspmax)
        w_norm = matplotlib.colors.Normalize(vmin = wmin, vmax = wmax)
        theta_norm = matplotlib.colors.Normalize(vmin = Tmin, vmax = Tmax)
        q_norm = matplotlib.colors.Normalize(vmin = qmin, vmax = qmax)
        
        wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = 'Oranges')
        w_map = matplotlib.cm.ScalarMappable(norm = w_norm, cmap = 'seismic')
        theta_map = matplotlib.cm.ScalarMappable(norm = theta_norm, cmap = 'plasma')
        q_map = matplotlib.cm.ScalarMappable(norm = q_norm, cmap = 'Blues')
        
        #%%
        cube_times = ucube.coords('time')[0]
        dates = cube_times.units.num2date(cube_times.points)
        timepoints = []
        for t in range(len(dates)):
            timepoint = dates[t].strftime('%H:%M')
            timepoints.append(timepoint)
        times = np.array(timepoints)
        print(times)
        
        #for tidx in range(48):
        #%% plot the xsections
        levels = 20
        cbar_keywargs = {'orientation':'vertical'}#, 'pad' : 0.01}
        fig = plt.figure(figsize = (16,16))
        gs = gridspec.GridSpec(nrows = 4, ncols = 2, width_ratios = [1,0.04])
        gs.update(hspace = 0.05, wspace = 0.03)
        ax0 = fig.add_subplot(gs[0,0])
        iplt.contourf(wspcube[tidx,:,latidx,:], cmap = 'Oranges', levels = levels, vmin = wspmin, vmax = wspmax)
        ax0.set_ylim(top = 8000)
        ax0.set_xlim(right = 364)
        ax0.plot(gridlon, xorog.data, color = 'k')
        ax0.set_xticks([])
        ax0.set_yticks([2000, 4000, 6000, 8000])
        ax0.set_title(f'{res} Atmostate cross-island x-section {times[tidx]} UTC')
        
        ax1 = fig.add_subplot(gs[1,0])
        iplt.contourf(wcube[tidx,:,latidx,:], cmap = 'seismic', vmin = wmin, vmax = wmax, levels = levels)
        ax1.set_ylim(top = 8000)
        ax1.set_xlim(right = 364)
        ax1.plot(gridlon, xorog.data, color = 'k')
        ax1.set_xticks([])
        ax1.set_yticks([2000, 4000, 6000, 8000])
        
        ax2 = fig.add_subplot(gs[2,0])
        iplt.contourf(thetacube[tidx,:,latidx,:], cmap = 'plasma', levels = levels, vmin = Tmin, vmax = Tmax)
        ax2.set_ylim(top = 8000)
        ax2.set_xlim(right = 364)
        ax2.plot(gridlon, xorog.data, color = 'k')
        ax2.set_xticks([])
        ax2.set_yticks([2000, 4000, 6000, 8000])
        
        ax3 = fig.add_subplot(gs[3,0])
        iplt.contourf(qcube[tidx,:,latidx,:]*1000, cmap = 'Blues', levels = levels , vmin = qmin, vmax = qmax)
        ax3.set_ylim(top = 8000)
        ax3.set_xlim(right = 364)
        ax3.plot(gridlon, xorog.data, color = 'k')
        xticklabels = []
        for o in range(5):
            xticklabels.append(f'{np.abs(lons[latidx, ::100][o]):.1f}W')
        ax3.set_xticks(gridlon[::100])
        ax3.set_xticklabels(xticklabels)
        ax3.set_yticks([2000, 4000, 6000, 8000])
        ax3.set_ylabel('Degrees longitude')
        
        ## colorbars
        cax0 = fig.add_subplot(gs[0,1])
        cbar0 = plt.colorbar(wsp_map, cax = cax0, **cbar_keywargs)
        cbar0.ax.set_ylabel('wsp $ms^{-1}$')
        cax1 = fig.add_subplot(gs[1,1])
        cbar1 = plt.colorbar(w_map, cax = cax1, **cbar_keywargs)
        cbar1.ax.set_ylabel('w $ms^{-1}$')
        cax2 = fig.add_subplot(gs[2,1])
        cbar2 = plt.colorbar(theta_map, cax = cax2, **cbar_keywargs)
        cbar2.ax.set_ylabel('$\theta$ $K$')
        cax3 = fig.add_subplot(gs[3,1])
        cbar3 = plt.colorbar(q_map, cax = cax3, **cbar_keywargs)
        cbar3.ax.set_ylabel('q $gkg^{-1}$')
        print('Complete')