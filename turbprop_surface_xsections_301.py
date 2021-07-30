# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 12:27:55 2021

@author: kse18nru
"""

#%%

import iris
import numpy as np
import iris.plot as iplt
import matplotlib.pyplot as plt
#%%
for res in ['0p5km']:#,'1p5km','4p4km']:
    tidx = 24
    directory = 'D:/Project/Model_Data/u-bu807/nc/Control/'
    lhfile = f'RA1M_{res}_um_surface_upward_latent_heat_flux_24hrs_g_301.nc'
    shfile = f'RA1M_{res}_um_surface_upward_sensible_heat_flux_24hrs_g_301.nc'
    eaststressfile = f'RA1M_{res}_um_surface_downward_eastward_stress_24hrs_g_301.nc'
    northstressfile = f'RA1M_{res}_um_surface_downward_northward_stress_24hrs_g_301.nc'
    lhcube = iris.load_cube(f'{directory}{lhfile}', 'surface_upward_latent_heat_flux')
    landcube = iris.load_cube(f'{directory}land_binary_mask/{res}_land_binary_mask_flt301.nc', 'land_binary_mask')
    shcube = iris.load_cube(f'{directory}{shfile}', 'surface_upward_sensible_heat_flux')
    estresscube = iris.load_cube(f'{directory}{eaststressfile}', 'surface_downward_eastward_stress')
    nstresscube = iris.load_cube(f'{directory}{northstressfile}', 'surface_downward_northward_stress')
    ufile = f'RA1M_{res}_um_x_wind_24hrs_g_301.nc'
    vfile = f'RA1M_{res}_um_y_wind_24hrs_g_301.nc'
    ucube = iris.load_cube(f'{directory}{ufile}', 'x_wind')
    vcube = iris.load_cube(f'{directory}{vfile}', 'y_wind')
    mspfile = f'RA1M_{res}_um_air_pressure_at_sea_level_24hrs_g_301.nc'
    mspcube = iris.load_cube(f'{directory}{mspfile}', 'air_pressure_at_sea_level')
    
    #print(ucube, vcube[:,:-1])
    
    #%% extract data arrays
    lhdata = lhcube.data
    shdata = shcube.data
    estress = estresscube.data
    nstress = nstresscube.data[:,:-1,:]
    stressdata = (estress**2+nstress**2)**0.5
    glon = lhcube.coords('grid_longitude')[0].points
    glat = lhcube.coords('grid_latitude')[0].points
    wspdata = (ucube.data**2 + vcube.data[:,:-1]**2)**0.5
    mspdata = mspcube.data
    #msplevels = 
    
    landmask = landcube.data
    landlon = landcube.coords('grid_longitude')[0].points
    landlat = landcube.coords('grid_latitude')[0].points
    #%% dynamically adjust figure size
    size_ratio = len(glat)/len(glon)
    print(size_ratio)
    base_length = 10
    xscale = base_length
    yscale = base_length*size_ratio
    figsize = (xscale, yscale)
    
    for tidx in [20,24]:
        #%%
        fig, ((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2,figsize = figsize)
        
        ax0.pcolormesh(glon, glat, stressdata[tidx],cmap = 'Blues', vmin = 0, shading = 'auto')
        ax0.contour(landlon, landlat, landmask, colors = 'k')
        
        
        ax1.pcolormesh(glon, glat, shdata[tidx], shading = 'auto', cmap = 'seismic', vmin = -200, vmax = 200)
        ax1.contour(landlon, landlat, landmask, colors = 'k')
        
        ax2.pcolormesh(glon, glat,lhdata[tidx], cmap = 'seismic', vmin = -200, vmax = 200, shading = 'auto')
        ax2.contour(landlon, landlat, landmask, colors = 'k')
        
        ax3.pcolormesh(glon, glat, wspdata[tidx], cmap = 'Oranges', vmin = 0, shading = 'auto')
        ax3.contour(landlon, landlat, landmask, colors = 'k')