# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 14:19:08 2020

@author: kse18nru
"""

#%% Import modules and scripts
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import iris
import iris.coord_categorisation
import iris.plot as iplt
import model_foundry as foundry



#%% global definitions
res = '0p5km'
variable = 'surface_altitude'
Domain = 'fulldomain'
flight = 306
suite = 'u-cf117'
experiment = 'LONGTAIL'
config = ''
plt.close()
#try:
ocube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{experiment}_{res}_umpa1_flt{flight}.nc', variable)
lmcube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{experiment}_{res}_umpa1_flt{flight}.nc', 'land_binary_mask')
#except:
#    ocube = iris.load_cube(f'../../Model_Data/{suite}/nc/{experiment}/{variable}/{res}_{variable}_301.nc', variable)
#    lmcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{experiment}/land_binary_mask/{res}_land_binary_mask_301.nc')
lmdata = lmcube.data
lmlon = lmcube.coord('grid_longitude').points
lmlat = lmcube.coord('grid_latitude').points
print(np.shape(ocube.data))
odata = ocube.data
olat = ocube.coord('grid_latitude').points
olon = ocube.coord('grid_longitude').points

polelat = ocube.coord('grid_latitude').coord_system.grid_north_pole_latitude
polelon = ocube.coord('grid_longitude').coord_system.grid_north_pole_longitude

#%%
if Domain == 'fulldomain':
    
    lons, lats = foundry.modf.unrotate_coords(olon, olat, polelon, polelat)
    ## filter out sea points
    odata[np.where(lmdata != 1)] = np.nan
    ## define contour levels for each resolution
    if res == '0p5km':
        levels = [0, 30,60,90,120,150,180,210,240,270,300,350,400,450,500,600,700,800,900,1000]
        levels = [0,500,1000,1500]
        figsize = (10,14)
    elif res == '1p5km':
        levels = [0,50,100,200,300,500,600, 700,800,900,1000,1100,1200,1400,1600,1800]
        figsize = (17,16)
    elif res == '4p4km':
        levels = [0,50,100,200,300,500,600,800,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000]
        figsize = (17,16)
    
    fig, ax = plt.subplots(1,1, figsize = figsize)
    #q0 = ax.pcolormesh(olon,olat,odata)
    cmap = matplotlib.cm.terrain
    norm = matplotlib.colors.Normalize(vmin=-400, vmax = 1500)
    q0 = ax.pcolormesh(olon, olat, odata, cmap = cmap, norm = norm)
    #q0 = ax.contour(olon, olat, odata, colors = 'k', levels = levels )
    #q0 = ax.contour(olon, olat, odata, colors = 'k')
    q1a = ax.contour(olon,olat,lons, colors = 'k', linestyles = 'dotted')
    q1b = ax.contour(olon,olat,lats, colors = 'k', linestyles = 'dotted')
    
    label_keys = {'fontsize' : 10, 'inline' : 1}
    #q2a = ax.clabel(q0, **label_keys)
    # colourbar
    cbar0 = fig.colorbar(q0, ax = ax)
    matplotlib.rcParams['contour.negative_linestyle'] = 'dotted' # make contours for negative values solid, not dashed
    q2b = ax.clabel(q1a, **label_keys)
    q2c = ax.clabel(q1b, **label_keys)
    # hide axis labels
    ql = ax.set_yticklabels([]); q9 = ax.set_xticklabels([])
    
    # set title
    qt = ax.set(Title = f'Model orography for {res} with e = 1 smoothing')
    
    # save figure
    #plt.savefig(f'D:/Project/Figures/PDF/{flight}/{suite}/{Case}_surface_altitude_{Domain}_pcolor_{experiment}_cc134_flt{flight}.pdf')
    plt.savefig(f'D:/Project/Figures/PNG/{flight}/{experiment}/{suite}/{res}_surface_altitude_{Domain}_pcolor_{experiment}_cf117_flt{flight}.png')
    plt.show()
    
    
#%% Subset to fit peninsulas domain
if Domain == 'peninsulas':
    ## boundaries (in points):
    if res == '4p4km':
        South = 115; North = 200; West = 150; East = 185
    if res == '1p5km':
        South = 110; North = 350; West = 90; East = 185
    if res == '0p5km':
        South = 80; North = -1; West = 50; East = 350
    # subset 4d data
    orog = foundry.modf.subset_2d(odata, South, North, West, East)
    
           
    # two-dimensional subsetting of land-mask
    lm = foundry.modf.subset_2d(lmdata, South, North, West, East )
    
    ## filter out sea points
    orog[np.where(lm != 1)] = np.nan
    # subset coordinates and unrotate
    sub_olon = foundry.modf.subset_1d(olon, West, East)
    sub_olat = foundry.modf.subset_1d(olat, South, North)
    sub_lmlon = foundry.modf.subset_1d(lmlon, West, East)
    sub_lmlat = foundry.modf.subset_1d(lmlat, South, North)
    #%% unrotate pole
    lons, lats = foundry.modf.unrotate_coords(sub_olon, sub_olat, polelon, polelat)
    # make list for lons lats levels to be plotted
    lats_levels = np.arange(50,70,0.5); print(lats_levels)
    lons_levels = np.sort(-np.arange(10,30,0.5)); print(lons_levels)
    #%% plotting
    sub_figsize = (16,20)
    fig, ax = plt.subplots(1,1,figsize = sub_figsize)
    
    cmap = matplotlib.cm.terrain
    norm = matplotlib.colors.Normalize(vmin=-400, vmax = 1500)
    q0 = ax.contourf(sub_olon, sub_olat, orog, cmap = cmap, norm = norm, levels = 20)
    cbar = fig.colorbar(q0, ax = ax)
    #q1 = ax.contour(sub_lmlon, sub_lmlat, lm, colors = 'k', levels = 1)
    # set title
    qt = ax.set(Title = f'Model orography for {res} with e = 1 smoothing')
    # hide axis labels
    ql = ax.set_yticklabels([]); q9 = ax.set_xticklabels([])
    # save figure
    plt.savefig(f'../../Figures/PDF/301/{res}_surface_altitude_{Domain}_contf_{experiment}_bu807_flt301.pdf')
    plt.savefig(f'../../Figures/PNG/301/{res}_surface_altitude_{Domain}_contf_{experiment}_bu807_flt301.png')
    
    plt.show()

#%% Subset to fit snaesfellness domain
if Domain == 'snaesfellness':
    ## boundaries (in points):
    if res == '4p4km':
        South = 143; North = 155; West = 150; East = 178
    if res == '1p5km':
        South = 190; North = 220; West = 90; East = 160
    if res == '0p5km':
        South = 280; North = 370; West = 50; East = 250
    # subset 4d data
    orog = foundry.modf.subset_2d(odata, South, North, West, East)
    
           
    # two-dimensional subsetting of land-mask
    lm = foundry.modf.subset_2d(lmdata, South, North, West, East )
    
    ## filter out sea points
    #orog_spare = orog
    #orog_spare[np.where(orog_spare == np.nan)] = 0
    orog[np.where(lm != 1)] = np.nan
    # subset coordinates and unrotate
    sub_olon = foundry.modf.subset_1d(olon, West, East)
    sub_olat = foundry.modf.subset_1d(olat, South, North)
    sub_lmlon = foundry.modf.subset_1d(lmlon, West, East)
    sub_lmlat = foundry.modf.subset_1d(lmlat, South, North)
    #%% unrotate pole
    lons, lats = foundry.modf.unrotate_coords(sub_olon, sub_olat, polelon, polelat)
    # make list for lons lats levels to be plotted
    lats_levels = np.arange(50,70,0.25); print(lats_levels)
    lons_levels = np.sort(-np.arange(10,30,0.25)); print(lons_levels)
    label_keys = {'fontsize' : 10, 'inline' : 1}
    #%% take cross-section values
    orog_mean = np.nanmean(orog, axis = 0)
    print(np.isnan(orog_mean))
    orog_mean[np.isnan(orog_mean)] = 0
    print(orog_mean)
    print('hibijibees')
    orog_max = np.nanmax(orog, axis = 0)
    
    orog_max[np.isnan(orog_max)] = 0
    #orog_min = np.nanmin(orog, axis = 1)
    
    #%% plotting
    sub_figsize = (16,12)
    fig, axes = plt.subplots(2,1,figsize = sub_figsize)
    
    cmap = matplotlib.cm.terrain
    norm = matplotlib.colors.Normalize(vmin=-400, vmax = 1500)
    q0 = axes[0].contourf(sub_olon, sub_olat, orog, cmap = cmap, norm = norm, levels = 20)
    cbar = fig.colorbar(q0, ax = axes[0],  orientation = 'horizontal', label = 'Surface altitude [m]', pad = 0.1)
    q1a = axes[0].contour(sub_olon,sub_olat,lons, colors = 'k', linestyles = 'dotted', levels = lons_levels)
    q1b = axes[0].contour(sub_olon,sub_olat,lats, colors = 'k', linestyles = 'dotted', levels = lats_levels)
    q2b = axes[0].clabel(q1a, **label_keys)
    q2c = axes[0].clabel(q1b, **label_keys)
    q1 = axes[1].plot(orog_mean, color = 'k', label = 'Mean land point altitude \nalong axis')
    #q2 = axes[1].plot(orog_min, linestyle = ':', color = 'k')
    q3 = axes[1].plot(orog_max, linestyle = '--', color = 'k', label = 'Max. altitude along axis')
    axes[1].legend()
    #q1 = ax.contour(sub_lmlon, sub_lmlat, lm, colors = 'k', levels = 1)
    # set title
    qt = axes[0].set(Title = f'Model orography for {res} with e = 1 smoothing')
    # hide axis labels
    axes[1].set(ylabel = 'Surface altitude [m]')
    axes[0].set_yticklabels([]); axes[0].set_xticklabels([]); 
    axes[1].set_xticklabels([]);
    axes[1].set_ylim(top = 1400)
    plt.tight_layout()
    # save figure
    #plt.savefig(f'../../Figures/PDF/301/{Case}_surface_altitude_{Domain}_contf_{experiment}_bu807_flt301.pdf')
    #plt.savefig(f'../../Figures/PNG/301/{Case}_surface_altitude_{Domain}_contf_{experiment}_bu807_flt301.png')
    
    plt.show()

if Domain == 'xsections':
   
    
    #%% Case 4p4km
    South = 143; North = 155; West = 150; East = 178
    # subset 4d data
    ocube = iris.load_cube(f'../../Model_Data/{suite}/nc/{experiment}/{variable}/4p4km_{variable}_flt301.nc', variable)
    odata = ocube.data
    lmcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{experiment}/land_binary_mask/1p5km_land_binary_mask_flt301.nc')
    lmdata = lmcube.data
    olat = ocube.coord('grid_latitude').points
    olon = ocube.coord('grid_longitude').points
    orog = foundry.modf.subset_2d(odata, South, North, West, East)
    
           
    # two-dimensional subsetting of land-mask
    lm = foundry.modf.subset_2d(lmdata, South, North, West, East )
    
    ## filter out sea points
    #orog_spare = orog
    #orog_spare[np.where(orog_spare == np.nan)] = 0
    orog[np.where(lm != 1)] = np.nan
    # subset coordinates and unrotate
    sub_olon = foundry.modf.subset_1d(olon, West, East)
    sub_olat = foundry.modf.subset_1d(olat, South, North)
    
    
    
    
    
    orog_mean = np.nanmean(orog, axis = 0)
    
    orog_mean[np.isnan(orog_mean)] = 0
    
    orog_max = np.nanmax(orog, axis = 0)
    
    orog_max[np.isnan(orog_max)] = 0
    lon_4p4 = sub_olon
    mean_4p4 = orog_mean
    max_4p4 = orog_max
    
    #%% Case 1p5km
    South = 190; North = 220; West = 90; East = 160
    # subset 4d data
    ocube = iris.load_cube(f'../../Model_Data/{suite}/nc/{experiment}/{variable}/1p5km_{variable}_flt301.nc', variable)
    odata = ocube.data
    lmcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{experiment}/land_binary_mask/1p5km_land_binary_mask_flt301.nc')
    lmdata = lmcube.data
    olat = ocube.coord('grid_latitude').points
    olon = ocube.coord('grid_longitude').points
    orog = foundry.modf.subset_2d(odata, South, North, West, East)
    
           
    # two-dimensional subsetting of land-mask
    lm = foundry.modf.subset_2d(lmdata, South, North, West, East )
    
    ## filter out sea points
    #orog_spare = orog
    #orog_spare[np.where(orog_spare == np.nan)] = 0
    print('!!!',np.mean(orog))
    orog[np.where(lm != 1)] = np.nan
    # subset coordinates and unrotate
    sub_olon = foundry.modf.subset_1d(olon, West, East)
    sub_olat = foundry.modf.subset_1d(olat, South, North)
    
    
    
    
    
    
    orog_mean = np.nanmean(orog, axis = 0)
    
    orog_mean[np.isnan(orog_mean)] = 0
    
    orog_max = np.nanmax(orog, axis = 0)
    
    orog_max[np.isnan(orog_max)] = 0
    
    lon_1p5 = sub_olon
    mean_1p5 = orog_mean
    max_1p5 = orog_max
    print(odata)
    
    #%% Case 0p5km
    
    South = 280; North = 370; West = 50; East = 250
    # subset 4d data
    ocube = iris.load_cube(f'../../Model_Data/{suite}/nc/{experiment}/{variable}/0p5km_{variable}_flt301.nc', variable)
    lmcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{experiment}/land_binary_mask/0p5km_land_binary_mask_flt301.nc')
    odata = ocube.data
    lmdata = lmcube.data
    olat = ocube.coord('grid_latitude').points
    olon = ocube.coord('grid_longitude').points
    orog = foundry.modf.subset_2d(odata, South, North, West, East)
    
           
    # two-dimensional subsetting of land-mask
    lm = foundry.modf.subset_2d(lmdata, South, North, West, East )
    
    ## filter out sea points
    #orog_spare = orog
    #orog_spare[np.where(orog_spare == np.nan)] = 0
    orog[np.where(lm != 1)] = np.nan
    # subset coordinates and unrotate
    sub_olon = foundry.modf.subset_1d(olon, West, East)
    sub_olat = foundry.modf.subset_1d(olat, South, North)
    
    
    
    
    orog_mean = np.nanmean(orog, axis = 0)
    
    orog_mean[np.isnan(orog_mean)] = 0
    
    orog_max = np.nanmax(orog, axis = 0)
    
    orog_max[np.isnan(orog_max)] = 0
    
    lon_0p5 = sub_olon
    mean_0p5 = orog_mean
    max_0p5 = orog_max
    
    #%% plot figure
    fig, ax  = plt.subplots(1,1, figsize = (18,10))
    ax.plot(lon_0p5, mean_0p5, color = 'b', label = '0p5km mean')
    ax.plot(lon_0p5, max_0p5, color = 'b', linestyle = '--', label = '0p5km maximum')
    ax.plot(lon_1p5, mean_1p5, color = 'r', label = '1p5km mean')
    ax.plot(lon_1p5, max_1p5, color = 'r', linestyle = '--', label = '1p5km maximum')
    ax.plot(lon_4p4, mean_4p4, color = 'g', label = '4p4km mean')
    ax.plot(lon_4p4, max_4p4, color = 'g', linestyle = '--',label = '4p4km maximum')
    plt.legend()
    ax.set(ylabel = 'Surface altitude [m]', xlabel = r'Model rotated grid x coordinate [$^{\circ}$]')
    plt.savefig(f'../../Figures/PDF/301/Allres_surface_altitude_{Domain}_{experiment}_bu807_flt301.pdf')
    plt.savefig(f'../../Figures/PNG/301/Allres_surface_altitude_{Domain}_{experiment}_bu807_flt301.png')
    plt.show()