# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 14:05:22 2020

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
import pandas as pd
from iris.analysis.cartography import rotate_pole
from iris.analysis import trajectory
from model_foundry import toolkit
import datetime



#%% 
def dist_tick_generator(maincube, axis, ticknumber = 6, decprecision = 2, flat = False):
    ## calculate distance from coordiate differences using haversine formula
    if flat == False:
    
        if axis == 'grid_longitude':
            lats = maincube.coord('grid_latitude').points
            lons = maincube.coord('grid_longitude').points
            # print('lenlats', len(lats))
            # print('lenlons', len(lons))
            pointnum = len(lats)
            if ticknumber < pointnum:
                seperation = pointnum // ticknumber
            else: ticknumber = pointnum; seperation = pointnum // ticknumber
            #lons = np.arange(0,720)
            dist_steps = [0]
            for i in range(pointnum-1):
                dist = toolkit.haversine(lats[i],lons[i],lats[i+1],lons[i+1])
                dist_steps.append(dist)
            dist_cum = np.cumsum(dist_steps)
            return np.sort(lons)[::seperation], ((10*dist_cum[::seperation])//1)/10
        if axis == 'grid_latitude':
            lats = maincube.coord('grid_latitude').points
            lons = maincube.coord('grid_longitude').points
            pointnum = len(lons)
            if ticknumber < pointnum:
                seperation = pointnum // ticknumber
            else: ticknumber = pointnum; seperation = pointnum // ticknumber
            #lons = np.arange(0,720)
            dist_steps = [0]
            for i in range(pointnum-1):
                dist = toolkit.haversine(lats[i],lons[i],lats[i+1],lons[i+1])
                dist_steps.append(dist)
            dist_cum = np.cumsum(dist_steps)
            return np.sort(lats)[::seperation], ((10*dist_cum[::seperation])//1)/10
    if flat == True:
        if axis == 'grid_longitude':
            lats = maincube.coord('grid_latitude').points
            lons = maincube.coord('grid_longitude').points
            # print('lenlats', len(lats))
            # print('lenlons', len(lons))
            pointnum = len(lons)
            mono_coor = np.ones(pointnum)*lats[0]
            if ticknumber < pointnum:
                seperation = pointnum // ticknumber
            else: ticknumber = pointnum; seperation = pointnum // ticknumber
            #lons = np.arange(0,720)
            dist_steps = [0]
            for i in range(pointnum-1):
                dist = toolkit.haversine(mono_coor[i],lons[i],mono_coor[i+1],lons[i+1])
                dist_steps.append(dist)
            dist_cum = np.cumsum(dist_steps)
            return np.sort(lons)[::seperation], ((10*dist_cum[::seperation])//1)/10
        if axis == 'grid_latitude':
            lats = maincube.coord('grid_latitude').points
            lons = maincube.coord('grid_longitude').points
            pointnum = len(lats)
            mono_coor = np.ones(pointnum)*lons[0]
            if ticknumber < pointnum:
                seperation = pointnum // ticknumber
            else: ticknumber = pointnum; seperation = pointnum // ticknumber
            #lons = np.arange(0,720)
            dist_steps = [0]
            for i in range(pointnum-1):
                dist = toolkit.haversine(lats[i],mono_coor[i],lats[i+1],mono_coor[i+1])
                dist_steps.append(dist)
            dist_cum = np.cumsum(dist_steps)
            return np.sort(lats)[::seperation], ((10*dist_cum[::seperation])//1)/10
        
def obs_dist_tick_gen(lons, lats, axis, ticknumber = 6, decprecision = 2):
    if axis == 'grid_longitude':
            
            # print('lenlats', len(lats))
            # print('lenlons', len(lons))
            pointnum = len(lats)
            if ticknumber < pointnum:
                seperation = pointnum // ticknumber
            else: ticknumber = pointnum; seperation = pointnum // ticknumber
            #lons = np.arange(0,720)
            dist_steps = [0]
            for i in range(pointnum-1):
                dist = toolkit.haversine(lats[i],lons[i],lats[i+1],lons[i+1])
                dist_steps.append(dist)
            dist_cum = np.cumsum(dist_steps)
            return np.sort(lons)[::seperation], ((10*dist_cum[::seperation])//1)/10
    if axis == 'grid_latitude':
            
            pointnum = len(lons)
            if ticknumber < pointnum:
                seperation = pointnum // ticknumber
            else: ticknumber = pointnum; seperation = pointnum // ticknumber
            #lons = np.arange(0,720)
            dist_steps = [0]
            for i in range(pointnum-1):
                dist = toolkit.haversine(lats[i],lons[i],lats[i+1],lons[i+1])
                dist_steps.append(dist)
            dist_cum = np.cumsum(dist_steps)
            return np.sort(lats)[::seperation], ((10*dist_cum[::seperation])//1)/10
#%%
def plot_vertical_xsection(maincube, thetacube, orogcube,index, Case, axes = ('grid_latitude', 'grid_longitude'), vmin = -250, vmax = 250, cmap = 'viridis', levels = 10, figsize = (22,8),
                           theta_contour = True, theta_color = 'g', alt_limit = 1500, unit = None):
    ## extract Case information
    res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
    fileleg = Case['fileleg']
    print(Case.values())
    
    lon_ticks, dist_cum = dist_tick_generator(maincube, axis = axes[1])
    ## extract time coords as datetime
    time = maincube.coord('time')
    dates = time.units.num2date(time.points)
    ## collapse orography across latitude
    orogmax = orogcube.collapsed(axes[0], iris.analysis.MAX)
    
    ## create list of values for colorbar
    cbar_values = np.linspace(vmin,vmax,10)
    ## plot figure
    matplotlib.rcParams.update({'font.size' : 18})
    cont_keywargs = {'coords' : [axes[1],'altitude'], 'levels' : levels, 'vmin' : vmin, 'vmax' : vmax, 'cmap' : cmap}#'figsize' : (25,10)}
    fig = plt.figure(figsize = figsize)
    q0 = iplt.contourf(maincube[index], **cont_keywargs)
    iplt.plot(maincube.coord('surface_altitude'))
    plt.xticks(ticks = lon_ticks, labels = dist_cum)
    if unit == None:
        cbar = plt.colorbar(q0, pad = 0.005, boundaries= [vmin,vmax]); cbar.ax.set(ylabel = maincube.units)# r'Heat flux [Wm$^{-2}$]')
    else:
        cbar = plt.colorbar(q0, pad = 0.005, boundaries= [vmin,vmax]); cbar.ax.set(ylabel = unit)
    q1 = iplt.plot(orogmax, color = 'k', linestyle = 'dashed')
    
    if theta_contour == True:
        thetalevels = np.arange(100,350, 2)
        q2 = iplt.contour(thetacube[index], levels = thetalevels, coords = (axes[1], 'altitude'), colors = theta_color)
    plt.ylim(top = alt_limit, bottom = 0)
    plt.xlim(left = np.min(maincube.coord(axes[1]).points), right = np.max(maincube.coord(axes[1]).points))
    plt.title(f'{res} {varname} at {dates[index]} over leg {fileleg}', size = 18)
    plt.ylabel('Altitude [m]'); plt.xlabel('Horizontal distance [km]')
    #plt.tight_layout()
    return fig

def save_vert_xsection(fig, maincube, index, pngpath, Case, pdfpath = None, fileleg = 1, obs = False, obs_detail = None, xaxis = 'lon'):
    ## extract Case information
    res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
    timecoord = maincube.coord('time')
    datetimes = timecoord.units.num2date(timecoord.points)
    timepoint = datetimes[index].strftime('H%HM%M')
    print(timepoint)
    if obs:
        figname = f'{exp}_{res}_x{xaxis}_{varname}_vert_xsect_T{timepoint}_leg{fileleg}flt{flight}_obs{obs_detail}'
    else:
        figname = f'{exp}_{res}_x{xaxis}_{varname}_vert_xsect_T{timepoint}_leg{fileleg}flt{flight}'
    ## save figure
    try:
        fig.savefig(f'{pngpath}/{figname}.png')
    except: 
        print('PNG unsuccessful, attempting to save to Dump/{figname}.png')
        try:
            fig.savefig(f'Dump/{figname}.png')
        except: 
            print('PNG failed.')
    if pdfpath != None:
        try:
            fig.savefig(f'{pdfpath}/{figname}.pdf')
        except: 
            print('PDF unsuccessful')
            pass
        
def mplt_symmetric_normalisation(data, cmap = 'bwr', limit_type = 'maximum', percentile = 0.99):
    """
    

    Parameters
    ----------
    data : an ndarray of floats; your data.
    cmap : string or matplotlib cmap object. The default is 'bwr'.
    limit_type : string: 'maximum' r 'percentile'. If percentile is sued, check which is to be used (see argument below). The default is 'maximum'.
    percentile : float, the percentile limit t be used when limit_type 'percentile' is used.
         The default is 0.99.

    Returns
    -------
    dnorm : matplotlib norm object.
    dmap : matplotlib mappable.

    """
    
    if limit_type == 'maximum':
        ## find absolute maximum
        absolute_maximum = np.nanmax(np.absolute(data))
        dnorm = matplotlib.colors.Normalize(vmin = -absolute_maximum, vmax = absolute_maximum) # "data-norm"
        dmap = matplotlib.cm.ScalarMappable(norm = dnorm, cmap = cmap) # "data-map"
    elif limit_type == 'percentile':
        absolute_percentile = np.percentile(np.absolute(data), percentile)
        dnorm = matplotlib.colors.Normalize(vmin = -absolute_percentile, vmax = absolute_percentile) # "data-norm"
        dmap = matplotlib.cm.ScalarMappable(norm = dnorm, cmap = cmap) # "data-map"
    return dnorm, dmap