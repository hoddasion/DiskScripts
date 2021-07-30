# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 15:02:39 2020

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
def dist_tick_generator(maincube, axis, ticknumber = 6, decprecision = 2):
    ## calculate distance from coordiate differences using haversine formula
    if axis == 'grid_longitude':
        lats = maincube.coord('grid_latitude').points
        lons = maincube.coord('grid_longitude').points
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
    
#%%
def plot_vertical_xsection(maincube, thetacube, orogcube,index, Case, axes = ('grid_latitude', 'grid_longitude'), vmin = -250, vmax = 250, cmap = 'viridis', levels = 10, figsize = (22,8)):
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
    
    ## plot figure
    matplotlib.rcParams.update({'font.size' : 18})
    cont_keywargs = {'coords' : [axes[1],'altitude'], 'levels' : levels, 'vmin' : vmin, 'vmax' : vmax, 'cmap' : cmap}#'figsize' : (25,10)}
    fig = plt.figure(figsize = figsize)
    q0 = iplt.contourf(maincube[index], **cont_keywargs)
    plt.xticks(ticks = lon_ticks, labels = dist_cum)
    cbar = plt.colorbar(q0, pad = 0.005); cbar.ax.set(ylabel = maincube.units)# r'Heat flux [Wm$^{-2}$]')
    q1 = iplt.plot(orogmax, color = 'k', linestyle = 'dashed')
    thetalevels = np.arange(100,350, 2)
    q2 = iplt.contour(thetacube[index], levels = thetalevels, coords = (axes[1], 'altitude'), color = 'g')
    plt.ylim(top = 1500, bottom = 0)
    plt.xlim(left = np.min(maincube.coord(axes[1]).points), right = np.max(maincube.coord(axes[1]).points))
    plt.title(f'{varname} at {dates[index]} over leg {fileleg}', size = 18)
    plt.ylabel('Altitude [m]'); plt.xlabel('Horizontal distance [km]')
    #plt.tight_layout()
    return fig

def save_vert_xsection(fig, maincube, index, pngpath, Case, pdfpath = None, fileleg = 1):
    ## extract Case information
    res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
    timecoord = maincube.coord('time')
    datetimes = timecoord.units.num2date(timecoord.points)
    timepoint = datetimes[index].strftime('H%HM%M')
    print(timepoint)
    figname = f'{exp}_{res}_{varname}_vert_xsect_T{timepoint}_leg{fileleg}flt{flight}'
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
#%% execute function
legs = [1,2,3]
resolutions = ['0p5km', '1p5km', '4p4km']

for D in resolutions:
    for leg in legs:
        try:
            #%%
        
            Case = {'res' : D, 'flight' : 301, 'varname' : 'upward_heat_flux_in_air', 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
            res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
            fileleg = Case['fileleg']
            print(Case.values())
            #%%
            
            maincube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', varname)[:,:23]
            thetacube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
            orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
            print(maincube, '\n', thetacube, '\n', orogcube)
            print(maincube.units)
        except:
            print('File error')
            continue
        path = f'../../Figures/PNG/301/u-bu807/Vertical/{varname}/leg{fileleg}/{res}'
        for index in range(47):
            
    
            
            fig = plot_vertical_xsection(maincube, thetacube, orogcube,index, Case, cmap = 'seismic', vmin = -400,vmax = 400, levels = 20, alt_limit = 2000)
            save_vert_xsection(fig, maincube,index, path, Case, fileleg = leg)
            plt.close()
