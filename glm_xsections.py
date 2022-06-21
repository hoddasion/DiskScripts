# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 12:13:14 2021

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
import iris.analysis.cartography
import cartopy.crs as ccrs

#%% global variables
plt.rcParams.update({'font.size': 22})
flight = 306
suite = 'u-cc134'
variable = ''
kconstraints = {'longitude' :(-60,20), 'latitude' :(35,85)} # coordinate constraints for iris intersection function

def glm_load_and_concat(stashcode, kconstraints, fileno = 4, flight = 306, suite = 'u-cc134'):
    if fileno < 2:
        return iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_glm_{fileno}_flt{flight}.nc', stashcode).intersection(**kconstraints)
    else:   
       
        files = []
        for i in range(fileno - 1):
            filename = f'D:/Project/Model_Data/{suite}/RA1M_glm_{i+2}_flt{flight}.nc'
            files.append(filename)
        
        #from iris.experimental.equalise_cubes import equalise_attributes
        cubes =  iris.load(files, [stashcode])
        iris.util.unify_time_units(cubes)
        #print(equalise_attributes(cubes))
        #cubes[0] = cubes[0][1:]
        for i in range(3):
            print(cubes[i].coords('time')[0])
        conccube = cubes.concatenate()[0].intersection(**kconstraints) # a complete cube spanning full forecast
        
        #print(iris.load(files, [stashcode]).concatenate_cube())
        return conccube

def extract_data_and_coords(cube):
    return cube.data,cube.coords('longitude')[0].points,cube.coords('latitude')[0].points
#%%load data

cubes = iris.load(f'D:/Project/Model_Data/{suite}/RA1M_glm_1_flt{flight}.nc')
print(cubes)
cubes = None
#print(cubes[17])

#%% reduced multiplot
reduced_multiplot = True
if reduced_multiplot:
    sxwind_cube = glm_load_and_concat('m01s03i225', kconstraints)
    sywind_cube = glm_load_and_concat('m01s03i226', kconstraints)
    msp_cube = glm_load_and_concat('m01s16i222', kconstraints)
    lbm_cube = glm_load_and_concat('m01s00i030', kconstraints, fileno = 1)
    
    prec_cube = glm_load_and_concat('m01s05i216', kconstraints)
    #%%
    wind_cube = (sxwind_cube**2 + sywind_cube**2)**0.5
    #%%
    print(wind_cube.coords('time')[0])
    #%%
    
    wind_data, wind_lon, wind_lat = extract_data_and_coords(wind_cube)
    lbm_data, lbm_lon, lbm_lat = extract_data_and_coords(lbm_cube)
    sx_data = sxwind_cube.data
    sy_data = sywind_cube.data
    msp_data, msp_lon, msp_lat = extract_data_and_coords(msp_cube)
   
    prec_data, prec_lon, prec_lat = extract_data_and_coords(prec_cube)
    #%%  time frmatting
    
    cube_time = wind_cube.coord('time')
    dates = cube_time.units.num2date(cube_time.points)
    timepoints = []
    for t in range(len(dates)):
        timepoint = dates[t].strftime('%d-%m-18 %H:%M')
        timepoints.append(timepoint)
    times = np.array(timepoints)
    print(times)
    print(cube_time.points)
    
    #%% find min and max across day for normalisation
    wmin = 0
    wmax = np.nanpercentile(wind_data[:4],99.9)
    mmin = np.nanpercentile(msp_data[:4], 0.01)/100
    mmax = np.nanpercentile(msp_data[:4], 99.9)/100
    
    pmin = np.nanpercentile(prec_data, 0)*1000
    pmax = 0.06#np.nanpercentile(prec_data, 90)*1000
    
    wsp_norm = matplotlib.colors.Normalize(vmin = wmin, vmax = wmax)
    wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = 'Oranges')
    msp_norm = matplotlib.colors.Normalize(vmin = mmin, vmax = mmax)
    msp_map = matplotlib.cm.ScalarMappable(norm = msp_norm, cmap = 'viridis')
   
    
    prf_norm = matplotlib.colors.Normalize(vmin = pmin, vmax = pmax)
    prf_map = matplotlib.cm.ScalarMappable(norm = prf_norm, cmap = 'Blues')
    
    
    #%% mask out low precipitation values
    print(np.nanpercentile(prec_data, 75)*1000)
    print(pmax)
    #%%
    unmasked_prec_data = np.copy(prec_data)
    #prec_data[prec_data < np.nanpercentile(prec_data, 75)] = np.nan
    
    #%% plot figure
    data_crs = ccrs.PlateCarree()
    gs = gridspec.GridSpec(nrows = 3, ncols = 2,wspace=0.0, hspace=0.0)
    projtype = ccrs.NorthPolarStereo()
    #plevels = np.arange(-10,10)*5#[-50,-40,-30,-20,-10,0,10,20,30,40,50]
    splt_kwargs = {'projection' : projtype, 'frameon' : False} # subplot keyword arguments
    wind_kwargs = {'transform':data_crs,'cmap':'Oranges','vmin':wmin,'vmax':wmax}
    
    
    prec_kwargs = {'transform':data_crs,  'cmap':'Blues', 'vmin':pmin,'vmax':pmax}
    lbm_args = [lbm_lon,lbm_lat,lbm_data] # land mask arguments
    lbm_kwargs = {'colors':'k','transform':data_crs} # land mask keyword arguments
    qui_kwargs = {'transform':data_crs} # quiver keywords
    q = 15
    fig = plt.figure(figsize = (10,14))
    
    ### create and plot windspeed axes
    ax1 = fig.add_subplot(gs[0,0],**splt_kwargs)
    #iplt.pcolormesh(wind_cube[0], cmap = 'Oranges')
    ax1.pcolormesh(wind_lon,wind_lat,wind_data[0],**wind_kwargs)
    ax1.contour(*lbm_args, **lbm_kwargs)
    ax1.quiver(wind_lon[::q], wind_lat[:-20:q], sx_data[0,:-20:q,::q],sy_data[0,:-20:q,::q], **qui_kwargs)
    
    
    ax2 = fig.add_subplot(gs[1,0], **splt_kwargs)
    ax2.pcolormesh(wind_lon, wind_lat, wind_data[2], **wind_kwargs)
    ax2.contour(*lbm_args,**lbm_kwargs)
    ax2.quiver(wind_lon[::q], wind_lat[:-20:q], sx_data[1,:-20:q,::q],sy_data[1,:-20:q,::q], **qui_kwargs)
    
    
    ax3 = fig.add_subplot(gs[2,0], **splt_kwargs)
    ax3.pcolormesh(wind_lon, wind_lat, wind_data[4], **wind_kwargs)
    ax3.contour(*lbm_args,**lbm_kwargs)
    ax3.quiver(wind_lon[::q], wind_lat[:-20:q], sx_data[2,:-20:q,::q],sy_data[2,:-20:q,::q], **qui_kwargs)
    
    
    
    
    ### create and plot msp axes
    p_dev = 1020
    plevels = np.arange(-10,10)*5
    msp_kwargs = {'transform':data_crs, 'colors':'k','levels':plevels}
    
    ax5 = fig.add_subplot(gs[0,1], **splt_kwargs)
    ax5.contour(msp_lon, msp_lat, msp_data[0]/100 -p_dev, **msp_kwargs)
    ax5.pcolormesh(prec_lon, prec_lat, prec_data[0]*1000, **prec_kwargs)
    ax5.contour(*lbm_args, **lbm_kwargs)
    
    
    ax6 = fig.add_subplot(gs[1,1], **splt_kwargs)
    ax6.contour(msp_lon, msp_lat, msp_data[2]/100 -p_dev, **msp_kwargs)
    ax6.pcolormesh(prec_lon, prec_lat, prec_data[2]*1000, **prec_kwargs)
    ax6.contour(*lbm_args, **lbm_kwargs)
    
    ax7 = fig.add_subplot(gs[2,1], **splt_kwargs)
    ax7.contour(msp_lon, msp_lat, msp_data[4]/100 - p_dev, **msp_kwargs)
    ax7.pcolormesh(prec_lon, prec_lat, prec_data[4]*1000, **prec_kwargs)
    ax7.contour(*lbm_args, **lbm_kwargs)
    
    
    
    plt.text(0.15,0.86,times[0],transform=plt.gcf().transFigure, fontsize = 27)
    plt.text(0.15,0.608,times[2],transform=plt.gcf().transFigure, fontsize = 27)
    plt.text(0.15,0.355,times[4],transform=plt.gcf().transFigure, fontsize = 27)
    
    ### colorbars
    cbar1 = plt.colorbar(wsp_map,ax = [ax1,ax2,ax3],  orientation = 'horizontal', fraction = 0.012, pad = 0.01)
    cbar1.ax.set(xlabel = r'$ms^{-1}$')
    cbar1.set_ticks([0,5,10,15,20])
    cbar2 = plt.colorbar(prf_map,ax =[ax5,ax6,ax7] ,  orientation = 'horizontal', fraction = 0.012, pad = 0.01)
    cbar2.ax.set(xlabel = r'$gm^{-2}s^{-1}$')
    cbar2.set_ticks([0,0.02, 0.04, 0.06])
    
    fig.suptitle('UM GA6.1 n768 -  Case 2', y = 0.91)
    
    #plt.savefig(f'D:/Project/Figures/PDF/{flight}/{suite}/glm/GLM_um_surf_wsp_msp_prf_2colcompr{suite[2:]}_flt{flight}.pdf')
    plt.savefig(f'D:/Project/Figures/PNG/{flight}/{suite}/glm/GLM_um_surf_wsp_msp_prf_2colcompr{suite[2:]}_flt{flight}.png')
    plt.show()
    
    
#%%
surf_wind_plot = False
if surf_wind_plot:
    sxwind_cube = glm_load_and_concat('m01s03i225', kconstraints)
    sywind_cube = glm_load_and_concat('m01s03i226', kconstraints)
    msp_cube = glm_load_and_concat('m01s16i222', kconstraints)
    lbm_cube = glm_load_and_concat('m01s00i030', kconstraints, fileno = 1)
    temp_cube = glm_load_and_concat('m01s03i236', kconstraints)
    prec_cube = glm_load_and_concat('m01s05i216', kconstraints)
    #%%
    wind_cube = (sxwind_cube**2 + sywind_cube**2)**0.5
    #%%
    print(wind_cube.coords('time')[0])
    #%%
    
    wind_data, wind_lon, wind_lat = extract_data_and_coords(wind_cube)
    lbm_data, lbm_lon, lbm_lat = extract_data_and_coords(lbm_cube)
    sx_data = sxwind_cube.data
    sy_data = sywind_cube.data
    msp_data, msp_lon, msp_lat = extract_data_and_coords(msp_cube)
    temp_data, temp_lon, temp_lat = extract_data_and_coords(temp_cube)
    prec_data, prec_lon, prec_lat = extract_data_and_coords(prec_cube)
    #%%  time frmatting
    
    cube_time = wind_cube.coord('time')
    dates = cube_time.units.num2date(cube_time.points)
    timepoints = []
    for t in range(len(dates)):
        timepoint = dates[t].strftime('%d-%m-18 %H:%M')
        timepoints.append(timepoint)
    times = np.array(timepoints)
    print(times)
    print(cube_time.points)
    
    #%% find min and max across day for normalisation
    wmin = 0
    wmax = np.nanpercentile(wind_data,99.9)
    mmin = np.nanpercentile(msp_data, 0.01)/100
    mmax = np.nanpercentile(msp_data, 99.9)/100
    tmin = np.nanpercentile(temp_data, 0.01)
    tmax = np.nanpercentile(temp_data, 99.9)
    pmin = 0
    pmax = np.nanpercentile(prec_data, 99)*1000
    
    wsp_norm = matplotlib.colors.Normalize(vmin = wmin, vmax = wmax)
    wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = 'Oranges')
    msp_norm = matplotlib.colors.Normalize(vmin = mmin, vmax = mmax)
    msp_map = matplotlib.cm.ScalarMappable(norm = msp_norm, cmap = 'viridis')
    tem_norm = matplotlib.colors.Normalize(vmin = tmin, vmax = tmax)
    tem_map = matplotlib.cm.ScalarMappable(norm = tem_norm, cmap = 'plasma')
    prf_norm = matplotlib.colors.Normalize(vmin = pmin, vmax = pmax)
    prf_map = matplotlib.cm.ScalarMappable(norm = prf_norm, cmap = 'Blues')
    #%% plot figure
    data_crs = ccrs.PlateCarree()
    gs = gridspec.GridSpec(nrows = 3, ncols = 4)
    projtype = ccrs.NorthPolarStereo()
    plevels = np.arange(-10,10)*5#[-50,-40,-30,-20,-10,0,10,20,30,40,50]
    splt_kwargs = {'projection' : projtype, 'frameon' : False} # subplot keyword arguments
    wind_kwargs = {'transform':data_crs,'cmap':'Oranges','vmin':wmin,'vmax':wmax}
    msp_kwargs = {'transform':data_crs, 'colors':'k','levels':plevels}#'cmap':'viridis','vmin':mmin, 'vmax':mmax}
    temp_kwargs = {'transform':data_crs, 'cmap':'plasma', 'vmin':tmin,'vmax':tmax}
    prec_kwargs = {'transform':data_crs, 'cmap':'Blues', 'vmin':pmin,'vmax':pmax}
    lbm_args = [lbm_lon,lbm_lat,lbm_data] # land mask arguments
    lbm_kwargs = {'colors':'k','transform':data_crs} # land mask keyword arguments
    qui_kwargs = {'transform':data_crs} # quiver keywords
    q = 15
    fig = plt.figure(figsize = (18,13))
    
    ### create and plot windspeed axes
    ax1 = fig.add_subplot(gs[0,0],**splt_kwargs)
    #iplt.pcolormesh(wind_cube[0], cmap = 'Oranges')
    ax1.pcolormesh(wind_lon,wind_lat,wind_data[0],**wind_kwargs)
    ax1.contour(*lbm_args, **lbm_kwargs)
    ax1.quiver(wind_lon[::q], wind_lat[:-20:q], sx_data[0,:-20:q,::q],sy_data[0,:-20:q,::q], **qui_kwargs)
    ax1.set_title(times[0])
    
    ax2 = fig.add_subplot(gs[0,1], **splt_kwargs)
    ax2.pcolormesh(wind_lon, wind_lat, wind_data[1], **wind_kwargs)
    ax2.contour(*lbm_args,**lbm_kwargs)
    ax2.quiver(wind_lon[::q], wind_lat[:-20:q], sx_data[1,:-20:q,::q],sy_data[1,:-20:q,::q], **qui_kwargs)
    ax2.set_title(times[1])
    
    ax3 = fig.add_subplot(gs[0,2], **splt_kwargs)
    ax3.pcolormesh(wind_lon, wind_lat, wind_data[2], **wind_kwargs)
    ax3.contour(*lbm_args,**lbm_kwargs)
    ax3.quiver(wind_lon[::q], wind_lat[:-20:q], sx_data[2,:-20:q,::q],sy_data[2,:-20:q,::q], **qui_kwargs)
    ax3.set_title(times[2])
    
    ax4 = fig.add_subplot(gs[0,3], **splt_kwargs)
    ax4.pcolormesh(wind_lon, wind_lat, wind_data[3], **wind_kwargs)
    ax4.contour(*lbm_args,**lbm_kwargs)
    ax4.quiver(wind_lon[::q], wind_lat[:-20:q], sx_data[3,:-20:q,::q],sy_data[3,:-20:q,::q], **qui_kwargs)
    ax4.set_title(times[3])
    
    ### create and plot msp axes
    p_dev = 1020
    ax5 = fig.add_subplot(gs[1,0], **splt_kwargs)
    ax5.pcolormesh(prec_lon, prec_lat, prec_data[0]*1000, **prec_kwargs)
    ax5.contour(msp_lon, msp_lat, msp_data[0]/100 -p_dev, **msp_kwargs)
    ax5.contour(*lbm_args, **lbm_kwargs)
    
    ax6 = fig.add_subplot(gs[1,1], **splt_kwargs)
    ax6.pcolormesh(prec_lon, prec_lat, prec_data[1]*1000, **prec_kwargs)
    ax6.contour(msp_lon, msp_lat, msp_data[1]/100 -p_dev, **msp_kwargs)
    ax6.contour(*lbm_args, **lbm_kwargs)
    
    ax7 = fig.add_subplot(gs[1,2], **splt_kwargs)
    ax7.pcolormesh(prec_lon, prec_lat, prec_data[2]*1000, **prec_kwargs)
    ax7.contour(msp_lon, msp_lat, msp_data[2]/100 - p_dev, **msp_kwargs)
    ax7.contour(*lbm_args, **lbm_kwargs)
    
    ax8 = fig.add_subplot(gs[1,3], **splt_kwargs)
    ax8.pcolormesh(prec_lon, prec_lat, prec_data[3]*1000, **prec_kwargs)
    ax8.contour(msp_lon, msp_lat, msp_data[3]/100 -p_dev, **msp_kwargs)
    ax8.contour(*lbm_args, **lbm_kwargs)
    
    ### air temperature axes
    ax9 = fig.add_subplot(gs[2,0], **splt_kwargs)
    ax9.pcolormesh(temp_lon, temp_lat, temp_data[0], **temp_kwargs)
    ax9.contour(*lbm_args, **lbm_kwargs)
    
    ax10 = fig.add_subplot(gs[2,1], **splt_kwargs)
    ax10.pcolormesh(temp_lon, temp_lat, temp_data[1], **temp_kwargs)
    ax10.contour(*lbm_args, **lbm_kwargs)
    
    ax11 = fig.add_subplot(gs[2,2], **splt_kwargs)
    ax11.pcolormesh(temp_lon, temp_lat, temp_data[2], **temp_kwargs)
    ax11.contour(*lbm_args, **lbm_kwargs)
    
    ax12 = fig.add_subplot(gs[2,3], **splt_kwargs)
    ax12.pcolormesh(temp_lon, temp_lat, temp_data[3], **temp_kwargs)
    ax12.contour(*lbm_args, **lbm_kwargs)
    
    ### precipitation flux plots
# =============================================================================
#     ax13 = fig.add_subplot(gs[3,0], **splt_kwargs)
#     ax13.pcolormesh(prec_lon, prec_lat, prec_data[0]*1000, **prec_kwargs)
#     ax13.contour(*lbm_args, **lbm_kwargs)
#     
#     ax14 = fig.add_subplot(gs[3,1], **splt_kwargs)
#     ax14.pcolormesh(prec_lon, prec_lat, prec_data[1]*1000, **prec_kwargs)
#     ax14.contour(*lbm_args, **lbm_kwargs)
#     
#     ax15 = fig.add_subplot(gs[3,2], **splt_kwargs)
#     ax15.pcolormesh(prec_lon, prec_lat, prec_data[2]*1000, **prec_kwargs)
#     ax15.contour(*lbm_args, **lbm_kwargs)
#     
#     ax16 = fig.add_subplot(gs[3,3], **splt_kwargs)
#     ax16.pcolormesh(prec_lon, prec_lat, prec_data[3]*1000, **prec_kwargs)
#     ax16.contour(*lbm_args, **lbm_kwargs)
# =============================================================================
    
    
    ### colorbars
    cbar1 = plt.colorbar(wsp_map,ax = [ax1,ax2,ax3,ax4],  orientation = 'vertical', fraction = 0.012, pad = 0.1)
    cbar1.ax.set(ylabel = r'$ms^{-1}$')
    #cbar2 = plt.colorbar(msp_map,ax = [ax5,ax6,ax7,ax8],  orientation = 'vertical', fraction = 0.012, pad = 0.1)
    #cbar2.ax.set(ylabel = r'$hPa$')
    cbar3 = plt.colorbar(tem_map,ax = [ax9,ax10,ax11,ax12],  orientation = 'vertical', fraction = 0.012, pad = 0.1)
    cbar3.ax.set(ylabel = r'$K$')
    cbar4 = plt.colorbar(prf_map,ax =[ax5,ax6,ax7,ax8] ,  orientation = 'vertical', fraction = 0.012, pad = 0.1)
    cbar4.ax.set(ylabel = r'$gm^{-2}s^{-1}$')
    title = 'UM GA6.1 n768 - Case 2'
    fig.suptitle(title, y = 0.95)
    #cax,kw = matplotlib.colorbar.make_axes([ax for ax in (ax1,ax2,ax3,ax4,ax5,ax6)], orientation = 'horizontal',)
    #plt.colorbar(wsp_map, cax = cax, **kw)
    plt.subplots_adjust(right = 0.89,wspace=0.01, hspace=0.1)
    #plt.tight_layout()
    plt.savefig(f'D:/Project/Figures/PDF/{flight}/{suite}/glm/GLM_um_surf_wsp_msp_tem_prf_3rowcompr{suite[2:]}_flt{flight}.pdf')
    plt.savefig(f'D:/Project/Figures/PNG/{flight}/{suite}/glm/GLM_um_surf_wsp_msp_tem_prf_3rowcompr{suite[2:]}_flt{flight}.png')
    plt.show()


#%%
msp_plot = False
if msp_plot:
    msp_cube = glm_load_and_concat('m01s16i222', kconstraints)
    lbm_cube = glm_load_and_concat('m01s00i030', kconstraints, fileno = 1)
    
    lbm_data, lbm_lon, lbm_lat = extract_data_and_coords(lbm_cube)
    msp_data, msp_lon, msp_lat = extract_data_and_coords(msp_cube)
    
    #%%  time frmatting
    
    cube_time = wind_cube.coord('time')
    dates = cube_time.units.num2date(cube_time.points)
    timepoints = []
    for t in range(len(dates)):
        timepoint = dates[t].strftime('%d-%m-18 %H:%M')
        timepoints.append(timepoint)
    times = np.array(timepoints)
    print(times)
    print(cube_time.points)
    
    #%% find min and max across day for normalisation
    vmin = 0
    vmax = np.nanpercentile(wind_data,99.9)
    print(vmax)
    wsp_norm = matplotlib.colors.Normalize(vmin = vmin, vmax = vmax)
    wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = 'Oranges')
    #%% plot figure
    data_crs = ccrs.PlateCarree()
    gs = gridspec.GridSpec(nrows = 2, ncols = 2)
    projtype = ccrs.NorthPolarStereo()
    splt_kwargs = {'projection' : projtype, 'frameon' : False} # subplot keyword arguments
    wind_kwargs = {'transform':data_crs,'cmap':'Oranges','vmin':vmin,'vmax':vmax}
    lbm_args = [lbm_lon,lbm_lat,lbm_data] # land mask arguments
    lbm_kwargs = {'colors':'k','transform':data_crs} # land mask keyword arguments
    qui_kwargs = {'transform':data_crs} # quiver keywords
    q = 15
    fig = plt.figure(figsize = (16,16))
    
    
    #ax.set_global
    ax1 = fig.add_subplot(gs[0,0],**splt_kwargs)
    #iplt.pcolormesh(wind_cube[0], cmap = 'Oranges')
    ax1.pcolormesh(wind_lon,wind_lat,wind_data[0],**wind_kwargs)
    ax1.contour(*lbm_args, **lbm_kwargs)
    ax1.set_title(times[0])
    
    ax2 = fig.add_subplot(gs[0,1], **splt_kwargs)
    ax2.pcolormesh(wind_lon, wind_lat, wind_data[1], **wind_kwargs)
    ax2.contour(*lbm_args,**lbm_kwargs)
    ax2.set_title(times[1])
    
    ax3 = fig.add_subplot(gs[1,0], **splt_kwargs)
    ax3.pcolormesh(wind_lon, wind_lat, wind_data[2], **wind_kwargs)
    ax3.contour(*lbm_args,**lbm_kwargs)
    ax3.set_title(times[2])
    
    ax4 = fig.add_subplot(gs[1,1], **splt_kwargs)
    ax4.pcolormesh(wind_lon, wind_lat, wind_data[3], **wind_kwargs)
    ax4.contour(*lbm_args,**lbm_kwargs)
    ax4.set_title(times[3])
    
    
    
    #ax7 = fig.add_subplot(gs[2,0])
    cbar = plt.colorbar(wsp_map,ax = [ax1,ax2,ax3,ax4],  orientation = 'horizontal')
    cbar.ax.set(xlabel = r'$ms^{-1}$')
    title = '10m Surface Windspeed and Direction - UM GLM - Control'
    fig.suptitle(title, y = 0.96)
    #cax,kw = matplotlib.colorbar.make_axes([ax for ax in (ax1,ax2,ax3,ax4,ax5,ax6)], orientation = 'horizontal',)
    #plt.colorbar(wsp_map, cax = cax, **kw)
    plt.subplots_adjust(bottom = 0.25,wspace=0.01, hspace=0.1)
    
    plt.savefig(f'D:/Project/Figures/PDF/{flight}/{suite}/glm/GLM_um_10m_windspeed_{suite[2:]}_flt{flight}.pdf')
    plt.savefig(f'D:/Project/Figures/PNG/{flight}/{suite}/glm/GLM_um_10m_windspeed_{suite[2:]}_flt{flight}.png')
    plt.show()