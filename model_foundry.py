# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 14:53:49 2020

@author: Wilhelm Hodder

Foundry contains functions to produce complete figures with various presets.
"""
import numpy as np
import iris
from iris.analysis import trajectory
import iris.quickplot as qplt
import iris.plot as iplt
import iris.coord_categorisation
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import mlab,colors,cm
from matplotlib.font_manager import FontProperties
import matplotlib.patheffects as PathEffects
from scipy import interpolate,stats
#import scipy
import os
import sys
from copy import deepcopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from math import sin, cos, sqrt, atan2, log10
import datetime

import model_functions as modf
import model_foundry_toolkit as toolkit

def make_horizontal_model_contourf_symdiv(data, land_mask, lon, lat, lons, lats,time, level, data_label = '', res = '1p5km', domain = '',
                                   flight = 0,
                                   figure_path_pdf = '', figure_path_png = '', variable_in_file_name = '', savefig = True,
                                   variable_name = ''):
    """
    DEPRECATED.
    Produces horizontal 2D map at specified time and level, with symmetric diverging colourmap.
    """
    # Turn interactive plotting off
    plt.ioff()
    
    ## get min and max of vertical velocity data at eacht time and level
    day_wmin, day_wmax = modf.data_extrema(data[:][level], print_max_min = False) # over day
    
    
    # make list of contour levels
    
    #w_min = np.nanmin(data[i][j])
    levels = np.linspace(day_wmin - 0.01, day_wmax + 0.01)# 10)
    # formatting time string/datetime object
    import datetime
    hour = time//2
    minute = (time/2 - hour)*60
    minute = int(minute); hour = int(hour)
    TOD = datetime.time(hour, minute)  # time of day
    ## start figure
    fig, ax = plt.subplots(1,1, figsize = (18,15))
    q1 = ax.contourf(lon, lat, data[time][level], levels = levels,
                     vmin = day_wmin, vmax = day_wmax, cmap = plt.cm.seismic)
    # add coastlines from model data
    q2 = ax.contour(lon, lat, land_mask)
    # add real latitude and longitude lines
    q4 = ax.contour(lon, lat, lons, colors = 'k')
    q5 = ax.contour(lon, lat, lats, colors = 'k')
    # label latitude and longitude lines
    label_keys = {'fontsize' : 9, 'inline' : 1}
    matplotlib.rcParams['contour.negative_linestyle'] = 'solid' # make contours for negative values solid, not dashed
    q6 = ax.clabel(q4, **label_keys)
    q7 = ax.clabel(q5, **label_keys)
    # hide axis labels
    q8 = ax.set_yticklabels([]); q9 = ax.set_xticklabels([])
    # titles and text labels
    q3 = ax.set(Title = f'{variable_name} at model level {level + 1} at {TOD}UTC',
                ylabel = 'Latitude',
                xlabel = 'Longitude')
    # colourbar
    cbar0 = fig.colorbar(q1, ax = ax)
    cbar0.ax.set(ylabel = f'{data_label}')
    # save figures
    if savefig == True:
        try:
            if flight == 301:
                file_name = f'{res}_{domain}_{variable_in_file_name}_20180312-301_TH{hour}M{minute}_ML{level+1}'
            elif flight == 306:
                file_name = f'{res}_{domain}_{variable_in_file_name}_20180319-306_TH{hour}M{minute}_ML{level+1}'
            else:
                print('Flight not specified, figure not saved.')
                return fig
            fig.savefig(f'{figure_path_pdf}/{domain}/{level+1}/{file_name}.pdf')
            fig.savefig(f'{figure_path_png}/{domain}/{level+1}/{file_name}.png')   
        except:
            print('Saving figure unsuccessful; please check input arguments.')
            return fig
    
    return fig, ax



def plot_hoz_xsec_mass(data, land_mask, lon, lat, land_lon, land_lat, time, level, unrotated = True, cube = None,
                                          data_label = '', res = '1p5km', domain = '', contour_number = 30,
                                          flight = 0, time_norm = True,
                                          figure_path_pdf = None, figure_path_png = None, variable_in_file_name = '', savefig = True,
                                          variable_name = '', colourmap = 'plasma', quiver = 'False', quiver_u = [], quiver_v = [],
                                          quiver_n = 1, title = None,
                                          contour_type = 'contourf', obvs_on = False, msp_contour = False, msp_data = None, msp_lat = None, msp_lon = None,
                                          msp_levels = [], msp_style = 'solid', gph_contour = False, gph_data = None, gph_lon = None, gph_lat = None, 
                                          gph_style = 'dashed', gph_levels = [],
                                          obvs_z = [], obvs_x = [], obvs_y = [],
                                          figsize = (18,15), obvs_size = 150.0,
                                          coast = True, vmin = None, vmax = None, make_animation = False,
                                          sampling_rate = 'halfhourly', sym_cmap = False,
                                          surface_variable = False, pressure_levels = False,
                                          cscale_override = False, cscale_min = None, cscale_max = None, verstash = False,
                                          pressure_labels = ['200hPa','300hPa', '500hPa','800hPa','850hPa','950hPa'], PDF = True, PNG = True):
    """
    Produces horizontal cross-sections from various data arrays, not cubes. Also more suitable for mass-production of cross-sections.
    """
    
    # Turn interactive plotting off
    plt.ioff()
    # Global settings
    font = {'size' : 18}
    matplotlib.rc('font', **font)
    
    ## get min and max of vertical velocity data at each time and level
    if cscale_override == False:
        if time_norm == True:
            if obvs_on == True:
                model_min = np.nanmin(data[:,level]); obvs_min = np.nanmin(obvs_z)
                model_max = np.nanmax(data[:,level]); obvs_max = np.nanmax(obvs_z)
                day_wmin = np.nanmin(np.array([model_min, obvs_min]))
                day_wmax = np.nanmax(np.array([model_max, obvs_max]))
            else:   
                day_wmin = np.nanmin(data[:,level])
                day_wmax = np.nanmax(data[:,level])
            if vmin != None:
                day_wmin = vmin
            if vmax != None:
                day_wmax = vmax
            if sym_cmap == False:
                my_norm= matplotlib.colors.Normalize(vmin = day_wmin,vmax= day_wmax, clip=True)
            elif sym_cmap == True:
                # put min and max into list and calculate absolutes elementwise with numpy, then take maximum value
                sym_max = np.max(np.absolute([day_wmin, day_wmax])) 
                sym_min = -sym_max
                my_norm= matplotlib.colors.Normalize(vmin = sym_min,vmax= sym_max, clip=True)
        if time_norm == False: # if colourmap is to be normalised over each instance only
            if obvs_on == True:
                model_min = np.nanmin(data[time][level]); obvs_min = np.nanmin(obvs_z)
                model_max = np.nanmax(data[time][level]); obvs_max = np.nanmax(obvs_z)
                day_wmin = np.nanmin(np.array([model_min, obvs_min]))
                day_wmax = np.nanmax(np.array([model_max, obvs_max]))
            else:   
                day_wmin = np.nanmin(data[time][level])
                day_wmax = np.nanmax(data[time][level])
            if vmin != None:
                day_wmin = vmin
            if vmax != None:
                day_wmax = vmax
            if sym_cmap == False:
                my_norm = matplotlib.colors.Normalize(vmin = day_wmin,vmax= day_wmax, clip=True)
            elif sym_cmap == True:
                # put min and max into list and calculate absolutes elementwise with numpy, then take maximum value
                sym_max = np.max(np.absolute([day_wmin, day_wmax])) 
                sym_min = -sym_max
                my_norm = matplotlib.colors.Normalize(vmin = sym_min,vmax= sym_max, clip=True)
       
    elif cscale_override == True: # use this to set manual clour scale limits
        day_wmin = cscale_min; day_wmax = cscale_max
        my_norm = matplotlib.colors.Normalize(vmin = day_wmin,vmax= day_wmax, clip=True)
        
    
    
    # formatting time string/datetime object
    if sampling_rate == 'halfhourly':
        hour = time//2
        minute = (time/2 - hour)*60
    elif sampling_rate == 'hourly':
        hour = time//1
        minute = (time-hour)*60
    elif sampling_rate == 'threehourly':
        hour = (time//1)*3
        minute = 0
    minute = int(minute); hour = int(hour)
    TOD = datetime.time(hour, minute)  # time of day
    ## start figure
    fig, ax = plt.subplots(1,1, figsize = figsize)
    if surface_variable == False:
        if contour_type == 'contourf':
            # make list of contour levels
            contour_levels = np.linspace(day_wmin, day_wmax, contour_number)
            q1 = ax.contourf(lon, lat, data[time][level], norm = my_norm, cmap = colourmap, levels = contour_levels)
        elif contour_type == 'contour':
            # make list of contour levels
            contour_levels = np.linspace(day_wmin, day_wmax, contour_number)
            q1 = ax.contour(lon, lat, data[time][level], norm = my_norm, cmap = colourmap, levels = contour_levels)
        elif contour_type == 'pcolor':
            q1 = ax.pcolormesh(lon, lat, data[time][level], norm = my_norm, cmap = colourmap)
        else:
            print('Contour type not recognised.')
            return None
    elif surface_variable == True:
        if contour_type == 'contourf':
            q1 = ax.contourf(lon, lat, data[time], norm = my_norm, cmap = colourmap, levels = contour_levels)
        elif contour_type == 'contour':
            q1 = ax.contour(lon, lat, data[time], norm = my_norm, cmap = colourmap, levels = contour_levels)
        elif contour_type == 'pcolor':
            q1 = ax.pcolormesh(lon, lat, data[time], norm = my_norm, cmap = colourmap)
        else:
            print('Contour type not recognised.')
            return None
    # add coastlines from model data
    if coast == True:
        q2 = ax.contour(land_lon, land_lat, land_mask, colors = 'k')
    # add real latitude and longitude lines
    if unrotated == True:
        gridlon = cube.coords('grid_longitude')[0].points
        gridlat = cube.coords('grid_latitude')[0].points
        lons, lats = foundry.modf.unrotate_coords(gridlon, gridlat, polelon, polelat)
        q4 = ax.contour(lon, lat, lons, colors = 'k', linestyles = 'dotted')
        q5 = ax.contour(lon, lat, lats, colors = 'k', linestyles = 'dotted')
        # label latitude and longitude lines
        label_keys = {'fontsize' : 10, 'inline' : 1}
        matplotlib.rcParams['contour.negative_linestyle'] = 'dotted' # make contours for negative values solid, not dashed
        q6 = ax.clabel(q4, **label_keys)
        q7 = ax.clabel(q5, **label_keys)
    # optional quiverplot
    if quiver == True:
        n = quiver_n
        if surface_variable == True:
            u_n = quiver_u[time][::n]
            v_n = quiver_v[time][::n]
        elif surface_variable == False:
            u_n = quiver_u[time][level][::n]
            v_n = quiver_v[time][level][::n]
        
        u_n2 = []; v_n2 = []
        for j in range(len(u_n)):
            u_n2.append(u_n[j][0::n])
            v_n2.append(v_n[j][0::n])
        print(len(lon[0::n]),len(lat[0::n]),np.shape(u_n2),np.shape(v_n2))
        qq = ax.quiver(lon[0::n], lat[0::n], u_n2, v_n2, pivot = 'middle')
    # hide axis labels
    q8 = ax.set_yticklabels([]); q9 = ax.set_xticklabels([])
    ## second contour: e.g geopotential height or mean-sealevel pressure
    if msp_contour == True:
        second_ax = ax.contour(msp_lon, msp_lat, msp_data, colors = 'k', linestyles = msp_style, levels = msp_levels)
        # Add contour labels based on the contour we have just created.
        plt.clabel(second_ax, inline=False, fontsize = 14)
    ## third contour: e.g geopotential height or mean-sealevel pressure
    if gph_contour == True:
        third_ax = ax.contour(gph_lon, gph_lat, gph_data[time][level], colors = 'k', linestyles = gph_style, levels = gph_levels)
        # Add contour labels based on the contour we have just created.
        plt.clabel(third_ax, inline=False, fontsize = 14)
    # optional observational data (scatter plot)
    if obvs_on == True:
        q9 = plt.scatter(x = obvs_x, y = obvs_y, c = obvs_z, cmap = colourmap, vmin = day_wmin, vmax = day_wmax, s = obvs_size )
    # titles and text labels
    if title == None:
        if surface_variable == False:
            if pressure_levels == False:
                q3 = ax.set(Title = f'{variable_name} at {res} resolution at model level {level + 1} and {TOD}UTC')#,
            elif pressure_levels == True:
                # special titeling if we are iterating over pressure levels
                # levels specific to my choice in suite settings
                q3 = ax.set(Title = f'{variable_name} at {res} resolution at {pressure_labels[level]} and {TOD}UTC')
        elif surface_variable == True:
            q3 = ax.set(Title = f'{variable_name} at {res} resolution at {TOD}UTC')
    else:
        q3 = ax.set(Title = title)
    # colourbar
    cbar0 = fig.colorbar(q1, ax = ax)
    cbar0.ax.set(ylabel = f'{data_label}')
    # save figures
    if savefig == True:
        toolkit.save_xsection(fig, res, domain, variable_in_file_name, flight, level = level + 1, verstash = verstash, PDF = PDF, PNG = PNG,
                      figure_path_pdf = figure_path_pdf, figure_path_png = figure_path_png, 
                      time = time, sampling_rate = sampling_rate, surface_variable = surface_variable, pressure_levels = pressure_levels)
    
    return fig, ax, my_norm

def plot_vert_xsec_on_model_axis(data_cube, og_cube, time, wide_coor_start, wide_coor_end, point_coor, axis, bottom = 0, top = -1,
                              figsize = (9,4), norm_offset = 0.5, cmap = plt.cm.viridis, plot_type = 'contourf', contour_levels = 21, sym_cmap = False,
                              isentropes_on = False, isentropes_step = 1, theta_cube = None, isentr_cont = 21, isentr_labels = False, cbar_orientation = 'horizontal',
                              var_name_in_title = None, line_on_map = False, land_mask_cube = None,
                              ylim_top = 6000, ylim_bottom = 0, divisor = 50):
    
    """
    Plot vertical slice through a given cube.
    colour scale normalisation done cross-temporally (equal colour scale across all time steps.)
    
    Positional args:
    cube - your iris cube of the desired variable
    time - specifically the time index of your desired moment of time.
    wide_coor_start - the first index point of the wide section of your grid that you want to plot
    wide_coor_end - the last index point of the above
    point_coor - the specific index point you want to extract data from opn the other coordinate axis.
    axis - specify which axis you are plotting along (longitudinal or latitudinal) Note: this does not take into account the real Earth grid, and is hence only advisory.
    
    Keyword args:
    bottom - index of bottom level: default 0
    top - index of top level: default -1 (the highest value index)
    figsize - tuple of size of final figure in inches
    norm_offset - offset from min/max values in coloursale normalisation
    cmap - colour map
    plot_type - 'pcolor', contourf', or 'contour' 
    contour_levels - integer number of contour levels or list of specified levels to be passed to contour/contourf functions
    sym_cmap - toggle for normalising colourmap symmetrically about zero or not
    """
    
    ## translate coordinates into index 
    # not implemented in this function, as this one chooses the exact axis
    wide_index_start = wide_coor_start; wide_index_end = wide_coor_end; point_index = point_coor
    ## make new cube from spaital slice of main cube, but containing all time steps - needed for cross-temporal normalisation
    if axis == 'longitudinal':
        cube_slice_wta = data_cube[:,bottom:top,wide_index_start:wide_index_end,point_index] # wt suffix refers to 'with time axis'
        #orography and theta (which will be plotted in black contours) don't need to retain time dimension
        og_cube_slice = og_cube[0,0,wide_index_start:wide_index_end,point_index]
        if isentropes_on == True: # use theta and slice it 
            th_cube_slice = theta_cube[time,bottom:top,wide_index_start:wide_index_end,point_index]
        
        
    if axis == 'latitudinal':
        cube_slice_wta = data_cube[:,bottom:top,point_index, wide_index_start:wide_index_end]
        #orography and theta (which will be plotted in black contours) don't need to retain time dimension
        og_cube_slice = og_cube[0,0,point_index, wide_index_start:wide_index_end]
        if isentropes_on == True: # use theta and slice it 
            th_cube_slice = theta_cube[time,bottom:top,point_index, wide_index_start:wide_index_end]
        
    
    ## distance calculation for x-axis labelling
    select_coords, select_dist = toolkit.singleaxis_distance_xticklists(cube_slice_wta, axis, point_index, divisor)    
    print(' select coords =\n',select_coords)
    ## set colourmap and normalise colourmap to data
    my_cmap = cmap
    my_cmap.set_under('w')
    ## normalise the colourmap to be between a certain min and max
    data_min = np.nanmin(cube_slice_wta.data[:]) - norm_offset;
    data_max = np.nanmax(cube_slice_wta.data[:]) + norm_offset
    if sym_cmap == False:
        my_norm= matplotlib.colors.Normalize(vmin = data_min,vmax= data_max, clip=True)
    elif sym_cmap == True:
        # put min and max into list and calculate absolutes elementwise with numpy, then take maximum value
        sym_max = np.max(np.absolute([data_min, data_max])) 
        sym_min = -sym_max
        my_norm= matplotlib.colors.Normalize(vmin = sym_min,vmax= sym_max, clip=True)
    
    # generate final version of w cube
    cube_slice = cube_slice_wta[time] # single out point in time
    
    fig0=plt.figure(figsize=figsize)
    if plot_type == 'pcolor':
        ax = iplt.pcolormesh(cube_slice, norm = my_norm, cmap = my_cmap)
    elif plot_type == 'contourf':
        ax = iplt.contourf(cube_slice, norm = my_norm, cmap = my_cmap, levels = contour_levels)
    elif plot_type == 'contour':
        ax = iplt.contour(cube_slice, norm = my_norm, cmap = my_cmap, levels = contour_levels)
    cbar = plt.colorbar(ax, orientation=cbar_orientation); cbar.ax.set_xlabel('Windspeed [ms$^{-1}$]')
    if isentropes_on == True: # plot isentropes
        thcontlevs=np.linspace(250,1000,750//isentropes_step + 1) # sample theta at every full degree Kelvin
        isenax = iplt.contour(th_cube_slice, levels = thcontlevs, colors = 'k', linewidths = 0.75)
        if isentr_labels == True:
            # Add contour labels based on the contour we have just created.
            plt.clabel(isenax, inline=False, fontsize = 14)
    # plot orogrpahy
    iplt.plot(og_cube_slice, color ='k', linewidth = 1.5)
    # set y limit - cut off over pre-ploted area for clean cut
    plt.ylim(top = ylim_top, bottom = ylim_bottom)
    # trim front of plot as well
    plt.xlim(left = select_coords[0], right = select_coords[-1])
    # override ticks and tick labes on x-axis
    plt.xticks(select_coords, select_dist) # at positions of select_coords, set labels as select_dist
    # set axis labels and title
    plt.ylabel('Model altitude [m]')
    plt.xlabel('Distance scale [km]')
    plt.title(f'{var_name_in_title} cross-section on {axis} axis')
    
    plt.tight_layout()
    
    ## plot line of cross-section onto horizontal map
    if line_on_map == True:
        ## take bottom line of data cube
        latitude = cube_slice.coord('grid_latitude').points
        longitude = cube_slice.coord('grid_longitude').points
        latitude = np.array([np.min(latitude), np.max(longitude)])
        longitude = np.array([np.min(longitude),np.max(longitude)])
        for i in [0,1]:
            if latitude[i] > 180:
                latitude[i] -= 360
            if longitude[i] > 180:
                longitude[i] -= 360
        print(latitude, longitude)
        landmask = land_mask_cube.data
        land_lon = np.array(land_mask_cube.coord('grid_longitude').points)
        land_lat = np.array(land_mask_cube.coord('grid_latitude').points)
        for i in range(len(land_lon)):
            if land_lon[i] > 180:
                land_lon[i] -= 360
        for i in range(len(land_lat)):
            if land_lat[i] > 180:
                land_lat[i] -= 360
        print(len(landmask), len(land_lon), len(land_lat))
        fig1 = plt.figure(figsize=figsize)
        #iplt.pcolormesh(bottom_line)
       # qplt.contour(data_cube[0][0])
        plt.contour(land_lon, land_lat, landmask, color = 'k')
        plt.plot(longitude, latitude, color = 'r')
        
        return fig0, fig1
        
        
    return fig0



def plot_vert_xsec_on_variant_axis(data_cube, og_cube, time, coor_1, coor_2, bottom = 0, top = -1,
                              figsize = (9,4), norm_offset = 0.5, cmap = plt.cm.viridis, plot_type = 'contourf', contour_levels = 21, sym_cmap = False,
                              isentropes_on = False, theta_cube = None, isentr_cont = 21, isentr_labels = False, cbar_orientation = 'horizontal',
                              var_name_in_title = None, divisor = 50,ylim_top = 6000, ylim_bottom = 0, points_in_line = 100,
                              obs_on = False, obs_lon = None, obs_lat = None, obs_z = None, obs_values = None,
                              obs_leg_startend_file = None, res = ''):
    
    """
    Plot vertical slice through a given cube.
    colour scale normalisation done cross-temporally (equal colour scale across all time steps.)
    Interpolates between two given coordinates
    
    Positional args:
    cube - your iris cube of the desired variable
    time - specifically the time index of your desired moment of time.
    coor_1 and coor_2 - the coordinates between whcih to interpolate the x-section, given as (lon,lat) tuples
    axis - specify which axis you are plotting along (string: 'longitudinal' or 'latitudinal') Note: this does not take into account the real Earth grid, and is hence only advisory.
    
    Keyword args:
    bottom - index of bottom level: default 0
    top - index of top level: default -1 (the highest value index)
    figsize - tuple of size of final figure in inches
    norm_offset - offset from min/max values in coloursale normalisation
    cmap - colour map
    plot_type - 'pcolor', contourf', or 'contour' 
    contour_levels - integer number of contour levels or list of specified levels to be passed to contour/contourf functions
    sym_cmap - toggle for normalising colourmap symmetrically about zero or not
    """
    ### Rotate coordinates to model frame
    # extract rotated north pole coordinates
    NPoleLon = data_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    NPoleLat = data_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    # rotate database coordinates
    lons = np.array([coor_1[0], coor_1[1]])
    lats = np.array([coor_2[0], coor_2[1]])
    rot_lons, rot_lats = iris.analysis.cartography.rotate_pole(lons, lats, NPoleLon, NPoleLat)
    
    ### We need to generate a line of sample points between two specified (lat,lon) coordinates in rotated frame
    ### and interpolate our data and orog cube given those sample points
    if points_in_line == 'auto':
        span = toolkit.haversine(rot_lats[0], rot_lons[0], rot_lats[1], rot_lons[1])
        if res == '0p5km': resn = 0.5
        if res == '1p5km': resn = 1.6
        if res == '4p4km': resn = 4.4
        cell_number = span / resn
        auto_points = cell_number//1; print('number of points =', auto_points)
        slice_lons = np.linspace(rot_lons[0], rot_lons[1], auto_points)
        slice_lats = np.linspace(rot_lats[0], rot_lats[1], auto_points)
    else:
        slice_lons = np.linspace(rot_lons[0], rot_lons[1], points_in_line)
        slice_lats = np.linspace(rot_lats[0], rot_lats[1], points_in_line)
    # interpolate cubes into slices using nearest neighbour method
    cube_slice = trajectory.interpolate(data_cube, [('grid_longitude', slice_lons), ('grid_latitude', slice_lats)], method='nearest')
    og_cube_slice = trajectory.interpolate(og_cube, [('grid_longitude', slice_lons), ('grid_latitude', slice_lats)], method='nearest')
    print('cube_slice:\n',cube_slice)    
    
    ## set colourmap and normalise colourmap to data
    my_cmap = cmap
    my_cmap.set_under('w')
    ## normalise the colourmap to be between a certain min and max
    data_min = np.nanmin(cube_slice[time].data[:]) - norm_offset;
    data_max = np.nanmax(cube_slice[time].data[:]) + norm_offset
    if sym_cmap == False:
        my_norm= matplotlib.colors.Normalize(vmin = data_min,vmax= data_max, clip=True)
    elif sym_cmap == True:
        # put min and max into list and calculate absolutes elementwise with numpy, then take maximum value
        sym_max = np.max(np.absolute([data_min, data_max])) 
        sym_min = -sym_max
        my_norm= matplotlib.colors.Normalize(vmin = sym_min,vmax= sym_max, clip=True)
    
    ## distance calculation for x-axis labelling
    select_lons, select_lats, select_dist, select_index = toolkit.variableaxis_distance_xticklists(slice_lons, slice_lats, divisor = divisor)
    sample_indices = np.arange(1,auto_points + 1)
    sample_indices = sample_indices[select_index]
    print(select_index, sample_indices)
    ## format observational data coordinates
    if obs_on == True:
        ptpdist = toolkit.haversine(obs_lat[:-1],obs_lon[:-1],obs_lat[1:],obs_lon[1:])
        obsdist = np.cumsum(ptpdist) # elementwise cumulative sum
        obsdist = np.concatenate(([0], obsdist)) # add zero to beginning of array; now same length as coords
        
    ## now start plotting and formatting the figure
    fig=plt.figure(figsize=figsize)
    if plot_type == 'pcolor':
        ax = iplt.pcolormesh(cube_slice[time], norm = my_norm, cmap = my_cmap)
    elif plot_type == 'contourf':
        ax = iplt.contourf(cube_slice[time], norm = my_norm, cmap = my_cmap, levels = contour_levels)
    elif plot_type == 'contour':
        ax = iplt.contour(cube_slice[time], norm = my_norm, cmap = my_cmap, levels = contour_levels)
    cbar = plt.colorbar(ax, orientation=cbar_orientation); cbar.ax.set_xlabel('Windspeed [ms$^{-1}$]')
    if isentropes_on == True: # plot isentropes
        th_cube_slice = trajectory.interpolate(theta_cube, [('grid_longitude', slice_lons), ('grid_latitude', slice_lats)], method='nearest')
        thcontlevs=np.linspace(250,1000,751) # sample theta at every full degree Kelvin
        isenax = iplt.contour(th_cube_slice[time], levels = thcontlevs, colors = 'k', linewidths = 0.75)
        if isentr_labels == True:
            # Add contour labels based on the contour we have just created.
            plt.clabel(isenax, inline=False, fontsize = 14)
    # plot orogrpahy
    iplt.plot(og_cube_slice, color ='k', linewidth = 1.5)
    # set y limit - cut off over pre-ploted area for clean cut
    plt.ylim(top = ylim_top, bottom = ylim_bottom)
    # trim front of plot as well
    #plt.xlim(left = select_coords[0], right = select_coords[-1])
    # override ticks and tick labes on x-axis
    plt.xticks(select_index, select_dist) # at positions of select_coords, set labels as select_dist
    # set axis labels and title
    plt.ylabel('Model altitude [m]')
    plt.xlabel('Distance scale [km]')
    plt.title(f'{var_name_in_title} cross-section between\n({coor_1[1]},{coor_1[0]}) and ({coor_2[1]},{coor_2[0]})\n at {res} resolution at ')
    
    #plt.tight_layout()
    
    return fig

def plot_vert_xsec_single(wvelppname,uvelppname, vvelppname,thppname,oppname,slicetype,slicedttm,lonstart,latstart,lonend,latend,ylims,
                   velcontmax,prefix = None,outdir = None,loccoords=None,resolvedir=None, savefig = False, res = None, 
                   domain = None, variable_in_file_name = None, time = None, 
                   figure_path_pdf = None, figure_path_png = None, flight = None, sampling_rate = 'halfhourly', show = False):
    '''
    ADAPTED from Peter Sheridan (Met Office) - CubeCrossSectioner.py
    Most suitable for single/small numbers of plots
    
    Plot specified vertical slice through 3D model data
    Vertical or hoz velocities, or potentially other field, 
    plotted with relative distance coords and terrain outline  
    (hoz vectors using quiver not yet implemented)
    
    "plot_vert_xsec(velppname,thppname,oppname,slicetype,slicedttm,lonstart,latstart,lonend,latend,ylims,\
                    velcontmax,prefix,outdir,loccoords=None,resolvedir=None)"
    velppname,thppname,oppname: name of data files for velocities, theta and orography
    slicetype: vertvel (vertical velocity), hozcompt (horizontal wind component in-plane or 
                resolved along direction given by resolvedir) or hozspeed (horizontal wind speed)
    resolvedir: direction into which winds are resolved for plotting (default in-plane), 
                given as a vector orientation (as opposed to a normal angular wind direction definition) 
    slicedttm: list containing [year,month,day,hourminutes]
    lonstart,latstart,lonend,latend: start and end coordinates of vertical slice in unrotated lon/lat
    ylims: limits of y axis
    velcontmax: colour contour maximum
    prefix: prefix to add to image filenames
    outdir: directory to output image files
    loccoords: coordinates of a location to plot (projected) in the depicted plane
    '''
    
    npointsinspan=100 #no. points along the section 
    #analdttm=[filepart for filepart in velppname.split('/')[1:] if filepart[-1] == 'Z'][0]# and filepart[-6] == 'T'][0]
    #validdttm=str(slicedttm[0])+''.join(['%02i' % dt for dt in slicedttm[1:]])#''.join(slicedttm[1:])
    #print validdttm
    slicedttm=datetime.datetime(slicedttm[0],slicedttm[1],slicedttm[2],slicedttm[3],slicedttm[4])
    #print 'slicedttm',slicedttm
    
    affix='%06.2f_%06.2f_to_%06.2f_%06.2f_%05i-%05im_contmax%02i' % ( lonstart,latstart,lonend,latend,ylims[0],ylims[1],velcontmax ) # !!!outdated string formatting

    if loccoords: 
        loclatlondict={}
        for loc in loccoords:
            loclatlondict[loc[0]]=[loc[1],loc[2]]

    thcube=iris.load_cube(thppname,iris.AttributeConstraint(STASH='m01s00i004'))
    wcube=iris.load_cube(wvelppname,iris.AttributeConstraint(STASH='m01s00i150')) 
    ucube=iris.load_cube(uvelppname,iris.AttributeConstraint(STASH='m01s00i002')) 
    vcube=iris.load_cube(vvelppname,iris.AttributeConstraint(STASH='m01s00i003')) 
    ocube=iris.load_cube(oppname,iris.AttributeConstraint(STASH='m01s00i033'))
    ocube=ocube.regrid(thcube,iris.analysis.Linear())
    thcube,wcube,ucube,vcube=toolkit.add_hybrid_height(ocube,[thcube,wcube,ucube,vcube])
    
    polelon=ocube.coord_system('RotatedGeogCS').grid_north_pole_longitude
    polelat=ocube.coord_system('RotatedGeogCS').grid_north_pole_latitude
    if ocube.coord_system('RotatedGeogCS') != None and polelon != 0 and polelat != 90: rotated=True

    #print 'slice location',(lonstart,latstart),(lonend,latend)
    
    if rotated: 
        (lonstart,lonend),(latstart,latend)=iris.analysis.cartography.rotate_pole(np.array([lonstart,lonend]),np.array([latstart,latend]),polelon,polelat)
        if loccoords:
            for location in loclatlondict.keys():
                loclatlondict[location][0],loclatlondict[location][1]=iris.analysis.cartography.rotate_pole(
                 np.array([loclatlondict[location][0]]),np.array([loclatlondict[location][1]]),polelon,polelat)
                if loclatlondict[location][0] < 0: loclatlondict[location][0]+=360.

    if lonstart < 0: lonstart+=360.
    if lonend < 0: lonend+=360.
    slicelons = np.linspace(lonstart,lonend, num=npointsinspan)
    slicelats = np.linspace(latstart,latend, num=npointsinspan)
    
    with iris.FUTURE.context(cell_datetime_objects = True):
        timeconstraint=iris.Constraint(time=lambda cell: cell.point == slicedttm)
        thcube=thcube.extract(timeconstraint)
        wcube=wcube.extract(timeconstraint)
        ucube=ucube.extract(timeconstraint)
        vcube=vcube.extract(timeconstraint)

        #thcube=thcube.extract(timeconstraint).collapsed('time',iris.analysis.MEAN)
        #ucube=ucube.extract(timeconstraint).collapsed('time',iris.analysis.MEAN)
        #vcube=vcube.extract(timeconstraint).collapsed('time',iris.analysis.MEAN)
        #wcube=wcube.extract(timeconstraint).collapsed('time',iris.analysis.MEAN)

    wslice = trajectory.interpolate(wcube, (('grid_longitude', slicelons), ('grid_latitude', slicelats)), method='nearest')
    uslice = trajectory.interpolate(ucube, (('grid_longitude', slicelons), ('grid_latitude', slicelats)), method='nearest')
    vslice = trajectory.interpolate(vcube, (('grid_longitude', slicelons), ('grid_latitude', slicelats)), method='nearest')
    thslice = trajectory.interpolate(thcube, (('grid_longitude', slicelons), ('grid_latitude', slicelats)), method='nearest') 
    orogslice = trajectory.interpolate(ocube, (('grid_longitude', slicelons), ('grid_latitude', slicelats)), method='nearest')

    #print "obtained vertical data slice"
    
    # plot figure
    fig=plt.figure(figsize=(9,4))

    affix+='_'+slicetype
    if slicetype == 'hozcompt' and resolvedir: affix+='_%05.1f' % (resolvedir)
    thcontlevs=np.linspace(250,1000,376)

    if slicetype == "vertvel":
        wcontlevs=np.linspace(-velcontmax,velcontmax,4*velcontmax+1)

        plotobj=iplt.contourf(wslice,levels=wcontlevs,cmap=toolkit.wcmap(),extend="both")

    elif slicetype == "hozcompt": 

        comptslice=toolkit.getslicecompt(uslice,vslice,lonstart,latstart,lonend,latend,resolvedir=resolvedir)
        if velcontmax == None:
            iviewedlevs=np.where(comptslice.coord('level_height').points <= ylims[1])[0]
            velcontmax=( int(np.max(comptslice[iviewedlevs,:].data)) / 2 ) * 2
        velcontlevs=np.linspace(-velcontmax/2,velcontmax,1+3*velcontmax/2)

        plotobj=iplt.contourf(comptslice,levels=velcontlevs,cmap=toolkit.velcmap(),extend="both")

    elif slicetype == "hozspeed": 

        spdslice=deepcopy(uslice)
        spdslice.data=(uslice.data**2 + vslice.data**2)**0.5
        if velcontmax == None:
            velcontmax=( int(np.max(spdslice.data)) / 2 ) * 2
        velcontlevs=np.linspace(-velcontmax/2,velcontmax,1+3*velcontmax/2)

        plotobj=iplt.contourf(spdslice,levels=velcontlevs,cmap=toolkit.velcmap(),extend="both")

    else: raise Exception("unsupported slice type")

    # plot isentropes        
    iplt.contour(thslice,levels=thcontlevs,colors='k',linewidths=.75)

    # plot vectors
    extentkm=toolkit.haversine(lonstart,latstart,lonend,latend)
    griddx=extentkm/(npointsinspan-1)
    xvals=np.arange(npointsinspan)
    if not(slicetype == "hozcompt"): comptslice=toolkit.getslicecompt(uslice,vslice,lonstart,latstart,lonend,latend)
    #ilevs=np.where(comptslice.coord('level_height').points <= ylims[1])

    #xgrid=np.tile( xvals , [comptslice.data.shape[0],1] )[ilevs,:]
    #ygridu=comptslice.coord('altitude').points[ilevs,:]
    #ygridw=wslice.coord('altitude').points[ilevs,:]

    #nvect=[5,5]
    #nskipx=max([xgrid.shape[1]/nvect[0],1])
    #print 'nskipx',nskipx
    #print 'xgrid.shape,ygridu.shape',xgrid.shape,ygridu.shape
    #xvect,yvect = np.meshgrid(xvals[::nskipx] , np.linspace(ylims[1]/nvect[1],ylims[1],nvect[1]))
    #uvect=mlab.griddata(xgrid.flatten(),ygridu.flatten(),comptslice.data[ilevs,:].flatten(),xvect,yvect,interp='nn')
    #wvect=mlab.griddata(xgrid.flatten(),ygridw.flatten(),wslice.data[ilevs,:].flatten(),xvect,yvect,interp='nn')

    #ax=plt.gca()
    #spdscale=300.
    #arrowobj=ax.quiver(xvect, yvect, uvect, 1000*wvect,minshaft=0.,\
    #                   angles='xy',width=0.002,scale_units='width',scale=nvect[0]*spdscale,\
    #                   headwidth=10,headlength=14) 
    #plt.quiverkey(arrowobj,0.05,-0.15,20.,'20 m/s',color='k')

    # plot orography
    iplt.plot(orogslice,color='k')
    plt.fill_between(np.arange(npointsinspan),orogslice.data,color='k')#,zorder=zorder)
    plt.ylim(ylims)

    # plot locations
    if loccoords: 
        #print 'adding station locations....'
        ilocs,xlocs=toolkit.find_closest_point(loclatlondict,xvals,orogslice)
        for i,xloc in enumerate(xlocs):
            #print 'orogslice.data[ilocs[i]]',orogslice.data[ilocs[i]]
            plt.plot(xloc,orogslice.data[ilocs[i]]+ylims[1]/80.,'k*',markeredgewidth=1,markeredgecolor='w',markersize=14)

    # plot sonde path



    xtickvals=toolkit.tickspacecalc(extentkm)
    xticklocs=xvals[0]+(xtickvals/griddx)*(xvals[-1]-xvals[0])/(len(xvals)-1)
    plt.xticks(xticklocs,xtickvals)
    plt.xlabel('x (km)')
    plt.ylabel('z (m)')
    
    bar=plt.colorbar(plotobj,orientation='vertical')
    if slicetype == 'vertvel': barlabel='vertical velocity (ms$^{-1})$'
    elif slicetype == 'hozcompt': 
        if not resolvedir: barlabel='in-plane horizontal velocity (ms$^{-1})$'
        else: barlabel='resolved horizontal velocity component (ms$^{-1})$'
    elif slicetype == 'hozspeed': barlabel='horizontal wind speed (ms$^{-1})$'
    bar.set_label(barlabel)
    
    plt.tight_layout()
    if savefig == True:
        
        toolkit.save_xsection(fig, res, domain, variable_in_file_name, flight, 
                      figure_path_pdf = figure_path_pdf, figure_path_png = figure_path_png, 
                      time = time, sampling_rate = sampling_rate)
    if show == False:    
        plt.close()
    elif show == True:
        plt.show()
    else:
        print('Please use boolean for \'show\' keyword. Closing figure.')
        plt.close()

