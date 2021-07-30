#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 17:37:05 2020

@author: kse18nru

contains home written functions to help with model data interpolation
"""

import numpy as np
#import numbers
from iris.analysis import trajectory

def basic_interp(xa, xb, xc, ya, yb):
    if np.isnan(xa):# check if bounds are valid
        return np.nan # return non-number
    else: # otherwise proceed with interpolation
        grad = (yb - ya)/(xb - xa) # gradient of interpolant
        dx = xc - xa # height of target above lower point
        yc = grad*dx + ya # carry out interpolation
        return yc

def find_boundpoints(data, alts, lon_idx = None, lat_idx = None, t_idx= None, xc = None,coloumns = False, c_idx = None):
    
    if coloumns:
        ## if data cube comes pre interpolated into a series of coloumns without multiple time coords
   
        ypillar = data[:,c_idx] # coloumn of diagnostic data points
    
        xpillar = alts[:,c_idx] # coloumn of cube altitudes
        
        dpillar = xpillar - xc # transform alt coloumn to set target altitude to zero
        
        neg_section = dpillar[np.where(dpillar < 0)] # subset negative values (linear, will always be together at array bottom )
        if len(neg_section) == 0: # if no negative values are present
            return np.nan, np.nan, np.nan, np.nan # return non-numbers
        else:
            low_idx = len(neg_section) - 1 # otherwise use subset length to find lower boundary index
            
            xa = xpillar[low_idx] 
            xb = xpillar[low_idx + 1]
            ya = ypillar[low_idx]
            yb = ypillar[low_idx + 1]
            return xa,xb,ya,yb # return boundary points
    else:
        ypillar = data[t_idx,:,lat_idx,lon_idx] # coloumn of diagnostic data points
        
        xpillar = alts[:,lat_idx,lon_idx] # coloumn of cube altitudes
        
        dpillar = xpillar - xc # transform alt coloumn to set target altitude to zero
        
        neg_section = dpillar[np.where(dpillar < 0)] # subset negative values (linear, will always be together at array bottom )
        if len(neg_section) == 0: # if no negative values are present
            return np.nan, np.nan, np.nan, np.nan # return non-numbers
        else:
            low_idx = len(neg_section) - 1 # otherwise use subset length to find lower boundary index
            
            xa = xpillar[low_idx] 
            xb = xpillar[low_idx + 1]
            ya = ypillar[low_idx]
            yb = ypillar[low_idx + 1]
            return xa,xb,ya,yb # return boundary points
    
def find_series_boundpoints(data, alts, coor_idx, t_idx, xc, time_req = True):
    ypillar = data[t_idx,:,coor_idx] # coloumn of diagnostic data points
    
    xpillar = alts[:,coor_idx] # coloumn of cube altitudes
    
    dpillar = xpillar - xc # transform alt coloumn to set target altitude to zero
    
    neg_section = dpillar[np.where(dpillar < 0)] # subset negative values (linear, will always be together at array bottom )
    if len(neg_section) == 0: # if no negative values are present
        return np.nan, np.nan, np.nan, np.nan # return non-numbers
    else:
        low_idx = len(neg_section) - 1 # otherwise use subset length to find lower boundary index
        
        xa = xpillar[low_idx] 
        xb = xpillar[low_idx + 1]
        ya = ypillar[low_idx]
        yb = ypillar[low_idx + 1]
        return xa,xb,ya,yb # return boundary points
    
def generate_interpolated_4Dfield(cube, target_heights, time_idx):
    """
    Parameters
    ----------
    cube : Iris cube with derived altitude auxiliary required
    target_heights : list or tuple of desired altitudes to be interpolated on
    time_idx : list of desired time indices
    Returns interpolated diagnostic field as numpy array.
    
    """
    data = cube.data
    alts = cube.coords('altitude')[0].points
    times = []
    for t_idx in time_idx:
        layers = []
        print('T=',t_idx)
        for xc in target_heights:
            equ_lat_row = []
            print('H=',xc)
            for y_idx in range(np.shape(data)[2]): # looping throw rows of equal grid latitude
                equ_lon_col = []
                for x_idx in range(np.shape(data)[3]): # looping hrough coloumns of equal grid longitude
                    xa,xb,ya,yb = find_boundpoints(data, alts, x_idx,y_idx, t_idx, xc) # find boundary values in model grid
                    yc = basic_interp(xa,xb,xc,ya,yb) # perform linear interpolation based on altitude for model data
                    equ_lon_col.append(yc) # append interpolated result to coloumn of equal grid longitude
                equ_lat_row.append(equ_lon_col) # append values to rows to form 2D field of interpolated values
            layers.append(equ_lat_row) # append 2D field to list to form 3D space of interpolated values
        times.append(layers)
    print(np.shape(times))
    return np.array(times)

def coloumn_alt_interpolation(cube, target_heights):
    data = cube.data
    alts = cube.coords('altitude')[0].points
    interp_values = []
    for idx in range(len(data[0,:]-1)):
        xa,xb,ya,yb = find_boundpoints(data, alts, xc = target_heights[idx],coloumns = True, c_idx = idx)
        yc = basic_interp(xa,xb,target_heights[idx],ya,yb) 
        interp_values.append(yc)
    return np.array(interp_values)
        
        
def interpolate_simple_series(cube, target_height, time_idx, grid_axis = 'grid_longitude', mode = '2d'):
    """
    Parameters
    ----------
    cube : cube of single vertical slice required with derivedd altitude coordinate field.
    target_height : desired altitude to interpolate onto (float or integer).
    time_idx : desired time index.
    
    """
    data = cube.data
    alts = cube.coords('altitude')[0].points
    iterlen = len(cube.coords(grid_axis)[0].points)
    
    if np.isscalar(target_height): # do this if target_heights argument is a number, not array or list
        series = []
        for i in range(iterlen):
            xa,xb,ya,yb = find_series_boundpoints(data, alts, i, time_idx, target_height) # find boundary values in model grid
            yc = basic_interp(xa,xb,target_height,ya,yb) # perform linear interpolation based on altitude for model data
            series.append(yc)
        return np.array(series), cube.coords(grid_axis)[0].points

def interpolate_matched_series(cube, obs_alts, obs_lons, time_idx, obs_lats = None, grid_axis = 'grid_longitude', mode = '2d'):
    """
    

    Parameters
    ----------
    cube : TYPE
        DESCRIPTION.
    obs_alts : TYPE
        DESCRIPTION.
    obs_lons : TYPE
        DESCRIPTION.
    time_idx : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    if mode  == '2d':
        sample_points = [('grid_longitude', obs_lons)]
    elif mode == '3d':
        sample_points = [('grid_longitude', obs_lons), ('grid_latitude', obs_lats)]
    newcube = trajectory.interpolate(cube,sample_points, method =  'nearest')
    #print(newcube)
    data = newcube.data
    alts = newcube.coords('altitude')[0].points
    series = []
    for i in range(len(newcube.coords(grid_axis)[0].points)):
        target_height = obs_alts[i]
        xa,xb,ya,yb = find_series_boundpoints(data, alts, i, time_idx, target_height) # find boundary values in model grid
        yc = basic_interp(xa,xb,target_height,ya,yb) # perform linear interpolation based on altitude for model data
        series.append(yc)
    return np.array(series), newcube.coords(grid_axis)[0].points

def interpolate_matched_series_4D(cube, obs_alts, obs_lons, obs_lats, obs_times, t_multiplier = 0):
    

   
    cube.coords('time')[0].points = cube.coords('time')[0].points - cube.coords('time')[0].points[0] + 0.5*(1 + t_multiplier)
    sample_points = [('grid_longitude', obs_lons), ('grid_latitude', obs_lats)]
    tempcube = trajectory.interpolate(cube,sample_points, method =  'nearest')
    newcube = trajectory.interpolate(cube,[('time',obs_times/(60*60))], method =  'linear')
    return(newcube)
    
        
            
def interpolate_timeseries(cube, tpoint_x, tpoint_y, tpoint_h):
    
    sample_points = [('grid_longitude', [tpoint_x]),('grid_latitude', [tpoint_y])]
    newcube = trajectory.interpolate(cube, sample_points, method = 'nearest')
    data = newcube.data
    alts = newcube.coords('altitude')[0].points
    series = []
    for t in range(len(newcube.coords('time')[0].points)):
        xa,xb,ya,yb = find_series_boundpoints(data, alts, 0, t, tpoint_h)
        yc = basic_interp(xa,xb,tpoint_h,ya,yb)
        series.append(yc)
    return np.array(series), newcube.coords('time')[0].points  

def interpolate_timeseries_coloumns(cube, tpoints_x, tpoints_y, method = 'nearest'):
    """
    

    Parameters
    ----------
    cube : An iris cube
    tpoints_x : Must be list/array of floats equal in length to toints_y
    tpoints_y : Must be list/array of floats equal in length to toints_x
    method : TYPE, optional
        DESCRIPTION. The default is 'nearest'.

    Returns
    -------
    newcube : TYPE
        DESCRIPTION.

    """
    sample_points = [('grid_longitude', tpoints_x),('grid_latitude', tpoints_y)]
    newcube = trajectory.interpolate(cube, sample_points, method = method)
    return newcube
        
def interpolate_timeseries_coloumns1D(cube, tpoints, axis, method = 'nearest'): 
     
    sample_points = [(axis, tpoints)]
    newcube = trajectory.interpolate(cube, sample_points, method = method)
    return newcube