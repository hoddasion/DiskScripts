# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 15:21:31 2020

@author: Wilhelm Hodder

Supporting functions for foundry adapted from Peter Sheridan's CubeCrossSectioner_UK + more that I may write.
These functions can also  be used outside of foundry in scripts.
"""

import numpy as np
import iris
from iris.analysis import trajectory
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
from matplotlib import mlab,colors,cm
from matplotlib.font_manager import FontProperties
import matplotlib.patheffects as PathEffects
from scipy import interpolate,stats
import scipy
import os
import sys
from copy import deepcopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from math import sin, cos, sqrt, atan2, log10
import datetime
from scipy.interpolate import interp2d




def roundit(floatNumber,updown):
    """ 
    rounds numbers up or down as specified
    """
    if updown == 'up':
        if round(floatNumber) < floatNumber:
            return int(round(floatNumber)+1)
        else:
            return int(round(floatNumber))
    elif updown == 'down':
        if round(floatNumber) > floatNumber:
            return int(round(floatNumber - 1))
        else:
            return int(round(floatNumber))

def find_closest_point(loclatlondict,xvals,orogslice):
    ''' 
    Crude minimum distance of a gridded line from a point, returns position along the line in 
    distance and point number  
    '''
    iout,xout=[],[] 
    for location in loclatlondict.keys():
        lldist2=(orogslice.coord('grid_longitude').points-loclatlondict[location][0])**2 + \
                (orogslice.coord('grid_latitude').points-loclatlondict[location][1])**2
        iclosest=np.argmin(lldist2)
        xout+=[xvals[iclosest]]
        iout+=[iclosest]
    
    return iout,xout

def getslicecompt(uslice,vslice,lonstart,latstart,lonend,latend,resolvedir=None):
    '''
    Convert horizontal velocities into the component in-plane using the initial bearing of 
    a vertical slice, or defined by resolvedir 
    ''' 
    dlong=(lonend-lonstart)*np.pi/180.
    if resolvedir: slicedir=np.pi*resolvedir/180.
    else: slicedir=np.arctan2( np.sin(dlong)*np.cos(latend*np.pi/180.) , \
                               np.cos(latstart*np.pi/180.)*np.sin(latend*np.pi/180.) -
                               np.sin(latstart*np.pi/180.)*cos(latend*np.pi/180.)*cos(dlong) )

    udir=90.*np.pi/180.
    vdir=0.
    comptslice=deepcopy(uslice)
    comptslice.data=uslice.data*np.cos(udir-slicedir) + vslice.data*np.cos(vdir-slicedir)

    return comptslice

def tickspacecalc(extent):
    '''calculate tick separation'''
    extent=float(extent)
    nmin=4
    nmax=8#5*nmin/2
    powr=min([1+roundit(log10(nmin/extent),'down'),1+roundit(log10(nmax/extent),'down')])
    facs=(1.,2.,3.,5.)
    fac=facs[[i for i in range(len(facs)) if nmin <= extent*10**powr/facs[i] <= nmax][-1]]
    if fac == 3.: fac=facs[[i for i in range(len(facs)) if nmin <= extent*10**powr/facs[i] <= nmax][0]]
    div=fac/10**powr
    mindat=0.
    maxdat=extent
    mindat=div*roundit(mindat/div,'up')
    maxdat=div*roundit(maxdat/div,'down')
    ticks=np.arange(mindat,maxdat+div,div)

    return ticks

def haversine(lat1,lon1,lat2,lon2):
    '''distance between points on a great circle'''

    R = 6371.0088 # mean radius of earth (km) 
    
    lat1 = np.pi/180*lat1
    lat2 = np.pi/180*lat2
    lon1 = np.pi/180*lon1
    lon2 = np.pi/180*lon2
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = (np.sin(dlat/2))**2 + np.cos(lat1) * np.cos(lat2) * (np.sin(dlon/2))**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    distance = R * c

    return distance

def generate_flattop_ridge(x,y,a,b,c,h0):
    """

    Parameters
    ----------
    x : 1D numpy array
        Symmetrical 1-D x coordinate array
    y : 1D numpy array
        Symmetrical 1-D y coordinate array
    a : float
        ridge slope half-width
    b : float
        ridge semi-circular end radius
    c : float
        major length; flat top extends over a long-axis distance of 2(b+c)
    h0 : float
        uniform height of flat top of ridge

    Returns
    -------
    numpy ndarray
    
    Calculates a 2-D surface altitude field for a flat-top ridge, using the formula adapted
    from Gabersek and Durran (2004), Gap Flows through Idealized Topography. 
    Part I: Forcing by Large-Scale Winds in the Nonrotating Limit, Journal of the Atmospheric Sciences.

    """
    # generate uniform ndarray grid of right dimensions to work with
    grid = np.ones((len(x),len(y))) 
    s_grid = np.ones((len(x),len(y))) 
    for i in range(len(x)):
        for j in range(len(y)):
            # compute s function elementwise
            if np.absolute(x[i]) <= c:
                s_grid[i,j] = s_grid[i,j]*np.maximum(0, np.absolute(y[j]) - b)
            else:
                s_grid[i,j] = s_grid[i,j]*np.maximum(0, ((y[j]**2)+(np.absolute(x[i])-c)**2)**0.5 - b)
    
    for i in range(len(x)):
        for j in range(len(y)):
            if s_grid[i,j] <= 4*a:
                #print('skipped', grid[i,j])
                grid[i,j] = (h0/16)*((1 + np.cos((np.pi*s_grid[i,j])/(4*a)))**4)
            else:
                grid[i,j] = 0
    return grid

def generate_ridgegap(x, y, d, e, floor_fraction = 0, recenter = 0):
    """
    

    Parameters
    ----------
    x : 1-D numpy array
        x coordinate array
    y : 1-D numpy array
        y coordinate array
    d : float
        Width of the floor of the gap.
    e : float
        The horizontal distance over which  the sidewalls rise from the floor
        of the gap to the ridgeline.
    floor_fraction : float
        The fractional height of the gap floor compared to the ridge height. 
        Must be positive and below 1. The default is 0.
    recenter : float
        The amount by which the center of the gap is moved along the long axis 
        of the ridge away from the center of the ridge. A -ve (+ve) value moves the center
        towards the left (right) (or into the negative (positive) domain of x).

    Returns
    -------
    grid : numpy ndarray
        Calculates a 2-D surface fractional altitude field for a single gap in a flat-top ridge, using the formula adapted
    from Gabersek and Durran (2004), Gap Flows through Idealized Topography. 
    Part I: Forcing by Large-Scale Winds in the Nonrotating Limit, Journal of the Atmospheric Sciences.

    """
    x = x - recenter
    grid = np.ones((len(x), len(y)))
    
    for i in range(len(x)):
        ## compute g-function per value of x
        if np.absolute(x[i]) <= d/2:
            g = floor_fraction
        elif (np.absolute(x[i]) > d/2 and np.absolute(x[i]) <= e +  d/2):
            g = np.sin((np.pi*(np.abs(x[i])-d/2))/(2*e))*(1-floor_fraction) + floor_fraction
        else:
            g = 1
        ## apply g-function to grid row
        grid[i] = grid[i]*g
    return grid


def velcmap():
    """
    defines a colour map for velocities, 2/3 of the map is for positive velocities, 1/3 for negative
    """
    getmap=cm.ScalarMappable(cmap='jet').get_cmap()
    poscols=[getmap(3*i/2) for i in range(2*255/3)]
    getmap=cm.ScalarMappable(cmap='bone').get_cmap()
    negcols=[getmap(255-2*i) for i in range(255-2*255/3)]

    colmap=colors.LinearSegmentedColormap.from_list('velocity',negcols+poscols)

    return colmap


def wcmap():
    """
    defines a colour map for vertical velocities
    """
    getmap=cm.ScalarMappable(cmap='RdBu_r').get_cmap()
    getwhite=cm.ScalarMappable(cmap='Greys').get_cmap()
    cols=[getmap(i) for i in range(255)]
    cols[256*7/16:256*9/16]=[getwhite(0)]*len(cols[256*7/16:256*9/16])

    colmap=colors.LinearSegmentedColormap.from_list('vertvel',cols)

    return colmap

def add_hybrid_height(orography_cube, cubes):
    ''' Add hybrid height to a cube where it's missing - thanks Becky Stretton '''
    regrid_cache = {}
    for cube in cubes:
        orog = iris.fileformats.rules._ensure_aligned(regrid_cache, orography_cube, cube)
        new_coord = iris.coords.AuxCoord(orog.data,
                                         orog.standard_name,
                                         orog.long_name,
                                         orog.var_name,
                                         orog.units,
                                         attributes=orog.attributes)
        dims = [cube.coord_dims(src_coord)[0] for src_coord in orog.dim_coords]
        cube.add_aux_coord(new_coord, dims)
        cube.add_aux_factory(iris.aux_factory.HybridHeightFactory(cube.coord('level_height'),
                                                                  cube.coord('sigma'),
                                                                  cube.coord('surface_altitude')))
    return cubes


def calc_time_from_index(index, sample = 'halfhourly' ):
    """
    Calculate hour and minute from integer time index.
    """
    if sample == 'halfhourly':
        hour = index//2
        minute = (index/2 - hour)*60
        minute = int(minute); hour = int(hour)
        return hour, minute
    if sample == 'hourly':
        hour = index
        minute = 0
        return hour, minute
    if sample == 'threehourly':
        hour = index*3
        minute = 0
        return hour, minute
    
def save_xsection(fig, res, domain, variable_in_file_name, flight, level = False, surface_variable = False,
                  figure_path_pdf = None, figure_path_png = None, time = None, sampling_rate = 'halfhourly', test = False,
                  pressure_levels = False, verstash = False, PDF = True, PNG = True):
    """
    Function to fomart figure file names and save figure in specified directories.
    Requires time as index time step (integer).
    """
    if test == True:
        import time
        t = time.localtime()
        current_time = time.strftime("%H%M%S", t)
        fig.savefig(f'test_figures/{current_time}.png')
        return 0
    hour, minute = calc_time_from_index(time, sample = sampling_rate)   # first calculate hour and minute values
    if surface_variable == False:
        if pressure_levels == False:
            file_name = f'{res}_{domain}_{variable_in_file_name}_20180312-{flight}_TH{hour}M{minute}_ML{level}'
            pdf_path = f'{figure_path_pdf}/{domain}/{level}/{file_name}.pdf'
            png_path = f'{figure_path_png}/{domain}/{level}/{file_name}.png'
        if pressure_levels == True:
            if verstash == False:
                pressure_labels = ['200hPa','300hPa', '500hPa','800hPa','850hPa','950hPa']
                file_name = f'{res}_{domain}_{variable_in_file_name}_20180312-{flight}_TH{hour}M{minute}_PL{pressure_labels[level-1]}'
                pdf_path = f'{figure_path_pdf}/{domain}/{pressure_labels[level-1]}/{file_name}.pdf'
                png_path = f'{figure_path_png}/{domain}/{pressure_labels[level-1]}/{file_name}.png'
            if verstash == True:
                pressure_labels = ['100hPa','150hPa','200hPa','250hPa','300hPa','400hPa','500hPa','600hPa','650hPa','700hPa','750hPa','800hPa','850hPa','925hPa','950hPa','1000hPa']
                file_name = f'{res}_{domain}_{variable_in_file_name}_20180312-{flight}_TH{hour}M{minute}_verstash_PL{pressure_labels[level-1]}'
                pdf_path = f'{figure_path_pdf}/{domain}/verstash/{pressure_labels[level-1]}/{file_name}.pdf'
                png_path = f'{figure_path_png}/{domain}/verstash/{pressure_labels[level-1]}/{file_name}.png'
    if surface_variable == True:
        file_name = f'{res}_{domain}_{variable_in_file_name}_20180312-{flight}_TH{hour}M{minute}_surface'
        pdf_path = f'{figure_path_pdf}/{file_name}.pdf'
        png_path = f'{figure_path_png}/{file_name}.png'
    try:
        
        if PDF == True:
            if level == False:
                file_name = f'{res}_{domain}_{variable_in_file_name}_20180312-{flight}_TH{hour}M{minute}'
                fig.savefig(f'{figure_path_pdf}/{file_name}.pdf')
                 
            else: 
                fig.savefig(pdf_path)
                 
        if PNG == True:
            if level == False:
                file_name = f'{res}_{domain}_{variable_in_file_name}_20180312-{flight}_TH{hour}M{minute}'
                
                fig.savefig(f'{figure_path_png}/{file_name}.png') 
            else: 
               
                fig.savefig(png_path) 
            
    except:
        print(f'Saving figure unsuccessful; please check input arguments.\nAttempted saving to {pdf_path} \nand\n{png_path}')
        
def singleaxis_distance_xticklists(cube, axis, point_coor, divisor):
    """
    Reads in cube and calculates distances from single axis coordinates via Haversine function.
    Outputs two arrays:
        1. array of coordinate values
        2. array of corresponding cumulative distance value
    """
    if axis == 'longitudinal':
        coords = cube.coord('grid_latitude').points
        # calculate point to point distance in kilometers with haversine function
        pointtopoint_dist = haversine(coords[:-1], point_coor, coords[1:], point_coor) 
        cumulative_dist = np.cumsum(pointtopoint_dist) # take cumulative sum
        cumulative_dist = np.concatenate(([0], cumulative_dist)) # add zero to beginning of array; now same length as coords
        unique_cum_dist, unique_index = np.unique(cumulative_dist.astype(int), return_index = True) # remove duplicate rounded integer values
        unique_coords = coords[unique_index] # reduce coord array to match unique_cum_dist
        select_dist = unique_cum_dist[np.where(unique_cum_dist%divisor == 0)] # select all points with integer values divisible by 50
        select_coords = unique_coords[np.where(unique_cum_dist%divisor == 0)] # select corresponding rotated coordinate values
    if axis == 'latitudinal':
        coords = cube.coord('grid_longitude').points
        # calculate point to point distance in kilometers with haversine function
        pointtopoint_dist = haversine(point_coor, coords[:-1], point_coor, coords[1:]) 
        cumulative_dist = np.cumsum(pointtopoint_dist) # take cumulative sum
        cumulative_dist = np.concatenate(([0], cumulative_dist)) # add zero to beginning of array; now same length as coords
        unique_cum_dist, unique_index = np.unique(cumulative_dist.astype(int), return_index = True) # remove duplicate rounded integer values
        unique_coords = coords[unique_index] # reduce coord array to match unique_cum_dist
        select_dist = unique_cum_dist[np.where(unique_cum_dist%divisor == 0)] # select all points with integer values divisible by 50
        select_coords = unique_coords[np.where(unique_cum_dist%divisor == 0)] # select corresponding rotated coordinate values
    return select_coords, select_dist

def variableaxis_distance_xticklists(lon_array, lat_array, divisor = 50):
    """
    Reads in equal length of lon and lat coordinates and returns both:
        1. two arrays of coordinate values, lon and lat
        2. array of corresponding cumulative distance value

    """
    # user haversine function to compute point-to-point distances
    ptpdist = haversine(lat_array[:-1],lon_array[:-1], lat_array[1:], lon_array[1:])
    cumdist = np.cumsum(ptpdist) # elementwise cumulative sum
    cumdist = np.concatenate(([0], cumdist)) # add zero to beginning of array; now same length as coords
    unique_cum_dist, unique_index = np.unique(cumdist.astype(int), return_index = True) # remove duplicate rounded integer values
    unique_lons = lon_array[unique_index] # reduce coord array to match unique_cum_dist
    unique_lats = lat_array[unique_index]
    print(f'unique distance :\n{unique_cum_dist}')
    select_dist = unique_cum_dist[np.where(unique_cum_dist%divisor == 0)] # select all points with integer values divisible by 50
    select_lats = unique_lats[np.where(unique_cum_dist%divisor == 0)] # select corresponding rotated coordinate values
    select_lons = unique_lons[np.where(unique_cum_dist%divisor == 0)]
    select_index = unique_index[np.where(unique_cum_dist%divisor == 0)]
    return select_lons, select_lats, select_dist, select_index

"""
Observations interpolator:
    
Email from Andy Elvidge sent on 30th March 2020:
    
    Regarding your problem getting your obs on the same grid as your model data. 
    If I understand you correctly, one simple good-enough method which is independent of IRIS is as follows…
Interpolate the x,y model grid coordinates from the model lat/lon to each obs point lat/lon. 
I.e. where model_lon, model_lat and model_x are the lon, lat and x data on your 2-dimensional model grid, 
in scipy.interpolate.interp2d, something like this 
(though depending on your input dimensions you may first have to use meshgrid on your input lon and lat coords)…

f = interpolate.interp2d(model_lon, model_lat, model_x)
f(obs_point_lon,obs_point_lat)

And then do the same for model_y to get both x and y coords of your obs point.

   """ 
    
def interp_obs_coords(lons, lats, model_x, model_y, obs_lon, obs_lat):
    """
    Uses unrotated Earth coordinate mapping derived from model grid to map observation coordinates onto model's x-y grid via interpolation.
    """
    ## meshgrid model x-y arrays
    model_xx, model_yy = np.meshgrid(model_x,model_y)
    ## build 2d array as if lat-lon grid were flat, not curved
    flat_lons = np.tile(lons[0], (len(lons[:,0]),1)) # taking edge of curved lons grid, and copy over all points in y
    ## lats need to be reshaped to create an array of distinct 1-element arrays before expanding
    lats_border = lats[:,0]  #take edge of curved lats grid
    lats_border = np.reshape(lats_border, (len(lats_border),1))
    flat_lats = np.tile(lats_border, (len(lats[0])))
    # make lwarg dict
    keywords = {}
    ## before interpolation, which reorders the arrays, 
    ## we need to save the original order of elements for later using numpy argsort()
    x_idx = np.argsort(np.argsort(obs_lon)) # return list of indices from before ordering
    y_idx = np.argsort(np.argsort(obs_lat))
    
    ## interpolate using lats and lons edges (i.e. interpolating onto a flat lat-lon grid)
    fx_flat = interp2d(lons[0], lats[:,0], model_xx)
    fy_flat = interp2d(lons[0], lats[:,0], model_yy)
    # after interpolation, also resort arrays back into original order
    dbx_flat = fx_flat(obs_lon,obs_lat, **keywords) # input database observation coordinates
    dby_flat = fy_flat(obs_lon,obs_lat, **keywords) # this produces two 2d arrays from which the first index is what I want
    dbx_flat = dbx_flat[:,x_idx]
    dby_flat = dby_flat[y_idx,:]
    
    
    ## interpolate onto deviation grid using newly found model x-y coordinates
    lats_dev = lats - flat_lats # take grid deviation in latitude between flat and curved grid
    lons_dev = lons - flat_lons # take grid deviation in longitude between flat and curved grid
    fx_dev = interp2d(model_x,model_y,lons_dev) # now we're flipping the mapping 
    fy_dev = interp2d(model_x,model_y,lats_dev)
    # after interpolation, also resort arrays back into original order
    dblon_dev = fx_dev(dbx_flat[0],dby_flat[:,0], **keywords) # and input the interpolated coordinates
    dblat_dev = fy_dev(dbx_flat[0],dby_flat[:,0], **keywords)
    dblon_dev = dblon_dev[:,x_idx]
    dblat_dev = dblat_dev[y_idx,:]
    #print(dblon_dev[:,0])
    ## now we use this to modify the measured observations coordinates
    ## and repeat the above process for one round of interpolations
    db_lon_crd = db_lon - dblon_dev[0] # calculate 'corrected' coordinate for database obs
    db_lat_crd = db_lat - dblat_dev[:,0]
    # after interpolation, also resort arrays back into original order
    dbx = fx_flat(db_lon_crd,db_lat_crd, **keywords)
    dby = fy_flat(db_lon_crd,db_lat_crd, **keywords)
    dbx = dbx[:,x_idx]
    dby = dby[y_idx,:]
    return dbx,dby
    