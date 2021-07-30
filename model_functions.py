# -*- coding: utf-8 -*-
"""
Created on Sun Nov 10 14:35:49 2019

@author: Wilhelm Hodder

Function script to import for analysis of model data
"""


import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import iris
import iris.coord_categorisation

#import CubeCrossSectioner_UK as ccs

# Global settings
font = {'size' : 18}
matplotlib.rc('font', **font)

def file_info(path):
    """
    Print cube list from file
    """
    cube_list = iris.load(path)
    print(cube_list)
    
def cube_info(path, name):
    """
    Loads single specified cube from file defined by path,
    and prints standard iris cube info.
    """
    cube = iris.load_cube(path, name)
    print(cube)
    
def level_info(path, name, level, time = 0):
    """
    Loads cube specified at particular pressure/model level and prints info.
    """
    cube = iris.load(path, name)
    print(cube[0][time][level])
    
def get_model_levels(path, name, coord_name = 'altitude'):
    """
    Load cube with iris and return array of hybrid/model heights
    """
    cube = iris.load(path, name)
    height = cube[0].coord(coord_name).points
    return height

def get_coord_units(path, name, coord_name, return_unit = False):
    """
    Load cube with iris and return units or just print them.
    """
    cube = iris.load(path, name)
    units = cube[0].coord(coord_name).units
    if return_unit == True:
        return units
    else:
        print(units)
        
    

def return_cube_components(path, stash, lat_name = 'grid_latitude', lon_name = 'grid_longitude', print_info = True, **kwargs):
    """
    Loads cubes from data file using iris, 
    and extracts all data as seperate numpy arrays.
    Returns n-dimensional data, one dimensional latidue and longitude, 
    in that order.
    
    path - string: path to file and name
    stash - string: the variable/stash name, as specified in its cube: 
                use file_info() to find out the correct variable names
    """
    cube_list = iris.load(path, stash)
    # iris.load loads a list of cubes with one element. To access cube, index 0.
    cube = cube_list[0] 
    # print information 
    if print_info == True:
        print(cube)
    # load components from cube
    latitude = cube.coord(lat_name).points
    longitude = cube.coord(lon_name).points
    data = cube.data
    
    # if requested in **kwargs, extract additional coordinate arrays/keywords
    keyword_arrays = []
    for key, value in kwargs.items():
        if value == True:
            keyword_arrays.append(cube.coord(key).points)
    # return primary data and coordinates
    if len(keyword_arrays) < 1:
        return data, latitude, longitude
    else:
        return data, latitude, longitude, keyword_arrays

def return_concatenated_cube_components(path, file_letter, resolution, stash, file_number = 4, axis = 0, 
                                        lat_name = 'grid_latitude', lon_name = 'grid_longitude', short_file = True, **kwargs):
    """
    Same as return_cube_components(), but loads multiple files and concatenates
    data along time axis.:
        
    Loads cubes from data file using iris, 
    and extracts all data as seperate numpy arrays.
    Returns n-dimensional data, one dimensional latidue and longitude, 
    in that order.
    
    path - string: here only the path upto the directory in which the desired file is contained. 
            Make sure the final forward slash is included:
            e.g. '../Model_Data/u-bk574/nc/'
    file_letter - string: the letter that identifies the file from the UM model output streams:
                    f, g, h, i, j
    resolution - string: resolution of grid as specified by file name, including units:
                    e.g. '1p5km', '4p4km', '300m'
    file_number - integer: number of files to concatenate
    axis - integer: used for the np.concatenate() function.
            only redefine this if axis = 0 does not concatenate along time axis.
    **kwargs can include: t, Hybrid height,... - set to True to receive addiotnal arrays contained within list;
                            t array will be concatenated and comes as last(!) element.
    
    """
    data = []; keyword_arrays = []; t_coor = []
    for i in range(file_number):
        print('i=',i)
        # define file to load
        time = 6*i
        if time < 12:
            file_time = f'0{time}'
        else:
            file_time = time
        
        complete_path = f'{path}20180312T0000Z_IcelandGreenlandSeas_{resolution}_RA1M_p{file_letter}0{file_time}.nc'
        if short_file == True:
            complete_path = f'{path}{resolution}_{stash}_{i+1}.nc'
        if i == 0: # first step extracts coordinates as well as first segment of data
            data_segment, latitude, longitude = return_cube_components(complete_path, stash, 
                                                                       lat_name = lat_name, 
                                                                       lon_name = lon_name, 
                                                                       print_info = False)
            data = data_segment
            
            for key, value in kwargs.items():  # if requested in **kwargs, extract additional coordinate arrays/keywords
                #print(key, value)
                try:
                    if value == True:
                        print(key, value)
                        # use above function to apply variable keywords 
                        dump0, dump1, dump3, temporary = return_cube_components(complete_path, stash, 
                                                                           lat_name = lat_name, 
                                                                           lon_name = lon_name, 
                                                                           print_info = False,
                                                                           key = value)
                        # forecast time coordinate requires concatenation; single out
                        if key == 't':
                            t_coor = temporary[0]
                        else:
                            keyword_arrays.append(temporary[0])
                except iris.exceptions.CoordinateNotFoundError:
                    try:
                        keyword_arrays.append(temporary)
                    except:
                        print(f'Warning: specified argument not found: {key}')
        else:
            #print('i =', i)
            data_segment, dump0, dump1 = return_cube_components(complete_path, stash, 
                                                                lat_name = lat_name, 
                                                                lon_name = lon_name, 
                                                                print_info = False)
            # now concatenate data ndarrays
            data = np.concatenate((data, data_segment), axis = 0)
            
            for key, value in kwargs.items():  # if requested in **kwargs, extract additional coordinate arrays/keywords
                # specific for time coordinate concatenation
                print(key)
                try:
                    if value == True: 
                            if key == 't':
                        # use above function to apply variable keywords 
                                dump0, dump1, dump3, temporary = return_cube_components(complete_path, stash, 
                                                                                        lat_name = lat_name, 
                                                                                        lon_name = lon_name, 
                                                                                        print_info = False,
                                                                                        t = True)
                                # concatenate time cooridnate
                                t_coor = np.concatenate((t_coor, temporary[0]))
                                if i == file_number - 1:
                                    keyword_arrays.append(t_coor)
                except iris.exceptions.CoordinateNotFoundError:
                    print(f'Warning: specified argument not found: {key}')
            
            if i == file_number - 1: # last file
                
                # two return instances: one for kwargs and one without
                if len(keyword_arrays) > 0:
                    return data, latitude, longitude, keyword_arrays
                else:
                    return data, latitude, longitude
                                
def data_extrema(data, print_max_min = True):
    """
    Get max and minimum values from any data array.
    """
    w_max = np.nanmax(data)
    w_min = np.nanmin(data)
    if print_max_min == True:
        print(w_max, w_min)
    # which is bigger in magnitude?
    if w_max ** 2 > w_min ** 2:
        # set_wmax and set_wmin to be used in normalising colourmap
        set_max = w_max
        set_min = -w_max
    else:
        set_max = -w_min
        set_min = w_min 
    return set_min, set_max

def subset_1d(data, lower, higher):
    """ 
    subset 1d array, bounds set in gridpoint integer
    """
    return data[lower:higher]

def unrotate_coords(rot_lon, rot_lat, NPoleLon, NPoleLat):
    """
    Takes in 1d arrays of rotated longitude and latitude,
    and unrotates it using iris. Output are two 2d coordinate fields.
    """
    from iris.analysis.cartography import unrotate_pole
    rot_lon_2d, rot_lat_2d = np.meshgrid(rot_lon, rot_lat)
    lons, lats =  unrotate_pole(rot_lon_2d, rot_lat_2d,
                                NPoleLon, NPoleLat)
    return lons, lats




def subset_2d(data, South, North, West, East):
    """
    Function to subset 2D data set.
    South, North, West, and East are integers indicating the grid bounds of the subset
    """
    sub_data = []
    sub_data_rows = data[South:North]
    for row in range(len(sub_data_rows)): # subset row by row
        sub_data.append(sub_data_rows[row][West:East])
        
    return np.array(sub_data)

def subset_3d(data, South, North, West, East):
    """
    """
    sub_data = []
    sub_data_rows = data[:, South:North]
    for time in range(len(data)):
        sub_data_time = subset_2d(data[time], South = South, North = North, West = West, East = East)
        sub_data.append(sub_data_time)
    return np.array(sub_data)

def subset_4d(data, South, North, West, East):
    """
    Function to subset 4D data set.
    South, North, West, and East are integers indicating the grid bounds of the subset
    """
    sub_data = []
    for time in range(len(data)):
        sub_data_time = [] # subset data at specified time
        for level in range(len(data[time])):
            # singled out horizontal slice, will subset row by row within
            sub_data_time_level = subset_2d(data[time][level], South = South, North = North,
                                            West = West, East = East) # subset data at specified time and level
            sub_data_time.append(sub_data_time_level)
        sub_data.append(sub_data_time)
    return np.array(sub_data)

def filefriendly_time_from_index(index, increment):
    # formatting time in filename friendly format as string from index
    if increment == 'hourly':
        hour = index
        minute = 0
    if increment == 'half-hourly':
        hour = index//2
        minute = (index/2 - hour)*60
    minute = int(minute); hour = int(hour)
    return f'H{hour}M{minute}'
    
    
    
    


