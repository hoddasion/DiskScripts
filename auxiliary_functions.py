# python file for auxiliary functions to use in the data analyis of my observational data
# Wilhelm Hodder, 10th Jan 2019

import numpy as np
#from netCDF4 import Dataset
#import matplotlib
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt
from time import strftime
import math as math
import pandas as pd
from math import trunc
import sys # for exit function
 

# set flight and time
datetime = strftime("%Y%m%d%H%M%S")
flight = '306'; flight_date = '20180319'


def filter_chunck(lower, upper, data_array, x_array):
    # functions to remove single defined chunck of data from array
    # returns filtered data plus the the adjusted array for the x-variable (e.g. time)
    lower_chunck = data_array[np.where(x_array < lower)]
    upper_chunck = data_array[np.where(x_array > upper)]
    new_data = np.concatenate((lower_chunck, upper_chunck))
    lower_remain_x = x_array[np.where(x_array < lower)]
    upper_remain_x = x_array[np.where(x_array > upper)]
    new_x = np.concatenate((lower_remain_x, upper_remain_x))
    return new_data, new_x

def filter_spikes_306_1Hz(data_array, time_array):
    # function that removes manulally selected spike from TAT data
    # flight 306, 1Hz file
    data_array, time_array = filter_chunck(779.55, 779.74, data_array, time_array)
    data_array, time_array = filter_chunck(782.2, 782.3, data_array, time_array)
    data_array, time_array = filter_chunck(789.583, 789.633, data_array, time_array)
    data_array, time_array = filter_chunck(790.967, 791.016, data_array, time_array)
    data_array, time_array = filter_chunck(803.07, 803.15, data_array, time_array)
    data_array, time_array = filter_chunck(803.33, 803.38, data_array, time_array)
    data_array, time_array = filter_chunck(815.07, 815.13, data_array, time_array)
    data_array, time_array = filter_chunck(815.60, 815.65, data_array, time_array)
    data_array, time_array = filter_chunck(816.10, 816.15, data_array, time_array)
    data_array, time_array = filter_chunck(825.467, 825.56, data_array, time_array)
    data_array, time_array = filter_chunck(825.65, 826.52, data_array, time_array)
    data_array, time_array = filter_chunck(828.12, 828.182, data_array, time_array)
    data_array, time_array = filter_chunck(876.22, 876.40, data_array, time_array)
    data_array, time_array = filter_chunck(903.50, 903.75, data_array, time_array)
    data_array, time_array = filter_chunck(939.375, 939.45, data_array, time_array)
    data_array, time_array = filter_chunck(985.82, 986.04, data_array, time_array)
    data_array, time_array = filter_chunck(1034.28, 1034.34, data_array, time_array)
    data_array, time_array = filter_chunck(1034.47, 1034.58, data_array, time_array)
    data_array, time_array = filter_chunck(1046.480, 1046.570, data_array, time_array)
    
    new_data = data_array
    return new_data, time_array

def extract_chunk(lower, upper, data_array, x_array):
    # function to extract single defined chunck of data from array
    # returns desired data chunk plus the the adjusted array for the x-variable (e.g. time)
    chunk = data_array[np.where(x_array > lower)]
    new_x = x_array[np.where(x_array > lower)]
    chunk = chunk[np.where(new_x < upper)]
    new_x = new_x[np.where(new_x < upper)]
    return chunk, new_x

def calc_specific_hum(f):
    # function to calculate specific humidity from H2O molar fraction
    # returns array if array is supplied
    m_h2o = 18.02 # water vapour molecular weight
    m_air = 28.97 # dry air average molecular weight
    
    specific_hum = np.array(m_h2o*f/(m_air*(1-f)+m_h2o*f)) 
    return specific_hum
    
def figcoords_convert(value, lower, upper):
    # function to convert from data to figure scaled coord 
    # e.g. from lat and long coords to figure scale 0->1.
    return (value - lower)/(upper - lower)

def calc_long_in_nm(latitude, delta_longitude):
    # function to convert degrees longitude in nautical miles
    lat = latitude * (math.pi/180)
    lon = delta_longitude * (math.pi/180)
    
    C = 21600 # Earth's polar circumference in nm
    return math.cos(lat) * (lon/(2*math.pi)) * C 


def calculate_scalar_distance(lat_temp, lon_temp):
 # resolve coordinates (vector) into scalar distance
    x_or = 18.0756; y_or = 65.6546 #origin coordinates  
    distance = []  
    for l in range(len(lat_temp)):
        if l == 0:
            distance.append(0)
#         #distance.append(np.sqrt(pow(calc_long_in_nm(lat_temp[0], lon_temp[0] - x_or), 2) + pow((lat_temp[0] - y_or)*60, 2)))
        else:
            distance.append(distance[l-1] + np.sqrt(pow(calc_long_in_nm(lat_temp[l], lon_temp[l] - lon_temp[l - 1]), 2) + pow((lat_temp[l] - lat_temp[l - 1]) * 60, 2)))
    return distance   

def make_ticklist(data, N, dprec = 1):
    """
    to produce a list of evenly spaced floats from data to be used as tick values in plotting
    with a specified number of ticks
    data = an array/series of floats or integers
    N = number of desired ticks (elements)
    dprec = decimal precision of output values.
    """
    ticklist = []
    start = np.nanmin(data) # first value
    end = np.nanmax(data) # last value
    if dprec == 0: # ensure even spacing
        start = np.around(start); end = np.around(end)
    gapsize = (end - start) / (N - 1)
    # now start putting stuff in ticklist
    ticklist.append(round(start,dprec))
    for i in range(N-1):
        ticklist.append(round(start + (i + 1)*gapsize,dprec))
    #ticklist.append(round(end,dprec))
    return ticklist


def argnearest(items,pivot):
    """
    A function to find the nearest item to a given point. Very handy with coordinates
    Courstesy of Callum Rollo
    """
    near_item=min(items, key=lambda x: abs(x - pivot))
    for i in range(len(items)): 
        if items[i]==near_item:
            return i

def north_angle(u,v, radians = True):
    """
    tl:dr it returns the 4 quadrant arctan  of u and v
    but corrected as wind direction.
    """
    dir = np.arctan2(v,u)
    
    # turn around for convention of wind direction
    dirnew = np.copy(dir)
    dirnew[dir < 0] = dir[dir < 0] + np.pi 
    dirnew[dir >= 0] = dir[dir >= 0] - np.pi 
    if radians == True:
        return dirnew
    elif radians == False:
        return 180/np.pi * dirnew
    else:
        print('Value of \'radians\' keyword argument must be boolean.')
        
def variance_distrn(data, iterations, step_factor = 1, disclast = True):
    """
    Takes a 1D array of data and iteratively slices it, takes the variance 
    of each slice, and then the mean of the variances for each step in the loop.
    
    Returns a 1D numpy array of these variance means.
    
    data = 1D numpy array of floats or integers
    iterations = number of iterations to run the loop for
    step_factor = the scaling for which to increase the number of slices at each step. Default is 1; leads to linear increase (1,2,3,4,...)
    NOTE: step_factor and iterations must be integers
    disclast = discard last variance/slice at each step
    """
    mean_variances = []
    for i in range(iterations):
        w = data
        counter = 0 # for counting number of nans added (affects accuracy of last variance in run)
        while len(w) % ((i + 1)*step_factor) != 0: # number of elements has to be evenly split into equal chunks
            newvalue = np.array(np.nan)
            w = np.concatenate((w, newvalue), axis = None) # append a non-value
            counter += 1
        temp_sliced = w.reshape((i + 1)*step_factor, len(w)//((i+1)*step_factor)) # slice the data
        print('Number of slices =', (i+1)*step_factor)
        print('Counter =', counter)
        variances_temp = np.nanstd(temp_sliced, axis = 1)**2 # element-wise square of standard deviation
        if disclast == True and i != 0:
            variances = variances_temp[0:-1]
        else:
            variances = variances_temp
        mean_of_vars = np.nanmean(variances)
        mean_variances.append(mean_of_vars)
    return np.array(mean_variances)