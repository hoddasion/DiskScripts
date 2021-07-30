
"""
Created on Thu Jan  2 16:38:21 2020

Functions to assist with the generation of animations.

@author: Wilhelm Hodder
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import assist_functions.model_functions as modf
import itertools
import datetime
import matplotlib.animation as animation

def animate_pcolor(i, data, lon, lat, land_mask, level, norm_offset, cmap, ax, variable_name):
    """
    Generates a frame of pcolor when called in matplotlib's FuncAnimation(*args, **kwargs).
    First argument 'i' is passed to the iterator, while all following positional arguments must be passed in a tuple to the 'fargs' argument.
    """
    #data = sub_data_region2; lon = sub_lon2; lat = sub_lat2
    ## normalise the colourmap to be between a certain min and max
    data_min = np.nanmin(data[:][level]) - norm_offset;
    data_max = np.nanmax(data[:][level]) + norm_offset
    my_norm= matplotlib.colors.Normalize(vmin = data_min,vmax= data_max, clip=True)
    # formatting time string/datetime object
    time = i
    hour = time//2
    minute = (time/2 - hour)*60
    minute = int(minute); hour = int(hour)
    TOD = datetime.time(hour, minute)  # time of day
    # clear axis object, and plot new over it
    ax.clear()
    ax.contour(lon,lat,land_mask, colors = 'k', linewidths = 0.5)
    ax.pcolor(lon, lat, data[i][level], cmap = cmap, norm = my_norm, rasterized = True)
    ax.set(Title = f'{variable_name} at model level {level + 1} at {TOD}UTC') # set title with time stamp
    # return axis object
    return ax
    
def animate_pcolor_quiver(i, u, v, n, lon, lat, land_mask, level, norm_offset, cmap, ax):
    """
    Generates a frame of quiver plot with pcolor base when called in matplotlib's FuncAnimation(*args, **kwargs).
    First argument 'i' is passed to the iterator, while all following positional arguments must be passed in a tuple to the 'fargs' argument.
    """
    magnitude = np.sqrt(u[:,level]**2 + v[:,level]**2) # standard magnitude calculation
    ## sub sample u and v for quiver every nth data point in both horizontal dimensions
    
    u_level = u[i][level]#/np.sqrt(u[:][level]**2 + v[:][level]**2)
    v_level = v[i][level]#/np.sqrt(u[:][level]**2 + v[:][level]**2)
    u_n = u_level[::n]
    v_n = v_level[::n]
    u_n2 = []; v_n2 = []
    for j in range(len(u_n)):
        u_n2.append(u_n[j][0::n])
        v_n2.append(v_n[j][0::n])
    # normalise u and v now that they are in the correct shape
    u_n2 = np.array(u_n2); v_n2 = np.array(v_n2)
    #u_n2 = u_n2/np.sqrt(u_n2**2 + v_n2**2)
    #v_n2 = v_n2/np.sqrt(u_n2**2 + v_n2**2)
    
    ## normalise the colourmap to be between a certain min and max
    data_min = np.nanmin(magnitude) - norm_offset;
    data_max = np.nanmax(magnitude) + norm_offset
    my_norm= matplotlib.colors.Normalize(vmin = data_min,vmax= data_max, clip=True)
    # formatting time string/datetime object
    time = i
    hour = time//2
    minute = (time/2 - hour)*60
    minute = int(minute); hour = int(hour)
    TOD = datetime.time(hour, minute)  # time of day
    # clear axis object, and plot new over it
    ax.clear()
    ax.contour(lon,lat,land_mask, colors = 'k', linewidths = 0.5)
    ax.pcolor(lon, lat, magnitude[i], cmap = cmap, norm = my_norm, rasterized = True) # pcolor base
    ax.quiver(lon[::n], lat[:-1:n], u_n2, v_n2, pivot = 'middle', norm = my_norm) # plot quiver
    ax.set(Title = f'Wind speed and direction at model level {level + 1} at {TOD}UTC') # set title with time stamp
    return ax
    