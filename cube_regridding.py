# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 15:41:14 2022

@author: kse18nru
"""
#%%
import iris
import iris.plot as iplt
from iris.analysis.cartography import rotate_pole
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from iris.analysis import trajectory
from iris.coords import DimCoord, AuxCoord
import timeit

#%% globals
flight = 306

#%% function defintions
def rotate_xy(x,y,phi, angular_units = 'rad'):
    
    if angular_units == 'degrees':
        phi = (np.pi/180)*phi
    x0 = 0
    y0 = 0
    x0 = np.cos(phi)*x - np.sin(phi)*y
    y0 = np.sin(phi)*x + np.cos(phi)*y
    return x0,y0

def map_grid(vertexA, phi, x_length, y_length,res_step = 0.005, angular_units = 'degrees', output_old_vectors = False, output_new_vectors = False):
    
    vertexAr = rotate_xy(vertexA[0], vertexA[1], phi, angular_units = angular_units)
    print(vertexAr)
    
    ## create vector A-B in rotated space
    # use np.ones() to fix x coord in rotated space
    # use linspace to evenly space points in y direction in rotated space
    vectorABr = np.array([np.ones(y_length)*vertexAr[0], np.linspace(vertexAr[1], vertexAr[1]+res_step*y_length, num = y_length)])

    ## create vector A-C in rotated space
    # use np.ones() to fix x coord in rotated space
    # use linspace to evenly space points in y direction in rotated space
    vectorACr = np.array([np.linspace(vertexAr[0], vertexAr[0]+res_step*x_length, num = x_length), np.ones(x_length)*vertexAr[1]])
    
    ## meshgrid vectors
    xvector = vectorACr[0]
    yvector = vectorABr[1]
    
    xmeshr, ymeshr = np.meshgrid(xvector,yvector)
    
    
    xmesh, ymesh = rotate_xy(xmeshr, ymeshr, -phi, angular_units = angular_units)
    if output_old_vectors:
        vectorAC = rotate_xy(vectorACr[0], vectorACr[1], -phi, angular_units = angular_units)
        vectorAB = rotate_xy(vectorABr[0], vectorABr[1], -phi, angular_units = angular_units)
        return xmesh, ymesh, vectorAC, vectorAB
    elif output_new_vectors:
        return xmesh, ymesh, vectorACr, vectorABr
    else:
        return xmesh, ymesh

def rotate_and_extract_subcube(cube, vertexA, phi, x_length, y_length,res_step = 0.005, angular_units = 'degrees'):
    
    xmesh, ymesh, vectorACr, vectorABr = map_grid(vertexA, phi,x_length, y_length, output_new_vectors = True )
    
    transect_cubes = iris.cube.CubeList([])
    
    for i in range(y_length):
        sample_points = [('grid_longitude', xmesh[i]),('grid_latitude', ymesh[i])]
        transect = trajectory.interpolate(cube, sample_points, method = 'nearest')
        
        rot_latitude = AuxCoord([vectorABr[1][i]], standard_name = 'latitude', units = 'degrees')
        rot_longitude = DimCoord(vectorACr[0], standard_name = 'longitude', units = 'degrees')
    
        transect.add_dim_coord(rot_longitude, 2)
        transect.add_aux_coord(rot_latitude)
        transect_cubes.append(transect)
    
    merged_cube = transect_cubes.merge()
    return merged_cube
#%%
if flight == 306:
    suite = 'u-cc134'
    config = 'RA1M'
    variable = 'air_potential_temperature'
    stream = 'pi'
    res = '0p5km'
    cube_path = f'D:/Project/Model_Data/{suite}/'
    filename = f'{config}_{res}_um_{variable}_24hrs_{stream}_{flight}.nc'
    
    
    cube = iris.load_cube(f'{cube_path}{filename}', variable)
    lbm_cube = iris.load_cube(f'{cube_path}{config}_{res}_um_land_binary_mask_24hrs_pi_{flight}.nc', 'land_binary_mask')
    print(cube)

    #%% extract current cube grid coords
    grid_lon = cube.coords('grid_longitude')[0].points
    grid_lat = cube.coords('grid_latitude')[0].points
    print(grid_lon, grid_lat)
    
    #%% plot something for reference
    tidx = 30
    hidx = 5
    
    fig, ax  = plt.subplots(1,1, figsize = (6,10))
    
    ax.contourf(grid_lon, grid_lat, cube[tidx, hidx].data, levels = 10)
    ax.contour(grid_lon, grid_lat,lbm_cube.data, colors = 'k')
    ax.grid()
    #plt.scatter(x,y, c = 'k', s = 100)
    #plt.contour(gridx)
    
    #%% make new coord grid as two 2D arrays
    phi = 40
    res_step = 0.005
    y_length = 200
    x_length = 170
    
    
    vertexA = (360.075, 0.77)
    
    
    xmesh, ymesh, vectorAC, vectorAB = map_grid(vertexA, phi,x_length, y_length, output_old_vectors = True )
    
    
    ##%% plot something for reference
    tidx = 30
    hidx = 5
    
    fig, ax  = plt.subplots(1,1, figsize = (6,10))
    
    #ax.contourf(grid_lon, grid_lat, cube[tidx, hidx].data, levels = 10)
    ax.contour(grid_lon, grid_lat,lbm_cube.data, colors = 'k')
    ax.grid()
    ax.scatter(vertexA[0], vertexA[1], label = 'A')
    ax.plot(vectorAB[0], vectorAB[1], color = 'k')
    ax.plot(vectorAC[0], vectorAC[1], color = 'k')
    ax.legend()
    #ax.scatter(meshAC, meshAB)
    #plt.scatter(x,y, c = 'k', s = 100)
    #plt.contour(gridx)
    
    #%% interpolation testing
    
    
    print(np.shape(xmesh))
    #new_cube = rotate_and_extract_subcube(cube, vertexA, phi, x_length, y_length,res_step = 0.005, angular_units = 'degrees')