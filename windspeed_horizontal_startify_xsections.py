# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 12:02:38 2020

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
import xsection_workshop as workshop
import stratify
from iris.experimental import iris_stratify as Stratify
from scipy.interpolate import interp1d
from functools import partial
#%% universals
heights = [1000] # in meters
#Stratify.
#%% concuct our own interpolator

#interpolator = partial(stratify.in)
#%% load data cube
variable1 = 'm01s15i002'
ucube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_0p5km_um_{variable1}_24hrs_h_301.nc', variable1)[:,:30]
print(ucube)

print(np.shape(ucube.coord('altitude')))
#u_geom = Stratify.relevel(ucube, ucube.coord('altitude'), heights, axis =1)
#print(u_geom)
variable2 = 'm01s15i003'
vcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_0p5km_um_{variable2}_24hrs_h_301.nc', variable2)[:,:30]

#maincube = (ucube**2 + vcube **2)**0.5
#print(maincube)

#%% stratify (interpolate vertical)
data = ucube.data
alts = ucube.coords('altitude')[0].points
print(np.shape(data))

print(np.shape(alts))
#%%
x_idx = 200; y_idx = 200; t_idx = 24

#%%
import persinterpolate as pert
#%%
target_heights = [500,1000]

layers = []
for xc in target_heights:
    equ_lat_row = []
    for y_idx in range(np.shape(data)[2]): # looping throw rows of equal grid latitude
        equ_lon_col = []
        for x_idx in range(np.shape(data)[3]): # looping hrough coloumns of equal grid longitude
            xa,xb,ya,yb = pert.find_boundpoints(ucube, x_idx,y_idx, t_idx, xc) # find boundary values in model grid
            yc = pert.basic_interp(xa,xb,xc,ya,yb) # perform linear interpolation based on altitude for model data
            equ_lon_col.append(yc) # append interpolated result to coloumn of equal grid longitude
        equ_lat_row.append(equ_lon_col) # append values to rows to form 2D field of interpolated values
    layers.append(equ_lat_row) # append 2D field to list to form 3D space of interpolated values
layers = np.array(layers) # convert to numpy array
