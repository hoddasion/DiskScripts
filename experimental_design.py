# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 11:35:38 2021

@author: kse18nru
"""

#%% Imports
import sys
import iris
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#import model_foundry_toolkit as toolkit
import model_foundry_toolkit as toolkit
#%% Global settings
plt.rcParams.update({'font.size': 28})
#%% Load UM orography data
path_to_file = '../../Model_Data/u-bu807/nc/Control/surface_altitude/'
filename = '0p5km_surface_altitude_flt301.nc'

cube = iris.load(f'{path_to_file}{filename}')[0]
print(cube)

#%% subset model data
North = 370; South = 280
West = 100; East = 280
sub_cube = cube[South:North, West:East]

#%% extract model data
model_data = sub_cube.data
model_lon = sub_cube.coords('grid_longitude')[0].points

model_lat = sub_cube.coords('grid_latitude')[0].points

#%% convert grid coordinates to sclar distance coordinates (in km)
const_lon = np.ones(len(model_lat))*model_lon[0]
const_lat = np.ones(len(model_lon))*model_lat[0]
x_steps = toolkit.haversine(const_lat[:-1], model_lon[:-1], const_lat[1:], model_lon[1:])
y_steps = toolkit.haversine(model_lat[:-1], const_lon[:-1], model_lat[1:], const_lon[1:])
x_coor = np.cumsum(x_steps); x_asym = np.cumsum(x_steps)
x_coor = x_coor - x_coor[-1]/2
y_coor = np.cumsum(y_steps); y_asym = np.cumsum(y_steps)
y_coor = y_coor - y_coor[-1]/2
#%% take maximum along x-axis
zonal_maxima = sub_cube.collapsed('grid_longitude', iris.analysis.MAX)
model_zonal = zonal_maxima.data
print(model_zonal[20:-20])
mean_max_height = np.mean(model_zonal[20:-20])
print(mean_max_height)
h0 = mean_max_height
print(h0)
#%% generate flat top ridge
plainridge = toolkit.generate_flattop_ridge(x_coor,y_coor,a = 3,b = 5, c = 35, h0 = 650)
print(np.shape(plainridge))
print(np.shape(model_data[1:,1:]))
print(np.max(plainridge))
#%% generate ridge with gap
gap_params = {'d':3,'e':10,'floor_fraction':0.4}
gap_function = toolkit.generate_ridgegap(x_coor, y_coor, recenter = 0, **gap_params)
gapridge_one = plainridge*gap_function

#%% generate ridge with two gaps
gapf_left = toolkit.generate_ridgegap(x_coor, y_coor, recenter = -15, **gap_params)
gapf_right = toolkit.generate_ridgegap(x_coor, y_coor, recenter = 15, **gap_params)
gapf_comb = gapf_left*gapf_right
gapridge_two = plainridge*gapf_comb
#%% create contour value list
contour_values = (50, 150, 250, 350, 450, 550, 650)
kparams = {'levels' : contour_values, 'colors' : 'k'}
#%% plot figure with plain ridge
fig,ax = plt.subplots(1,1, figsize = (25,5))
ax.contour(x_coor,y_coor, np.rot90(plainridge),**kparams)
ax.set_xlabel('Horizontal West-East Distance [km]')
ax.set_ylabel('Horz. S-N Dist [km]')
plt.savefig('D:/Project/Figures/PNG/Idealised/plain_ridge_example.png')
plt.show()
#%% plot the figure with gapped ridges

fig,axes = plt.subplots(6,1, figsize = (25,25))
axes[0].contour(x_coor, y_coor, np.rot90(gapridge_one), label = 'A.', **kparams)
axes[0].set_xlim(left = -48, right = 48)
axes[0].set_xticks([])
axes[0].text(-45,15,'A.')
axes[0].set_ylabel('[km]')
axes[1].plot(x_coor, np.max(np.rot90(gapridge_one), axis = 0), label = 'B.', color = 'k')
axes[1].set_xlim(left = -48, right = 48)
axes[1].set_xticks([])
axes[1].set_ylim(top = 1250)
axes[1].text(-45,1000,'B.')
axes[1].set_ylabel('[m]')
axes[2].contour(x_coor, y_coor, np.rot90(gapridge_two), label = 'C.', **kparams)
axes[2].set_xlim(left = -48, right = 48)
axes[2].set_xticks([])
axes[2].text(-45,15,'C.')
axes[2].set_ylabel('[km]')
axes[3].plot(x_coor, np.max(np.rot90(gapridge_two), axis = 0), label = 'D.', color = 'k')
axes[3].set_xlim(left = -48, right = 48)
axes[3].set_xticks([])
axes[3].set_ylim(top = 1250)
axes[3].text(-45,1000,'D.')
axes[3].set_ylabel('[m]')
axes[4].contour(x_coor, y_coor, model_data[1:,1:], label = 'E.',**kparams )
axes[4].set_xlim(left = -48, right = 48)
axes[4].set_xticks([])
axes[4].text(-45,15,'E.')
axes[4].set_ylabel('[km]')
axes[5].plot(x_coor, np.max(model_data[1:,1:], axis = 0),label = 'F.', color = 'k')
axes[5].set_xlim(left = -48, right = 48)
axes[5].set_ylim(top = 1250)
axes[5].text(-45,1000,'F.')
axes[5].set_ylabel('[m]')
axes[5].set_xlabel('Horizontal West-East distance [km]')
plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)

plt.savefig('D:/Project/Figures/PNG/Idealised/gapped_ridge_examples.png')
plt.show()