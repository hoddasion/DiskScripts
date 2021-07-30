# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 14:01:25 2019

@author: Wilhelm Hodder

Vertical_altitude_profile
"""
#%% module imports
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import iris
import iris.coord_categorisation
#import CubeCrossSectioner_UK as ccs
import iris.quickplot as qplt
import cartopy
import cartopy.feature as cfeat
import cartopy.crs as ccrs
import model_functions as modf
import itertools
import datetime

#%% global configurations 

import warnings
warnings.filterwarnings("ignore") ##BAD! If something doesn't work, re-enable warnings to check
font = {'size' : 18}
matplotlib.rc('font', **font)

#%% load cubes and data, then detach cube
variable = 'air_potential_temperature'
seventy_fourty_cube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/{variable}/1p5km_{variable}_24hrs.nc', variable)
ninety_fourty_cube = iris.load_cube(f'../../Model_Data/u-bk574/nc/Control/{variable}/0p5km_{variable}_24hrs.nc', variable)

seventy_heights = seventy_fourty_cube.coord('altitude').points
ninety_heights = ninety_fourty_cube.coord('altitude').points
seventy_fourty_cube = None
ninety_fourty_cube = None
print(ninety_heights.shape)
#%%
levels71 = np.arange(1,72)
levels91 = np.arange(1,92)
# sample means, mins and maxima for both altitude-level regimes
mean_71 = np.mean(seventy_heights, axis = 2)
mean_71 = np.mean(mean_71, axis = 1)
mean_91 = np.mean(ninety_heights, axis = 2)
mean_91 = np.mean(mean_91, axis = 1)
min_71 = np.min(seventy_heights, axis = 2)
min_71 = np.min(min_71, axis = 1)
min_91 = np.min(ninety_heights, axis = 2)
min_91 = np.min(min_91 , axis = 1)
max_71 = np.max(seventy_heights, axis = 2)
max_71 = np.max(max_71, axis = 1)
max_91 = np.max(ninety_heights, axis = 2)
max_91 = np.max(max_91 , axis = 1)
print(len(mean_71))
#%%
fig, axes = plt.subplots(1,1, figsize = (12,12))
q71mean = axes.plot(levels71, mean_71/1000, color = 'b')
q71min = axes.plot(levels71, min_71/1000, color = 'b', linestyle = 'dashed')
q71max = axes.plot(levels71, max_71/1000, color = 'b', linestyle = 'dashed')
q91mean = axes.plot(levels71, mean_91/1000, color = 'r')
q91min = axes.plot(levels71, min_91/1000, color = 'r', linestyle = 'dashed')
q91max = axes.plot(levels71, max_91/1000, color = 'r', linestyle = 'dashed')
axes.set(xlabel = 'Model theta level',
         ylabel = 'Model altitude [km]')
plt.savefig('../../Figures/PDF/301/ninety_seventy_thetalevel_altitude_flt301.pdf')
plt.savefig('../../Figures/PNG/301/ninety_seventy_thetalevel_altitude_flt301.png')
plt.close()

#%% print emans with model level through loop
## for reference
for i in range(len(mean_71)):
    print(i+1, mean_71[i])
for i in range(len(mean_91)):
    print(i+1, mean_91[i])
#%%
fig, axes = plt.subplots(1,1, figsize = (12,12))
n = 21
q71mean = axes.plot(levels71[:n], mean_71[:n], color = 'b')
q71min = axes.plot(levels71[:n], min_71[:n], color = 'b', linestyle = 'dashed')
q71max = axes.plot(levels71[:n], max_71[:n], color = 'b', linestyle = 'dashed')
q91mean = axes.plot(levels71[:n], mean_91[:n], color = 'r')
q91min = axes.plot(levels71[:n], min_91[:n], color = 'r', linestyle = 'dashed')
q91max = axes.plot(levels71[:n], max_91[:n], color = 'r', linestyle = 'dashed')
axes.set(xlabel = 'Model theta level',
         ylabel = 'Model altitude [m]')
xlabels = np.arange(0,21); print(xlabels)
plt.xticks(xlabels)
ylabels = np.arange(0,16)*200
plt.yticks(ylabels)
plt.grid(True)
#plt.savefig('../../Figures/PDF/301/ninety_seventy_thetalevel_altitude_sub20_flt301.pdf')
#plt.savefig('../../Figures/PNG/301/ninety_seventy_thetalevel_altitude_sub20_flt301.png')
#plt.close()
plt.show()
