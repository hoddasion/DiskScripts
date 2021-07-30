# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 13:04:44 2020

@author: kse18nru

"""

#%% Imports

import numpy as np
from scipy import interpolate

#%%

x = [1,2,3,4,5]
y = [3,1,4,2,5]

#%% Generate sample data sets

x_points = np.arange(20)
y_points = np.arange(20)
y_sin = np.sin(y_points)

f_lin = interpolate.interp1d(x_points,y_points)
f_sin = interpolate.interp1d(x_points,y_sin)
#print(f_sin([3.14 * 0.5]))

#%%
x_mesh, y_mesh = np.meshgrid(x_points, y_points)

z = np.sin(y_mesh) + 1

f_mesh = interpolate.interp2d(x_mesh,y_mesh,z)
print(f_mesh([12.2,15.7,15.8,10.11],[1.7,8.5]))
