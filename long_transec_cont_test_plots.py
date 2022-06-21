# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 15:06:44 2022

@author: kse18nru
"""
#%% module imports
import iris
from iris.analysis.cartography import rotate_pole
import data_management as dmna
import pandas as pd
import numpy as np
import persinterpolate as pert
import matplotlib
import matplotlib.pyplot as plt
import xsection_workshop as workshop
import iris.coord_categorisation
import iris.plot as iplt
import matplotlib.gridspec as gridspec



#%%
plt.rcParams['font.size'] = 28
resolutions = ['0p5km']#,'1p5km', '4p4km']

cubes = iris.load('D:/Project/Model_Data/u-cf117/LONGTAIL_0p5km_umpg1_flt306.nc')
print(cubes)