# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 15:50:55 2021

@author: kse18nru
"""

#%% module imports
import sys
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
import model_foundry_toolkit as toolkit
import data_management as dmna
import time
import model_foundry as foundry
import matplotlib.gridspec as gridspec
import iris.analysis.cartography
import cartopy.crs as ccrs

#%% global variables


def glm_load_and_concat(stashcode, kconstraints, fileno = 4, flight = 306, suite = 'u-cc134'):
    if fileno < 2:
        return iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_glm_{fileno}_flt{flight}.nc', stashcode).intersection(**kconstraints)
    else:   
       
        files = []
        for i in range(fileno - 1):
            print('fileno =', i)
            filename = f'D:/Project/Model_Data/{suite}/RA1M_glm_{i+2}_flt{flight}.nc'
            print(filename)
            files.append(filename)
        
        from iris.experimental.equalise_cubes import equalise_attributes
        cubes =  iris.load(files, [stashcode])
        iris.util.unify_time_units(cubes)
        print(equalise_attributes(cubes))
        #cubes[0] = cubes[0][1:]
        for i in range(3):
            print(cubes[i].coords('time')[0])
        conccube = cubes.concatenate()[0].intersection(**kconstraints) # a complete cube spanning full forecast
        
       # print(iris.load(files, [stashcode]).concatenate_cube())
        return conccube

def extract_data_and_coords(cube):
    return cube.data,cube.coords('longitude')[0].points,cube.coords('latitude')[0].points