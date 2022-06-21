# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 11:53:18 2021

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
from iris.analysis import trajectory


#%% 
res = '1p5km'
flight = '306'
suite = 'u-cc134'#'u-bu807'
config = 'RA1M'
stream_prefix = 'p'
nc_directory = ''#'nc/Control/'

#%%
print('Loading cubes')
#theta_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_directory}RA1M_{res}_um_air_potential_temperature_24hrs_{stream_prefix}i_{flight}.nc', 'air_potential_temperature')
#q_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_directory}RA1M_{res}_um_specific_humidity_24hrs_{stream_prefix}i_{flight}.nc', 'specific_humidity')
#p_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_directory}RA1M_{res}_um_air_pressure_24hrs_{stream_prefix}i_{flight}.nc', 'air_pressure')
lh_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_directory}{config}_{res}_um_upward_heat_flux_in_air_24hrs_{stream_prefix}i_{flight}.nc', 'upward_heat_flux_in_air')
vap_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{nc_directory}{config}_{res}_um_upward_water_vapor_flux_in_air_24hrs_{stream_prefix}i_{flight}.nc', 'upward_water_vapor_flux_in_air')


#%% approximation
print('Calculating LH')
L = 2.25e6
lh_data = vap_cube.data * L  #dmna.calc_lhf(q_cube.data,theta_cube.data,p_cube.data,vap_cube.data, q_conversion = 1e-3, R = 287, L = 2.25e6)

#%%
print('Modding cube')
lh_cube.data = lh_data
lh_cube.standard_name = 'upward_latent_heat_flux_in_air'

#%%
print('Saving cube')
iris.save(lh_cube, f'D:/Project/Model_Data/{suite}/{nc_directory}{config}_{res}_um_upward_latent_heat_flux_in_air_24hrs_{stream_prefix}i_{flight}.nc')