# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 12:09:20 2021

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
from iris.analysis import trajectory
#%% Global settings
plt.rcParams.update({'font.size': 20})
res = '0p5km'

time_series_plot_atmostate = True
if time_series_plot_atmostate:
    df = pd.read_csv('../../Figures/PNG/301/u-bu807/SlowBubble/tpoint_coors_301.csv')
    u_3dcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_{res}_um_m01s15i002_24hrs_h_301.nc', 'm01s15i002')
    v_3dcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_{res}_um_m01s15i003_24hrs_h_301.nc', 'm01s15i003')
    wsp_3dcube = (u_3dcube**2 + v_3dcube**2)**0.5
    w_3dcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_{res}_um_upward_air_velocity_24hrs_h_301.nc', 'upward_air_velocity')
    th_3dcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_{res}_um_air_potential_temperature_24hrs_i_301.nc', 'air_potential_temperature')
    q_3dcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_{res}_um_specific_humidity_24hrs_i_301.nc', 'specific_humidity')
    for i in range(9):
        label = df['label'][i]
        glon = df['grid_lon'][i]
        glat = df['grid_lat'][i]
        alt = df['alt_m'][i]
        
        print(glon)
        #%%
        sample_points = [('grid_longitude', [glon]),('grid_latitude', [glat])]
        wsp_traj = trajectory.interpolate(wsp_3dcube, sample_points, method = 'nearest')
        print(wsp_traj)
        #%%
        
        
    