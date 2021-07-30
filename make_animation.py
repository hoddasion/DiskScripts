# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 16:07:44 2019

@author: Wilhelm Hodder

Script to load and read images in a directory and turn them into animations
"""
#%% module imports

from PIL import Image
import matplotlib.pyplot as plt
from matplotlib import animation
from model_functions import filefriendly_time_from_index
#%% directory path formatting
stash_directory = 'Theta'
flight = 301
level_directory = 'Model_levels'
domain = 'subset1'
level = 1
path = f'../Figures/PNG/{301}/{level_directory}/{stash_directory}/{domain}/{level}/' # make sure to includ efinal forward slash

#%% figure name formatting
resolution = '1p5km'
stash_name = 'air_potential_temperature'
flight_date = 20180312
level_type = 'ML'
time = 'H21M30'

#%% initialise figures

plt.close()
def animate():
    