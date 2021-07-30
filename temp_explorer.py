"""
Created on Mon Aug 10 09:16:05 2020

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


#%%

file_contents = iris.load('../../Model_Data/u-bu807/nc/Control/RA1M_0p5km_umh1_flt301.nc')

print(file_contents)