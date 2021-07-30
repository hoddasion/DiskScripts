# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 15:37:47 2019

@author: kse18nru
"""

#import model_functions as modf
import iris

#%% global configurations 

import warnings
warnings.filterwarnings("ignore") ##BAD! If something doesn't work, re-enable warnings to check
#%% load first variable to analyse: Theta on theta levels
#file_name = 'RA1M_0p5km_umh1_up1608_flt301.nc'
file_name = 'RA1M_0p5km_umi1_flt301.nc'
path1 = '../../Model_Data/u-bu807/nc/Control/'
stash = 'm01s15i002'
#modf.file_info(f'{path1}{file_name}')
#modf.cube_info(f'{path1}{file_name}', stash)

cubes = iris.load(f'{path1}{file_name}')
print(file_name)
print(cubes)