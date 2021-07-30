# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 15:45:19 2020

@author: kse18nru
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import iris
import iris.coord_categorisation
import iris.plot as iplt
import model_foundry as foundry
import pandas as pd
from iris.analysis.cartography import rotate_pole
from iris.analysis import trajectory

#%%
Case = {'res' : '0p5km', 'flight' : 301, 'varname' : 'upward_heat_flux_in_air', 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'j'}
res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
print(Case.values())



fulldomain = True
if fulldomain:
    for res in ('0p5km', '1p5km', '4p4km'):
        orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}1_flt{flight}.nc', 'surface_altitude')
        
        landcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}1_flt{flight}.nc', 'land_binary_mask')
        print(landcube)
        #iplt.contourf(landcube)
        
        #plt.show()
        #%%
        #%%   
        
        cubes = [orogcube, landcube]
        iris.save(cubes,f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_fulldomain_flt301.nc')

#%%
interp_field = True
if interp_field:
    
    for res in ('0p5km', '1p5km', '4p4km'):
        orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}1_flt{flight}.nc', 'surface_altitude')
        
        landcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}1_flt{flight}.nc', 'land_binary_mask')
        print(landcube)
        #iplt.contourf(landcube)
        
        #plt.show()
        #%%
        #%%   
        if res == '4p4km':
            South = 125; North = 160; West = 148; East = 200
        if res == '1p5km':
            South = 150; North = 235; West = 85; East = 200
        if res == '0p5km':
            South = 180; North = 400; West = 35; East = 350
        sublandcube = landcube[South:North, West:East]
        subcube = orogcube[South:North, West:East]   
        iplt.contourf(subcube) 
        plt.show()
        print(subcube)
        cubes = [subcube, sublandcube]
        iris.save(cubes,f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_southbay_flt301.nc')
#%%
snaesfellness = False
if snaesfellness:
    
    for res in ('0p5km', '1p5km', '4p4km'):
        orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}1_flt{flight}.nc', 'surface_altitude')
        print(orogcube)
        
        iplt.contourf(orogcube)
        plt.show()
        #%%
        if res == '4p4km':
            South = 143; North = 155; West = 150; East = 200
        if res == '1p5km':
            South = 190; North = 220; West = 90; East = 200
        if res == '0p5km':
            South = 280; North = 370; West = 50; East = 350
            
        subcube = orogcube[South:North, West:East]   
        iplt.contourf(subcube) 
        plt.show()
        print(subcube)
        
        #iris.save(subcube,f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
        
westcoast = False
if westcoast:
    
    for res in ('0p5km', '1p5km', '4p4km'):
        orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}1_flt{flight}.nc', 'surface_altitude')
        print(orogcube)
        
        iplt.contourf(orogcube)
        plt.show()
        #%%
        if res == '4p4km':
            South = 100; North = 220; West = 140; East = 210
        if res == '1p5km':
            South = 100; North = 360; West = 70; East = 240
        if res == '0p5km':
            South = 0; North = -1; West = 0; East = -1
            
        subcube = orogcube[South:North, West:East]   
        iplt.contourf(subcube) 
        plt.show()
        print(subcube)
        
        #iris.save(subcube,f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_westcoast_flt301.nc')