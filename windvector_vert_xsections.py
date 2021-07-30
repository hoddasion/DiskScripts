# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 09:05:12 2020

@author: kse18nru
"""

#%% Import modules and scripts


import matplotlib.pyplot as plt
import numpy as np
import iris
import iris.coord_categorisation

import xsection_workshop as workshop
from iris.analysis.maths import add
#%% execute function
legs = [1,2,3]
resolutions = ['0p5km', '1p5km', '4p4km']

#%%
D = '0p5km'; leg = 1
CaseWind = {'res' : D, 'flight' : 301, 'varname_w' : 'upward_air_velocity', 
            'varname_u' : 'm01s15i002', 'varname_v' : 'm01s15i003',
            'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 
            'filestream' : 'i', 'fileleg' : leg}
res = CaseWind['res']; suite = CaseWind['suite']; exp = CaseWind['experiment']; 
varname_w = CaseWind['varname_w']; varname_u = CaseWind['varname_u']; varname_v = CaseWind['varname_v']
config = CaseWind['config']; stream = CaseWind['filestream']; flight = CaseWind['flight']
fileleg = CaseWind['fileleg']

wcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname_w}/{config}_{res}_{varname_w}_flt{flight}_leg{fileleg}.nc', varname_w)
ucube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_Slices/{varname_u}/{config}_{res}_{varname_u}_flt{flight}_leg{fileleg}.nc', varname_u)
vcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_Slices/{varname_v}/{config}_{res}_{varname_v}_flt{flight}_leg{fileleg}.nc', varname_v)
thetacube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname_w}/{config}_{res}_{varname_w}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')

print(wcube,'\n', ucube,'\n', vcube)
#print(ucube.data)
#%%
coord = 'altitude'
print(np.shape(ucube.coord(coord).points))
print(wcube)
print(ucube.coord(coord).points[0,:5])
#%% updating v and u cube units
units = wcube.units
wcube.units = '1'
#wcube.convert_units('1')# = units, ucube.units = units
print(wcube.units)
#%%


#%% processing magnitude variable
hozmagcube = ( vcube**2 + ucube**2)**0.5
magcube = (add(hozmagcube**2,wcube**2))**0.5
print(magcube)
#%%
#fig = workshop.plot_vertical_xsection(maincube, thetacube, orogcube,index = 6, Case = Case, cmap = 'seismic', vmin = -10,vmax = 10, levels = 20, alt_limit = 2000)
#plt.show()