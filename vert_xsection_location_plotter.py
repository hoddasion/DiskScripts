# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 22:06:36 2020

@author: kse18nru
"""
import matplotlib.pyplot as plt
from iris.analysis.cartography import rotate_pole
import numpy as np
import iris
import iris.coord_categorisation
import iris.plot as iplt
import pandas as pd
import xsection_workshop as workshop
import model_foundry as foundry
#from iris.maths import add
#%% execute function
legs = [3]
resolutions = ['0p5km', '1p5km', '4p4km']
D = '0p5km'; leg = 3
Case = {'res' : D, 'flight' : 301, 'varname' : 'upward_air_velocity', 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
fileleg = Case['fileleg']

maincube1 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg1.nc', varname)
maincube2 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg2.nc', varname)
maincube3 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg3.nc', varname)
maincube4 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg8.nc', varname)


orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{res}_surface_altitude_flt301.nc')[50:400,:300]
orog = orogcube.data
polelat = orogcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
polelon = orogcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
oroglat = orogcube.coord('grid_latitude').points
oroglon = orogcube.coord('grid_longitude').points
#%% extract cube coords
lon1 = maincube1.coord('grid_longitude').points
lat1 = maincube1.coord('grid_latitude').points
lon2 = maincube2.coord('grid_longitude').points
lat2 = maincube2.coord('grid_latitude').points
lon3 = maincube3.coord('grid_longitude').points
lat3 = maincube3.coord('grid_latitude').points
lon4 = maincube4.coord('grid_longitude').points
lat4 = maincube4.coord('grid_latitude').points

#%%
lons, lats = foundry.modf.unrotate_coords(oroglon, oroglat, polelon, polelat)
#%%

fig = plt.figure(figsize = (16,16))
#iplt.contourf(orogcube)
plt.contour(oroglon, oroglat,orog, colors = 'k')
plt.plot(lon1,lat1, color = 'r', linewidth = 10, label = 'Leg 1')
plt.plot(lon2,lat2, color = 'g', linewidth = 10, label = 'Leg 2')
plt.plot(lon3,lat3, color = 'y', linewidth = 10, label = 'Leg 3,4,5,6')
plt.plot(lon4,lat4, color = 'm', linewidth = 10, label = 'leg 8,9')
q1 = plt.contour(oroglon,oroglat, lons, colors = 'k', linestyle = '--')
q2 = plt.contour(oroglon,oroglat, lats, colors = 'k', linestyle = '--')
plt.clabel(q1);plt.clabel(q2)
plt.legend(fontsize = 20)
plt.xticks([]); plt.yticks([])
plt.title('Location of vertical cross-sections for flight legs 1,2,3, and 8', fontsize = 28)
plt.savefig('Dump/vert_xsec_location_legs123_ft301.png')
