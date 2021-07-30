# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 19:03:49 2020

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
from iris.analysis.cartography import rotate_pole

#%%

Case = {'res' : '0p5km', 'flight' : 301, 'varname' : 'upward_air_velocity', 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'h'}
res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
print(Case.values())

#%% load observational data

## load quality controlled and interval-meaned observational data from Database
# using pandas
database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_60s_301.txt', delimiter = ' ')
db_legno = database['legno']; legno = 8
legcond = db_legno == legno
db_w = np.array(database['w'][legcond])
db_lon = np.array(database['lon'][legcond])
db_lat = np.array(database['lat'][legcond])
db_alt = np.array(database['altgps'][legcond])

db_mtime = np.array(database['meantime'][legcond])


#%% load model data as iris cube - main variable - then concatenate

datacube1 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}1_flt{flight}.nc', varname); cubes = [datacube1]
datacube2 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}2_flt{flight}.nc', varname); cubes.append(datacube2)
datacube3 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}3_flt{flight}.nc', varname); cubes.append(datacube3)
datacube4 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}4_flt{flight}.nc', varname); cubes.append(datacube4)
cube = iris.cube.CubeList(cubes).concatenate()[0] # a complete cube spanning the entire 24 hour forecast

## derive altitude coordinate seperately
try:
    delta = cube.coord('atmosphere_hybrid_height_coordinate') # load hybrid height
    sigma = cube.coord('sigma') # load sigma
    orography = cube.coord('surface_altitude') # load orography
    factory = iris.aux_factory.HybridHeightFactory(delta, sigma, orography)
    factory.rename('altitude') # rename derived coordinate from 'unknown' to 'altitude'
    cube.add_aux_factory(factory) # integrate into cube
    ## save new cube to file
    #print(cube)
except:
    print('Exception in altitude derivation, main var')
    sys.exit()

#%% load model theta data cubes, then concatenate
datacube1 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}1_flt{flight}.nc', varname); cubes = [datacube1]
datacube2 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}2_flt{flight}.nc', varname); cubes.append(datacube2)
datacube3 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}3_flt{flight}.nc', varname); cubes.append(datacube3)
datacube4 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}4_flt{flight}.nc', varname); cubes.append(datacube4)
thetacube = iris.cube.CubeList(cubes).concatenate()[0] # a complete cube spanning the entire 24 hour forecast

## derive altitude coordinate seperately
try:
    delta = thetacube.coord('atmosphere_hybrid_height_coordinate') # load hybrid height
    sigma = thetacube.coord('sigma') # load sigma
    orography = thetacube.coord('surface_altitude') # load orography
    factory = iris.aux_factory.HybridHeightFactory(delta, sigma, orography)
    factory.rename('altitude') # rename derived coordinate from 'unknown' to 'altitude'
    thetacube.add_aux_factory(factory) # integrate into cube
    ## save new cube to file
    print(thetacube)
except:
    print('Exception in altitude derivation, theta stage')
    sys.exit()
## dump excess memory
datacube1 = 0; datacube2 = 0; datacube3 = 0; datacube4 = 0
#%% subsample hourly/ every second time output starting at 1am
cube = cube[1::2]
thetacube = thetacube[1::2]
print(cube.coord('time'))

#%% get rotated pole coords
polelat = cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
polelon = cube.coord('grid_longitude').coord_system.grid_north_pole_longitude

#%% rotate obs coordinates onto model grid
rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
rot_db_lon = rot_db_lon + 360

#%% identify x-section vertices 
xstart = rot_db_lon[0]
ystart = rot_db_lat[0]
xend = rot_db_lon[-1]
yend = rot_db_lat[-1]

#%% generate coordinate sample line
if res == '4p4km': gap = 0.04
if res == '1p5km': gap = 0.015
if res == '0p5km': gap = 0.005
num =int(( np.absolute((xend**2 + yend**2)**0.5 - (xstart**2 + ystart**2)**0.5 )/gap)//1)
print(num)
xrow = np.linspace(xstart,xend, num = num)
yrow = np.linspace(ystart,yend, num = num)

sample_points = [('grid_latitude', yrow),('grid_longitude', xrow)]

#%% interpolate cubes
cubeslice = cube.interpolate(sample_points, iris.analysis.Nearest())
print(cubeslice)
#%%
print(cubeslice[0,:, 0].data)
iplt.contourf(cubeslice[0], coords = ['grid_longitude','model_level_number'])
iplt.show()



