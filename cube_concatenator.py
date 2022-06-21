#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 14:44:28 2020

@author: Wilhelm Hodder
For purposes of concatenating forecast data files into single 24hour forecast files
"""
import iris
import sys
## defining resolution domain, directory, and filenames
variable = 'toa_outgoing_shortwave_flux'
name = variable
suite = 'u-cc134'
nc_path = f'D:/Project/Model_Data/{suite}/'
experiment = 'RA1M'
flight = 306
stream = 'pg'
check_file = True
if check_file:
    filename = f'{experiment}_0p5km_um'
    cubes = iris.load(f'{nc_path}{filename}{stream}1_flt{flight}.nc')
    print(cubes)
## load model data cubes
mlevel = True
if mlevel:
    for res in ['0p5km']:#, '1p5km','4p4km']:
        filename = f'{experiment}_{res}_um'
        cubes = []
        cube1 = iris.load_cube(f'{nc_path}{filename}{stream}1_flt{flight}.nc', name); cubes.append(cube1) # fcst hours 0030-0600
        cube2 = iris.load_cube(f'{nc_path}{filename}{stream}2_flt{flight}.nc', name); cubes.append(cube2) # fcst hours 0630-1200
        cube3 = iris.load_cube(f'{nc_path}{filename}{stream}3_flt{flight}.nc', name); cubes.append(cube3) # fcst hours 1230-1800
        cube4 = iris.load_cube(f'{nc_path}{filename}{stream}4_flt{flight}.nc', name); cubes.append(cube4) # fcst hours 1830-2400
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
            iris.save(cube, f'{nc_path}{filename}_{variable}_24hrs_{stream}_{flight}.nc')
        except:
            print('Exception occured')
            iris.save(cube, f'{nc_path}{filename}_{variable}_24hrs_{stream}_{flight}.nc')
        
        
        
        print(f'File complete: {res} {variable}')
        
#%%        
plevel = False
if plevel:
    stream = 'j'
    for res in ['0p5km','1p5km','4p4km']:
        
        varnames = ['air_temperature', 'specific_humidity', 'upward_air_velocity', 'x_wind', 'y_wind', 'geopotential_height']
        nc_path = 'D:/Project/Model_Data/u-cc134/'
        #for varname in varnames:
        filename = f'RA1M_{res}_ump{stream}1_flt{flight}.nc'
        filesuffix = f'RA1M_{res}_ump'
        filenames = [f'{nc_path}{filesuffix}{stream}1_flt{flight}.nc',f'{nc_path}{filesuffix}{stream}2_flt{flight}.nc',
                     f'{nc_path}{filesuffix}{stream}3_flt{flight}.nc',f'{nc_path}{filesuffix}{stream}4_flt{flight}.nc']
        cubes = iris.load(filenames, varnames).concatenate()
        lm_cube = iris.load_cube(f'{nc_path}{filename}', 'land_binary_mask')
        alt_cube = iris.load_cube(f'{nc_path}{filename}', 'surface_altitude')
        cubes.append(lm_cube)
        cubes.append(alt_cube)
        print(cubes, '\n\n\n')
        
        iris.save(cubes, f'{nc_path}RA1M_{res}_cc134_24hrs_pressurevars_{flight}.nc')
           

glm_select = False
if glm_select:
    cube1 = iris.load_cube()
        