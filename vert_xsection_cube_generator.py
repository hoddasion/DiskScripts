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
from iris.analysis import trajectory

#%%

Case = {'res' : '0p5km', 'flight' : 306, 'varname' : 'y_wind', 'experiment' : 'Control', 
        'config' : 'RA1M', 'suite' : 'u-cc134', 'filestream' : 'ph'}
res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
print(Case.values())

#%%
leg_xsections = True
if leg_xsections:
    legs = [6]#[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
    res = ['0p5km','1p5km','4p4km']
    for res in res:
        #%%
                
        cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_{res}_um_{varname}_24hrs_{stream}_{flight}.nc', varname)
        #print(cube)
        for i in legs:
                print('Cycle', f'leg {i+1}', res)
                #%% load observational data
                
                ## load quality controlled and interval-meaned observational data from Database
                # using pandas
                database = pd.read_csv('D:/Project/Obvs_Data/Databases/IGP_flights_database_obs_60s_306.txt', delimiter = ' ')
                db_legno = database['legno']; legno = i + 1
                legcond = db_legno == legno
                db_w = np.array(database['w'][legcond])
                db_lon = np.array(database['lon'][legcond])
                db_lat = np.array(database['lat'][legcond])
                db_alt = np.array(database['altgps'][legcond])
                
                db_mtime = np.array(database['meantime'][legcond])

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
                if i == 0:
                    xstart -= 0.05
                    xend -= 0.05
                grad = (yend - ystart)/(xend - xstart)
                #%% generate coordinate sample line
                if res == '4p4km': gap = 0.04
                if res == '1p5km': gap = 0.015
                if res == '0p5km': gap = 0.005
                delta_factor = 1.5
                dx = delta_factor*(xend - xstart)
                dy = delta_factor*(yend - ystart)
                dxs =4.5*(xend - xstart)
                dys = 4.5*(yend - ystart)
                num =int(( np.absolute(((xend+dx)**2 + (yend+dy)**2)**0.5 - ((xstart-dxs)**2 + (ystart-dys)**2)**0.5 )/gap)//1)
                print(num)
                print(xstart, xend, ystart, yend)
                #%%
                
                xrow = np.linspace(xstart-dxs,xend+dx, num = num)
                yrow = np.linspace(ystart-dys,yend+dy, num = num)
                
                sample_points = [('grid_latitude', yrow),('grid_longitude', xrow)]
                
                #%% interpolate cubes
                NS_offsets = True; NSofst = 0.20; NSofst_str = '0p20'
                try:
                    cubeslice = trajectory.interpolate(cube,sample_points, method =  'nearest')
                    if NS_offsets:
                        sample_points_N = [('grid_latitude', yrow + NSofst),('grid_longitude', xrow)]
                        cubeslice_N = trajectory.interpolate(cube,sample_points_N, method =  'nearest')
                        sample_points_S = [('grid_latitude', yrow - NSofst),('grid_longitude', xrow)]
                        cubeslice_S = trajectory.interpolate(cube,sample_points_S, method =  'nearest')
                except:
                    print(f'Exception occured for leg {i+1}')
                    pass
                
                #print(cubeslice)
                #%%
                save = True
                if save:
                    if num > 10:
                        iris.save(cubeslice, f'D:/Project/Model_data/{suite}/Vertical/{varname}/{config}_vertslice_{res}_{varname}_flt{flight}_leg{legno}long.nc')
                        if NS_offsets:
                            iris.save(cubeslice_N, f'D:/Project/Model_data/{suite}/Vertical/{varname}/{config}_vertslice_{res}_{varname}_flt{flight}_leg{legno}long_ofts_N{NSofst_str}.nc')
                            iris.save(cubeslice_S, f'D:/Project/Model_data/{suite}/Vertical/{varname}/{config}_vertslice_{res}_{varname}_flt{flight}_leg{legno}long_ofts_S{NSofst_str}.nc')
                
                #%%
                ocube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{config}_0p5km_um_land_binary_mask_24hrs_pi_{flight}.nc', 'land_binary_mask')
                #print(ocube)
                odata = ocube.data; olat = ocube.coord('grid_latitude').points; olon = ocube.coord('grid_longitude').points
                slice_lon = cubeslice.coord('grid_longitude').points; slice_lat = cubeslice.coord('grid_latitude').points
                
                fig = plt.figure(figsize = (16,16))
                plt.contour(olon,olat,odata)
                plt.plot(slice_lon,slice_lat)
                plt.plot(slice_lon,slice_lat+NSofst)
                plt.plot(slice_lon,slice_lat-NSofst)
                #plt.close()
                plt.show()



#%%
leg_xsections_with_concat = False
if leg_xsections_with_concat:
    legs = [0,1,2,3,4,5,6,7,8,9,10]
    res = ['0p5km','1p5km','4p4km']
    for res in res:
        for i in legs:
                print('Cycle', f'leg {i+1}', res)
            #try:
                #%% load observational data
                
                ## load quality controlled and interval-meaned observational data from Database
                # using pandas
                database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_60s_301.txt', delimiter = ' ')
                db_legno = database['legno']; legno = i + 1
                legcond = db_legno == legno
                db_w = np.array(database['w'][legcond])
                db_lon = np.array(database['lon'][legcond])
                db_lat = np.array(database['lat'][legcond])
                db_alt = np.array(database['altgps'][legcond])
                
                db_mtime = np.array(database['meantime'][legcond])
                
                
                #%% load model data as iris cube - main variable - then concatenate
                appendix = ''
                datacube1 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}1{appendix}_flt{flight}.nc', varname)[:,:35]; cubes = [datacube1]
                datacube2 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}2{appendix}_flt{flight}.nc', varname)[:,:35]; cubes.append(datacube2)
                datacube3 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}3{appendix}_flt{flight}.nc', varname)[:,:35]; cubes.append(datacube3)
                datacube4 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}4{appendix}_flt{flight}.nc', varname)[:,:35]; cubes.append(datacube4)
                cube = iris.cube.CubeList(cubes).concatenate()[0] # a complete cube spanning the entire 24 hour forecast
                print(cube)
                #print(cube.coord('level_height').points)
                ## derive altitude coordinate seperately
                
                
                #%% load model theta data cubes, then concatenate
                tvarname = 'air_potential_temperature'
                datacube1 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_umi1_flt{flight}.nc', tvarname)[:,:35]; cubes = [datacube1]
                datacube2 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_umi2_flt{flight}.nc', tvarname)[:,:35]; cubes.append(datacube2)
                datacube3 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_umi3_flt{flight}.nc', tvarname)[:,:35]; cubes.append(datacube3)
                datacube4 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_umi4_flt{flight}.nc', tvarname)[:,:35]; cubes.append(datacube4)
                thetacube = iris.cube.CubeList(cubes).concatenate()[0] # a complete cube spanning the entire 24 hour forecast
                print(thetacube)
                ## derive altitude coordinate seperately
                #sys.exit()
                try:
                    delta = thetacube.coord('atmosphere_hybrid_height_coordinate') # load hybrid height
                    sigma = thetacube.coord('sigma') # load sigma
                    orography = thetacube.coord('surface_altitude') # load orography
                    factory = iris.aux_factory.HybridHeightFactory(delta, sigma, orography)
                    factory.rename('altitude') # rename derived coordinate from 'unknown' to 'altitude'
                    thetacube.add_aux_factory(factory) # integrate into cube
                    ## save new cube to file
                    #print(thetacube)
                except:
                    print('Exception in altitude derivation, theta stage')
                    sys.exit()
                ## dump excess memory
                datacube1 = 0; datacube2 = 0; datacube3 = 0; datacube4 = 0
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
                    #try:
                        new_coord = iris.coords.AuxCoord(thetacube.coord('surface_altitude').points, standard_name = 'surface_altitude', units = 'm')
                        cube.add_aux_coord(new_coord, np.shape(thetacube.coord('surface_altitude').points))
                        delta = cube.coord('level_height') # load hybrid height
                        sigma = cube.coord('sigma') # load sigma
                        orography = cube.coord('surface_altitude') # load orography
                        factory = iris.aux_factory.HybridHeightFactory(delta, sigma, orography)
                        factory.rename('altitude') # rename derived coordinate from 'unknown' to 'altitude'
                        cube.add_aux_factory(factory) # integrate into cube
                    #except:
                        #print('Exception in altitude derivation, main var:\n', sys.exc_info()[0])
                        #sys.exit()
                #%% regridding
                regridding = True
                if regridding:
                    print('regridding intitiated')
                    cube.regrid(thetacube, iris.analysis.Linear())
                    print('cycle regridding complete')
                print(cube)
                #%% subsample hourly/ every second time output starting at 1am
                #cube = cube[1::2]
                #thetacube = thetacube[1::2]
               # print(cube.coord('time'))
                
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
                
                grad = (yend - ystart)/(xend - xstart)
                #%% generate coordinate sample line
                if res == '4p4km': gap = 0.04
                if res == '1p5km': gap = 0.015
                if res == '0p5km': gap = 0.005
                delta_factor = 1.0
                dx = delta_factor*(xend - xstart)
                dy = delta_factor*(yend - ystart)
                num =int(( np.absolute(((xend+dx)**2 + (yend+dy)**2)**0.5 - ((xstart-dx)**2 + (ystart-dy)**2)**0.5 )/gap)//1)
                print(num)
                print(xstart, xend, ystart, yend)
                #%%
                
                xrow = np.linspace(xstart-dx,xend+dx, num = num)
                yrow = np.linspace(ystart-dy,yend+dy, num = num)
                
                sample_points = [('grid_latitude', yrow),('grid_longitude', xrow)]
                
                #%% interpolate cubes
                cubeslice = trajectory.interpolate(cube,sample_points, method =  'nearest')
                thetaslice = trajectory.interpolate(thetacube, sample_points, method = 'nearest')
                cubes = [cubeslice, thetaslice]
                #print(cubeslice)
                #%%
                save = True
                if save:
                    try:
                        iris.save(cubes, f'../../Model_data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{legno}.nc')
                    except:
                        print('Directory save error; cycle', i+1)
                        iris.save(cubes, f'Dump/{config}_{res}_{varname}_flt{flight}_leg{legno}.nc')
                #%%
                #print(cubeslice[0,:, 0].data)
                #iplt.contourf(cubeslice[0],coords = ['grid_longitude','altitude'])
                #plt.plot(cubslice[0])
                #iplt.show()
                #fig, ax = plt.subplots(1,1, figsize = (12,12))
                #q0 = ax.contourf(cubeslice.data[0,0])
                #%%
                ocube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_0p5km_umi1_flt{flight}.nc', 'land_binary_mask')
                print(ocube)
                odata = ocube.data; olat = ocube.coord('grid_latitude').points; olon = ocube.coord('grid_longitude').points
                slice_lon = cubeslice.coord('grid_longitude').points; slice_lat = cubeslice.coord('grid_latitude').points
                
                fig = plt.figure(figsize = (16,16))
                plt.contour(olon,olat,odata)
                plt.plot(slice_lon,slice_lat)
                plt.show()
                #%%
                #print(cubeslice.coord('grid_longitude'))
            #except:
            #    print(f'Cycle {i+1} failed. ')
             #   continue
    
#%%
tpoint_xsections = False
if tpoint_xsections:
    res = ['0p5km','1p5km','4p4km']
    for res in res:
        tpoint_db = pd.read_csv('D:/Project/Figures/PNG/301/u-bu807/SlowBubble/tpoint_coors_301.csv')
        print(tpoint_db)
        
        print('Owema mo shindaru')
