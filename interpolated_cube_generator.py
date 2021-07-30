#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 17:43:16 2020

@author: kse18nru
"""
#%%
import iris
import persinterpolate as pert
import sys
import numpy as np
import xarray as xr
import scipy.interpolate
#%% 
import os
os.environ["CUDA_VISIBLE_DEVICES"]="-1"
#%%

Case = {'res' : '0p5km', 'flight' : 301, 'varname' : 'm01s15i003', 'experiment' : 'Control', 
        'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'h'}
res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
print(Case.values())

resolutions = ['0p5km','1p5km','4p4km']
time_idx = np.arange(23,37)
print(time_idx)
target_heights = np.array([500,1000,1200])
#%% new code using iris
iris_code = False
if iris_code:
    #%% interpolate cubes onto selected altitudes
    if res == '4p4km':
        South = 136; North = 160; West = 150; East = 200
    if res == '1p5km':
        South = 155; North = 235; West = 90; East = 200
    if res == '0p5km':
        South = 190; North = 400; West = 50; East = 350
    time_idx = np.array([23,37])
    #, ('grid_longitude',p_3dcube.coords('grid_longitude')[0].points),
                     #('grid_latitude',p_3dcube.coords('grid_latitude')[0].points)]#[('altitude', [3000])]
    scheme = iris.analysis.Linear(extrapolation_mode='mask')
    cube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_{res}_um_{varname}_24hrs_{stream}_301.nc', varname)
    ## subset cube
    cube = cube[:,:,South:North, West:East] 
    xr_cube = xr.DataArray.from_iris(cube)
    print(xr_cube)
    
    #%%
    print(xr_cube.coords['altitude'].values)
    interp = scipy.interpolate.RegularGridInterpolator(
    tuple((xr_cube.coords['grid_longitude'].values, xr_cube.coords['grid_latitude'].values,
           xr_cube.coords['model_level_number'].values,xr_cube.coords['altitude'].values)),
    xr_cube.data[0][:],
    method='linear',
    bounds_error=np.nan)
    
    
    
    #%%
    cubes = []
    sample_points = [('altitude', [500,1000,1200])]
    xlen = len(cube.coords('grid_longitude')[0].points)
    ylen = len(cube.coords('grid_latitude')[0].points)
    
    #%%
    for y in range(ylen):
        rows = []
        for x in range(xlen):
            interp_point = cube[time_idx[0]:time_idx[-1],:,y,x].interpolate(sample_points, scheme)
            rows.append(interp_point)
            
            #print('coloumns concatenated at index',y, x)
        row = iris.cube.CubeList(rows).merge()
        cubes.append(row)
        row = 0
        print('completed y coloumn', y, 'out of ', ylen - 1)
    #%%
    interp_cube = iris.cube.CubeList(cubes).merge()
    cubes = []
    #print(cubes)
    print(interp_cube)
    #%%
    import iris.plot as iplt
    iplt.contourf(interp_cube[0][:,:,0,0])
    
    #%%
    save = True
    if save:
        try:
            iris.save([interp_cube], f'../../Model_data/{suite}/nc/Control/interp4d/{config}_{res}_{varname}_interp_southbay_flt{flight}_alt{target_heights[0]}to{target_heights[-1]}_tidx{time_idx[0]}to{time_idx[-1]}.nc')
            print('File successfully saved.')
        except:
            print('Directory save error')
            iris.save([interp_cube], f'../../Model_Data/Dump/{config}_{res}_{varname}_interp_southbay_flt{flight}_alt1000_tidx{time_idx[0]}to{time_idx[-1]}.nc')

#%% old code
domain = 'fulldomain'
legacy = True
if legacy:
    for res in resolutions:
           
        #%% load model data as iris cube - main variable - then concatenate
            appendix = ''
            datacube1 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}1{appendix}_flt{flight}.nc', varname); cubes = [datacube1]
            datacube2 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}2{appendix}_flt{flight}.nc', varname); cubes.append(datacube2)
            datacube3 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}3{appendix}_flt{flight}.nc', varname); cubes.append(datacube3)
            datacube4 = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{config}_{res}_um{stream}4{appendix}_flt{flight}.nc', varname); cubes.append(datacube4)
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
            #%%   
            subset = False
            if subset:
                if res == '4p4km':
                    South = 125; North = 160; West = 148; East = 200
                if res == '1p5km':
                    South = 150; North = 235; West = 85; East = 200
                if res == '0p5km':
                    South = 180; North = 400; West = 35; East = 350
                ## subset cube
                cube = cube[:,:,South:North, West:East]
            #%%
            print(cube)
            new_gridlat = iris.coords.DimCoord(cube.coords('grid_latitude')[0].points, standard_name='grid_latitude')
            new_gridlon = iris.coords.DimCoord(cube.coords('grid_longitude')[0].points, standard_name = 'grid_longitude')
            new_dimtime = iris.coords.DimCoord(cube.coords('time')[0].points[time_idx[0]:time_idx[-1]+1], standard_name = 'time')
            new_dimlevel = iris.coords.DimCoord(cube.coords('model_level_number')[0].points, standard_name = 'model_level_number')
            
            #%%
            # testcube = cube.copy()
            # testcube.remove_coord('grid_longitude')
            
            # print(testcube)
            # testcube.add_dim_coord(new_gridlon, 3)
            # print(testcube)
            #%%
            print(f'interpolation initiated: {res}')
            interp_data = pert.generate_interpolated_4Dfield(cube, target_heights, time_idx)
            
            #%% create iris cube
            
            newcube = iris.cube.Cube(interp_data, standard_name = cube.standard_name, long_name = cube.long_name,
                                     var_name = cube.var_name, units = cube.units, attributes = cube.attributes,
                                     )
            print('cube created')
            
            #%%
            newcube.add_dim_coord(new_gridlat, 2)
            newcube.add_dim_coord(new_gridlon, 3)
            newcube.add_dim_coord(new_dimtime, 0)
            #%%
            newcube.add_aux_coord(cube.coords('surface_altitude')[0], [2,3])
            altitude_auxcoord = iris.coords.AuxCoord(target_heights, standard_name = 'altitude', units = 'm')
            newcube.add_aux_coord(altitude_auxcoord, 1)
            print(newcube)
            #sys.exit()
            #%%
            save = True
            if save:
                
                
                try:
                    iris.save([newcube], f'../../Model_data/{suite}/nc/Control/interp4d/{config}_{res}_{varname}_interp_{domain}_flt{flight}_alt{target_heights[0]}to{target_heights[-1]}_tidx{time_idx[0]}to{time_idx[-1]}.nc')
                    print('File successfully saved.')
                except:
                    print('Directory save error')
                    iris.save([newcube], f'../../Model_Data/Dump/{config}_{res}_{varname}_interp_{domain}_flt{flight}_alt1000_tidx{time_idx[0]}to{time_idx[-1]}.nc')
                    
            #%%
            #print(newcube)