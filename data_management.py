#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 15:22:57 2020

@author: kse18nru
"""

import iris
import pandas as pd




def load_model_wsp(res, varname_u, varname_v, location,fileleg,flight, exp, suite = 'u-bu807', config = 'RA1M', level_range = -1):
    ucube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{location}/{varname_u}/{config}_{res}_{varname_u}_flt{flight}_leg{fileleg}.nc', varname_u)[:,:level_range]
    vcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{location}/{varname_v}/{config}_{res}_{varname_v}_flt{flight}_leg{fileleg}.nc', varname_v)[:,:level_range]
    maincube = (ucube**2 + vcube**2)**0.5 # get wind speed mangitude by quadriture
    #maincube.units = varname_u.units
    return maincube

def load_model_w(res, location,fileleg,flight, exp, varname = 'upward_air_velocity', suite = 'u-bu807', config = 'RA1M', level_range = -1):
    maincube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', varname)[:,:level_range]
    return maincube

def load_model_specific_hum(res, location,fileleg,flight, exp, varname = 'specific_humidity', suite = 'u-bu807', config = 'RA1M', level_range = -1):
    maincube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', varname)[:,:level_range]
    maincube = maincube * 1000
    maincube.units = 'g/kg'
    print(maincube.units)
    return maincube

def load_model_theta(res, location,fileleg,flight, exp, varname = 'air_potential_temperature', suite = 'u-bu807', config = 'RA1M', level_range = -1, secondary_var = 'air_potential_temperature'):
    maincube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{secondary_var}/{config}_{res}_{secondary_var}_flt{flight}_leg{fileleg}.nc', varname)[:,:level_range]
    return maincube

def load_model_windstress(res, varname_nws, varname_ews, location, fileleg, flight, exp, suite = 'u-bu807', config = 'RA1M', level_range = -1):
    ecube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{location}/{varname_ews}/{config}_{res}_{varname_ews}_flt{flight}_leg{fileleg}.nc', varname_ews)[:,:level_range]
    ncube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/{location}/{varname_nws}/{config}_{res}_{varname_nws}_flt{flight}_leg{fileleg}.nc', varname_nws)[:,:level_range]
    #print(ecube,'\n', ncube)
    newdata = (ecube.data**2 + ncube.data**2)**0.5 
    maincube = ecube
    maincube.shortname = 'downward_magnitude_stress'
    maincube.longname = 'atmosphere_downward_magnitude_stress'
    maincube.data = newdata
    # regrid one cube onto the grid of the other
    #ncube.regrid(ecube, iris.analysis.Linear())
    #maincube = (ecube**2 + ncube**2)**0.5 # get wind speed mangitude by quadriture
    #maincube.units = varname_ews.units
    return maincube

def load_model_pressure(res, location,fileleg,flight, exp, varname = 'air_pressure', suite = 'u-bu807', config = 'RA1M', level_range = -1):
    maincube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', varname)[:,:level_range]
    return maincube

def load_model_vapor_flux(res, location,fileleg,flight, exp, varname = 'upward_water_vapor_flux_in_air', suite = 'u-bu807', config = 'RA1M', level_range = -1):
    maincube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', varname)[:,:level_range]
    return maincube

def load_model_heat_flux(res, location,fileleg,flight, exp, varname = 'upward_heat_flux_in_air', suite = 'u-bu807', config = 'RA1M', level_range = -1):
    maincube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', varname)[:,:level_range]
    return maincube

def load_model_tke(res, location,fileleg,flight, exp, varname = 'm01s03i473', suite = 'u-bu807', config = 'RA1M', level_range = -1):
    maincube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', varname)[:,:level_range]
    return maincube

def calc_lhf(q,theta,p,vapor, q_conversion = 1e-3, R = 287, L = 2.25e6):
    print('Commencing with mixing ratio')
    mixing_ratio = q*q_conversion/(1-q*q_conversion)
    print('Calculated mixing ratio... commencing with virtual theta')
    virtual_temp = (1+0.61*mixing_ratio)*theta
    print('Calculated virtual theta... commencing with density')
    density = p/(R*virtual_temp)
    print('Calculated density... commencing with lh')
    LH_data =vapor*L# density*2.5*10e6*vapor
    print('Calculated lh... returning data')
    return LH_data

def unpack_testmeta(df):
    
    return 0
    