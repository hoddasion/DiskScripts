# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 13:40:22 2021

@author: kse18nru
"""

#%% module imports
import iris
from iris.analysis.cartography import rotate_pole
import data_management as dmna
import pandas as pd
import numpy as np
import persinterpolate as pert
import matplotlib
import matplotlib.pyplot as plt
import xsection_workshop as workshop
import iris.coord_categorisation
import iris.plot as iplt
import matplotlib.gridspec as gridspec
from iris.analysis import trajectory
import cf_units as unit
import datetime 

#%%

def lapse_rate(theta_rate_df, p_rate_df, p = 100000, T = 273.15, cp = 1004, R = 287, p0 = 100000, thtname = 'THETARATE', presname = 'PRESRATE', Tsurfcube = None, Tavr = False):
    """
    A function to calculate the adiabatic lapse rate for air temperature from measured values for 
    the lapse rate of theta and pressure. Derived from the definition of potential temperature.

    Parameters
    ----------
    theta_rate : Measured lapse rate of potential temperature, in K/m.
    p_rate : Measured rate of change of pressure with altitude z.
    p : Instantaneous pressure measured at altitude z, in Pa/m.
    T : Instananeous measured air temerature at z, in K.
    cp : Specific heat at constant pressure. The default is 1004 J/K/kg.
    R : Gas constant for dry air. The default is 287 J/K/kg.
    p0 : Standard/reference pressure. The default is 100000 Pa.

    Returns
    -------
    Air temperature measured adiabatic lapse rate in K/m.

    """
    Rcp = R/cp
    p = np.array(p_rate_df.SURFPRES)
    if Tsurfcube != None:
        if Tavr == True:
            cubedata = Tsurfcube.data
            meandata = np.nanmean(cubedata, axis = 1 )
            
    first_term = ((p/p0)**Rcp)*np.array(theta_rate_df[thtname]) # measured potential temperature lapse rate term
    second_term = (T/p)*Rcp*np.array(p_rate_df[presname]) # hydrostatic term
    lapserate = first_term + second_term
    
    #print(first_term, second_term)
    ## set up up new dataframe for results
    new_frame = pd.DataFrame(np.array(theta_rate_df), columns = ['TIME','LAPSERATE','GRIDY','GRIDX'])
    
    new_frame.LAPSERATE = lapserate
    return new_frame

def measure_theta_rate(theta_cube, sample_coords,max_alt, coloumn_names = ['TIMESTAMP', 'THETARATE','GRIDY', 'GRIDX']):
    """
    

    Parameters
    ----------
    theta_cube : 3D diagnostic cube of air potential temperature.
    sample_coords : list of sample points set up for iris trajectory function.
    max_alt : Maximum altitude to measure up to, in m.

    Returns
    -------
    Dataframe, containing time, mean lapse rate, grid latitude and grid longitude.

    """
    ## trajecotry interpolate cube to get coloumns
    coloumns = trajectory.interpolate(theta_cube,sample_coords, method =  'nearest') 
    ## identify closest altitude index
    df_contents = []
    for c in range(len(coloumns.data[0,0,:])): # loop through coloumns
        for t in range(len(coloumns.data[:,0,0])):
            altitudes = coloumns[t,:,c].coords('altitude')[0].points
            alt_idx = len(altitudes[altitudes < max_alt]) - 1
            altitudes = altitudes[altitudes < max_alt]
            coloumn_data = coloumns[t,:alt_idx,c].data
            
            rates = []
            for h in range(alt_idx - 1):
                val_diff = coloumn_data[h + 1] - coloumn_data[h]
                z_diff = altitudes[h+1] - altitudes[h]
                rates.append(val_diff/z_diff)
            mean_rate = np.nanmean(np.array(rates))
                
            df_rows = [coloumns.coords('time')[0].points[t],mean_rate, 
                       coloumns.coords('grid_latitude')[0].points[c],coloumns.coords('grid_longitude')[0].points[c]]
            df_contents.append(df_rows)
        
    df = pd.DataFrame(df_contents, columns = coloumn_names)
    df[coloumn_names[0]] = pd.to_datetime(df[coloumn_names[0]], unit = 'h')
    return df
    
def measure_pressure_rate(p_cube, sample_coords,max_alt, coloumn_names = ['TIMESTAMP', 'PRESRATE','GRIDY', 'GRIDX', 'SURFPRES', 'LOWESTZ']):
    """
    
      
    Parameters
    ----------
    theta_cube : 3D diagnostic cube of air potential temperature.
    sample_coords : list of sample points set up for iris trajectory function.
    max_alt : Maximum altitude to measure up to, in m.
      
    Returns
    -------
    Dataframe, containing time, mean pressure rate of change, grid latitude and grid longitude.
      
    """
    ## trajecotry interpolate cube to get coloumns
    coloumns = trajectory.interpolate(p_cube,sample_coords, method =  'nearest') 
    ## identify closest altitude index
    df_contents = []
    for c in range(len(coloumns.data[0,0,:])): # loop through coloumns
        for t in range(len(coloumns.data[:,0,0])):
            altitudes = coloumns[t,:,c].coords('altitude')[0].points
            alt_idx = len(altitudes[altitudes < max_alt]) - 1
            altitudes = altitudes[altitudes < max_alt]
            coloumn_data = coloumns[t,:alt_idx,c].data
            
            rates = []
            for h in range(alt_idx - 1):
                val_diff = coloumn_data[h + 1] - coloumn_data[h]
                z_diff = altitudes[h+1] - altitudes[h]
                rates.append(val_diff/z_diff)
            mean_rate = np.nanmean(np.array(rates))
                
            df_rows = [coloumns.coords('time')[0].points[t],mean_rate, 
                       coloumns.coords('grid_latitude')[0].points[c],coloumns.coords('grid_longitude')[0].points[c],
                       coloumn_data[0], altitudes[0]]
            df_contents.append(df_rows)
        
    df = pd.DataFrame(df_contents, columns = coloumn_names)
    df[coloumn_names[0]] = pd.to_datetime(df[coloumn_names[0]], unit = 'h')
    return df