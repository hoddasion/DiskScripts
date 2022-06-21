# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 15:28:15 2021

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
import glm_functions as glmf

#%%
def calc_bouancy_freq(theta, lapse_rate , g = 9.81):
    """
    Parameters
    ----------
    theta : float or array, potential temeprature values.
    lapse_rate : Float or array, strictly dtheta/dz. The default is 9.8*(10**(-3)) K/km, the dry adiabatic lapse rate.
    g : Gravitational acceleration. The default is 9.81 ms^(-2).

    Returns
    -------
    Float or array
        The Brunt-Vaeisaelae frequency.

    """
    theta = theta# + 0j
    lapse_rate = lapse_rate# + 0j
    return np.sqrt((g/theta)*(lapse_rate))

def calc_N_error(N, theta,theta_e, lapse_rate, lapse_e):
    
    return 0.5*N*np.sqrt((theta_e/theta)**2+(lapse_e/lapse_rate)**2)
def calc_Fr_error(Fr, N, N_e, U, U_e, h = 1000, h_e = 0):
    
    return Fr*np.sqrt( (N_e/N)**2 + (U_e/U)**2 + (h_e/h)**2 )

#%%
obs_times = [('20180319T0600Z',11),('20180319T1200Z',23),('20180319T1800Z',35)]
time_idcs = [11,23,35]

#%% load model data
suite = 'u-cc134'
res = '1p5km'
um_path = f'D:/Project/Model_Data/{suite}/'
## velocity diagnostics
nwind_cube = iris.load_cube(f'{um_path}RA1M_{res}_um_y_wind_24hrs_ph_306.nc', 'y_wind')[:,:36]
ewind_cube = iris.load_cube(f'{um_path}RA1M_{res}_um_x_wind_24hrs_ph_306.nc', 'x_wind')[:,:36]

magwind_cube = iris.load_cube(f'{um_path}RA1M_{res}_um_x_wind_24hrs_ph_306.nc', 'x_wind')[:,:36]
mag_data = (nwind_cube.data[:,:,:-1,:]**2 + ewind_cube.data**2)**0.5 
magwind_cube.data = mag_data
## theta diagnostic
theta_cube = iris.load_cube(f'{um_path}RA1M_{res}_um_air_potential_temperature_24hrs_pi_306.nc', 'air_potential_temperature')[:,:36]
#%% oad obs data
for info in obs_times:
    
    obs_time = info[0]
    print(obs_time)
    #%%
    obs_path = 'D:/Project/Obvs_Data/profiles/keflavik_04018BIKF_profiles_upto500mb.csv'
    obs_df = pd.read_csv(obs_path)
    obs_profiles = obs_df[obs_df.Time == obs_time]
    idx_1km = 22
    obs_h = np.array(obs_profiles.HGHT).astype(np.float)
    
    #%%
    sample_cond = (obs_h < 2000) *  (obs_h > 200)
    U_sample_cond = (obs_h < 1500)*(obs_h > 500)
    obs_h = obs_h[sample_cond]
    
    obs_p = np.array(obs_profiles.PRES).astype(np.float)[sample_cond]
    obs_theta = np.array(obs_profiles.THTA).astype(np.float)[sample_cond]
   
    obs_speed =np.array(obs_profiles.SKNT).astype(np.float)[U_sample_cond]*0.5144 # convert from knots to m/s
    
    
    
    
    
    #%% isolate coloumns in model data
    kef_lat = 63.96
    kef_lon = -22.60
    polelat = theta_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = theta_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    kef_x, kef_y = rotate_pole(np.array([kef_lon]), np.array([kef_lat]), polelon, polelat)
    kef_x = kef_x + 360
    #print(kef_y, kef_x)
    sample_points = [('grid_longitude', kef_x), ('grid_latitude', kef_y)]
    theta_coloumns = trajectory.interpolate(theta_cube,sample_points, method =  'nearest')
    wind_coloumns = trajectory.interpolate(magwind_cube,sample_points, method =  'nearest')
    #%%
    model_h = theta_coloumns.coord('altitude').points
    #print(np.shape(model_h))
    #print(np.shape(theta_coloumns.data[0]))
    
    #%% stat calculations
   
    ## obs
    obs_mean_theta = np.mean(obs_theta)
    obs_std_theta = np.std(obs_theta)
    obs_mean_speed = np.mean(obs_speed)
    obs_std_speed = np.std(obs_speed)
    obs_gradient = (obs_theta[-1]-obs_theta[0])/(obs_h[-1]-obs_h[0])
    obs_grad_err_p = np.abs((obs_theta[-1]-obs_theta[1])/(obs_h[-1]-obs_h[1])  - obs_gradient)
    obs_grad_err_m = np.abs((obs_theta[-2]-obs_theta[0])/(obs_h[-2]-obs_h[0])  - obs_gradient)
    obs_grad_err = (obs_grad_err_p + obs_grad_err_m)/2
    #print('Obs')
    #print('number of obs points =', len(obs_h))
    print(obs_time, 'obs theta',obs_mean_theta,'+/-', obs_std_theta)
    #print(obs_mean_speed, obs_std_speed)
    print(obs_time,'obs gradient',obs_gradient, '+/-',obs_grad_err)
    ## UM
    um_sample_cond = (model_h < 2000)*(model_h > 200)
    U_sample_cond = (model_h < 1500)*(model_h > 500)
    um_h = model_h[um_sample_cond]
    um_theta = theta_coloumns[info[1]].data[um_sample_cond]
    
    um_speed = wind_coloumns[info[1]].data[U_sample_cond]
    
    um_mean_theta = np.mean(um_theta)
    um_std_theta = np.std(um_theta)
    um_mean_speed = np.mean(um_speed)
    um_std_speed = np.std(um_speed)
    um_gradient = (um_theta[-1]-um_theta[0])/(um_h[-1]-um_h[0])
    um_grad_err_p = np.abs((um_theta[-1]-um_theta[1])/(um_h[-1]-um_h[1])  - um_gradient)
    um_grad_err_m = np.abs((um_theta[-2]-um_theta[0])/(um_h[-2]-um_h[0])  -um_gradient)
    um_grad_err = (um_grad_err_p + um_grad_err_m)/2
    #print('UM')
    #print('number of obs points =', len(um_h))
    print(obs_time, 'um theta',um_mean_theta,'+/-', um_std_theta)
    #print(um_mean_speed, um_std_speed)
    print(obs_time, 'um gradient',um_gradient,'+/-',um_grad_err)
    #%% perform Froude number calcualtion
    h = 1000
    
    obs_N = calc_bouancy_freq(um_mean_theta, obs_gradient)
    um_N = calc_bouancy_freq(um_mean_theta, um_gradient)
    
    obsNe = calc_N_error(obs_N, obs_mean_theta, obs_std_theta, obs_gradient, obs_grad_err)
    umNe = calc_N_error(um_N, um_mean_theta, um_std_theta, um_gradient, um_grad_err)
    
    obs_Fr = obs_mean_speed/(obs_N*h)
    um_Fr = um_mean_speed/(um_N*h)
    
    obsFre = calc_Fr_error(obs_Fr, obs_N, obsNe, obs_mean_speed, obs_std_speed)
    umFre = calc_Fr_error(um_Fr, um_N, umNe, um_mean_speed, um_std_speed)
    
    print(obs_time, 'N freq', obs_N, um_N)
    print(obsNe, umNe)
    print(obs_time,'Fr numbers',obs_Fr, um_Fr)
    print(obsFre,umFre)
    
    #%%
    fig, ax = plt.subplots(figsize = (4,10))
    
    ax.scatter(obs_theta, obs_h)
    #ax.set_xlim(right = 302)
    ax.set_ylim(bottom = 0,top = 5000)
    ax.scatter(theta_coloumns.data[35, :20], model_h[:20])