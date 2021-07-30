# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 16:46:37 2020

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
from model_foundry import toolkit
import datetime
import xsection_workshop as workshop

#%% execute function
legs = [9]
resolutions = ['0p5km', '1p5km', '4p4km']
vmin = -600; vmax = - vmin
varname = 'upward_heat_flux_in_air'
#%%
Test = False
if Test:
    D = '0p5km'; leg = 1
    Case = {'res' : D, 'flight' : 301, 'varname' : varname, 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
    res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
    fileleg = Case['fileleg']
    
    maincube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', varname)[:,:30]
    thetacube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
    orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
    print(maincube, '\n', thetacube, '\n', orogcube)
    print(maincube.data)
    #%%
    fig = workshop.plot_vertical_xsection(maincube, thetacube, orogcube,index = 6, Case = Case, cmap = 'seismic', vmin = vmin,vmax = vmax, levels = 20)
    plt.show()
#%%
run_loop = True
if run_loop:
    for D in resolutions:
        for leg in legs:
            try:
                #%%
            
                Case = {'res' : D, 'flight' : 301, 'varname' : varname, 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
                res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
                fileleg = Case['fileleg']
                print(Case.values())
                #%%
                
                maincube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', varname)[:,:23]
                thetacube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
                orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
                print(maincube, '\n', thetacube, '\n', orogcube)
                print(maincube.units)
            except:
                print('File error')
                continue
            path = f'../../Figures/PNG/301/u-bu807/Vertical/{varname}/Alt2000/leg{fileleg}/{res}'
            for index in range(47):
                
        
                
                fig = workshop.plot_vertical_xsection(maincube, thetacube, orogcube,index, Case, axes = ('grid_longitude', 'grid_latitude'), cmap = 'seismic', vmin = vmin,vmax = vmax, levels = 20,alt_limit = 2000)
                workshop.save_vert_xsection(fig, maincube,index, path, Case, fileleg = leg, xaxis = 'lat')
                #plt.close()
                plt.show()
                
obs_overlay = False
if obs_overlay:
    #%%
    D = '0p5km'; leg = 1
    Case = {'res' : D, 'flight' : 301, 'varname' : varname, 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
    res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
    fileleg = Case['fileleg']
    
    
    maincube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc',varname)[:,:30]
    
    thetacube = maincube#iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
    orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
    print(maincube, '\n', orogcube)
    print(maincube.long_name)
    print('maximum',np.nanmax(maincube.data))
    compcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/RA1M_0p5km_umi1_flt{flight}.nc', varname)
    polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    print(polelat)
    #%%
    ## load quality controlled and interval-meaned observational data from Database
    # using pandas
    database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_60s_301.txt', delimiter = ' ')
    
    db_legno = database['legno']
    #legcond = np.array(db_legno == 3) + np.array(db_legno == 4) + np.array(db_legno == 5) + np.array(db_legno == 6)
    legcond = db_legno == leg
    db_latent = np.array(database['lh'][legcond])
    db_sens = np.array(database['sh'][legcond])
    db_turb = np.array(database['turbh'][legcond])
    db_flux = db_latent + db_sens + db_turb
   
    db_lon = np.array(database['lon'][legcond])
    db_lat = np.array(database['lat'][legcond])
    db_alt = np.array(database['altgps'][legcond])
    
    db_mtime = np.array(database['meantime'][legcond])
    
    #%% rotate db coords

    rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
    rot_db_lon = rot_db_lon + 360
    
    
    #%%
    t_indices = [26]#[28,29,30]
    cmap = 'seismic'
    for index in t_indices:
        fig = workshop.plot_vertical_xsection(maincube, thetacube, orogcube,
                                              index = index, Case = Case, cmap = cmap, vmin = vmin,vmax = vmax, levels = 20, 
                                              theta_contour = False,alt_limit = 2000)
        plt.scatter(rot_db_lon, db_alt,c=db_sens, s = 300, cmap = cmap, 
                    vmin = vmin,vmax = vmax, edgecolor = 'k')
        path = f'../../Figures/PNG/301/u-bu807/Vertical/{varname}/Alt2000/Obs/{res}'
        #orkshop.save_vert_xsection(fig, maincube,index, path, Case, fileleg = leg, obs = True, obs_detail = 'l1')
        #plt.close()
        plt.show()