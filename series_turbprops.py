#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 15:03:19 2020

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
import sys
#%%
plt.rcParams['font.size'] = 26
resolutions = ['0p5km', '1p5km', '4p4km']


varname_tke = 'm01s03i473'
varname_wsn = 'atmosphere_downward_northward_stress'
varname_wse = 'atmosphere_downward_eastward_stress'
varname_uhf = 'upward_heat_flux_in_air'
varname_umf = 'upward_water_vapor_flux_in_air' # needed for latent heat
varname_p   = 'air_pressure' # needed for latent heat
varname_q   = 'specific_humidity' # needed for latent heat
varname_theta = 'air_potential_temperature' # needed for latent heat

lat_lon_ord = ('grid_latitude', 'grid_longitude')
lon_lat_ord = ('grid_longitude', 'grid_latitude')
leg_meta_obs = ((1,26,1, lon_lat_ord),(2,27,2,lon_lat_ord),(3,28,3,lon_lat_ord),(4,29,3,lon_lat_ord),
            (5,30,3, lon_lat_ord),(6,30,3, lon_lat_ord),(8,31,8, lat_lon_ord),(9,32,8, lat_lon_ord))
leg_meta_time = ((1,26,1, lon_lat_ord),(2,27,2,lon_lat_ord),(3,28,3,lon_lat_ord),(8,31,8, lat_lon_ord))
flight = 301
config = 'RA1M'
exp = 'Control'
suite = 'u-bu807'

for fileleg in leg_meta_obs:
    for res in resolutions:
        
        #%%
        #fileleg = (1,26,1, lon_lat_ord)
        #res = '0p5km'
        plt.close()
        wst_cube = dmna.load_model_windstress(res, varname_wse, varname_wsn, 'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
        p_cube = dmna.load_model_pressure(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
        q_cube = dmna.load_model_specific_hum(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
        theta_cube = dmna.load_model_theta(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35, secondary_var = 'upward_air_velocity')
        vapor_cube = dmna.load_model_vapor_flux(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
        sh_cube = dmna.load_model_heat_flux(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
        lh_cube = dmna.load_model_heat_flux(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
        #tke_cube = dmna.load_model_tke(res,'Vertical_slices', fileleg[2],301, 'Control')
        #print(vapor_cube)
        #print(tke_cube)
        
        #%% get coord system and surface altitude out for later
        cube_coor = sh_cube.coord(fileleg[3][0]).points
        surf_alt = sh_cube.coord('surface_altitude').points
        #print(surf_alt)
        #%% extract coordinates
        
        model_coors = p_cube.coord(fileleg[3][0])
        #%% load and process orography
        
        orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
        orogmax = orogcube.collapsed(fileleg[3][1], iris.analysis.MAX)
        orogdata = orogmax.data
        orogcoor = orogmax.coord(fileleg[3][0]).points
        
        #%% load and process databse shortrun obs data
        # using pandas
        database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_sr_301.txt', delimiter = ' ')
        leg = fileleg[0]
        db_legno = database['legno']
        print(database)
        legcond = db_legno == leg
        db_wst = np.array(database['windstress'][legcond]) 
        db_tke = np.array(database['tke'][legcond])
        db_sh = np.array(database['sh'][legcond])
        db_lh = np.array(database['lh'][legcond])
        
        db_lon = np.array(database['lon'][legcond])
        db_lat = np.array(database['lat'][legcond])
        db_alt = np.array(database['altgps'][legcond])
        
        #%% rotate db coords
        compcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/RA1M_0p5km_umh1_flt{flight}.nc', 'upward_air_velocity')[0,0]
        polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
        polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
        rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
        rot_db_lon = rot_db_lon + 360

        #%% interpolate model data onto obs altitudes
        time_index = fileleg[1]
        if fileleg[3][0] == 'grid_longitude':
            matched_wst, matched_coors = pert.interpolate_matched_series(wst_cube, db_alt, rot_db_lon, time_index)
            simple_wst, simple_coors_wst = pert.interpolate_simple_series(wst_cube, np.nanmean(db_alt), time_index)
            
            matched_sh, matched_coors = pert.interpolate_matched_series(sh_cube, db_alt, rot_db_lon, time_index)
            simple_sh, simple_coors_sh = pert.interpolate_simple_series(sh_cube, np.nanmean(db_alt), time_index)
            
            matched_lh, matched_coors = pert.interpolate_matched_series(lh_cube, db_alt, rot_db_lon, time_index)
            simple_lh, simple_coors_lh = pert.interpolate_simple_series(lh_cube, np.nanmean(db_alt), time_index)
            
            
            db_coor = rot_db_lon
        elif fileleg[3][0] == 'grid_latitude':
            matched_wst, matched_coors = pert.interpolate_matched_series(wst_cube, db_alt, rot_db_lat, time_index, grid_axis = 'grid_latitude')
            simple_wst, simple_coors_wst = pert.interpolate_simple_series(wst_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
            
            matched_sh, matched_coors = pert.interpolate_matched_series(sh_cube, db_alt, rot_db_lat, time_index, grid_axis = 'grid_latitude')
            simple_sh, simple_coors_sh = pert.interpolate_simple_series(sh_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
            
            matched_lh, matched_coors = pert.interpolate_matched_series(lh_cube, db_alt, rot_db_lat, time_index, grid_axis = 'grid_latitude')
            simple_lh, simple_coors_lh = pert.interpolate_simple_series(lh_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
            
            
            db_coor = rot_db_lat
            
        #%% create uniform altitude array for indication on plot
        uni_alt = np.ones(len(cube_coor))*np.nanmean(db_alt)
        #%% dynamically adjust horizontal bounds
        xmin = np.nanmin(simple_coors_sh)
        xmax = np.nanmax(simple_coors_sh)
        
        #%% get scalar distance ticklist for x axis labelling
        tick_coors, tick_dists = workshop.dist_tick_generator(sh_cube, fileleg[3][0], ticknumber = 6, decprecision = 2)
        #%% format date time for file name and figure title
        time = sh_cube.coord('time')
        dates = time.units.num2date(time.points)
        print(dates[fileleg[1]])
        #%% plot some stuff
        fig, (ax0,ax1,ax2,ax4) = plt.subplots(4,1, figsize = (20,17))
        
        ## orographic panel (bottom, ax4)
        ax4.plot(orogcoor,orogdata, label = 'Model background orography', color = 'k')
        ax4.plot(cube_coor, surf_alt, linestyle = '--', label = 'Model orography directly below', color = 'k')
        ax4.plot(cube_coor, uni_alt, color = 'g', label = 'Interpolation surface')
        ax4.plot(db_coor, db_alt, color = 'g', linestyle = '--', label = 'Obs gps')
        ax4.set_ylabel('m')
        ax4.set_ylim(bottom = 0, top = 1700)
        ax4.legend()
        ax4.set_xlim(left = xmin, right = xmax)
        ax4.set_xticks(tick_coors)
        ax4.set_xticklabels(tick_dists)
        ax4.set_xlabel('Scalar surface distance [km]')
        ax4.grid()
        ## wsp
        ax0.plot(simple_coors_wst, simple_wst, label = 'UM windstress magn.')
        ax0.plot(db_coor,db_wst, label = 'obs')
        ax0.set_ylabel(r'$Nm^{-1}$')
        ax0.set_ylim(bottom = 0, top = 4)
        ax0.set_xlim(left = xmin, right = xmax)
        ax0.set_xticks(tick_coors)
        ax0.set_xticklabels([])
        ax0.set_title(f'Atmospheric state variables, leg {fileleg[0]} xsection, {dates[fileleg[1]]}, {(np.nanmean(db_alt)//1)}m altitude')
        ax0.legend()
        ax0.grid()
        ## w
        ax1.plot(simple_coors_sh, simple_sh, label = 'UM upward SH')
        ax1.plot(db_coor,db_sh, label = 'obs')
        ax1.set_ylabel(sh_cube.units)
        ax1.set_ylim(bottom = -200, top = 200)
        ax1.set_xlim(left = xmin, right = xmax)
        ax1.set_xticks(tick_coors)
        ax1.set_xticklabels([])
        ax1.legend()
        ax1.grid()
        ## theta
        ax2.plot(simple_coors_lh, simple_lh, label = 'UM upward LH')
        ax2.plot(db_coor, db_lh, label = 'obs')
        ax2.set_ylabel(sh_cube.units)
        ax2.set_ylim(bottom = -200, top = 200)
        ax2.set_xlim(left = xmin, right = xmax)
        ax2.set_xticks(tick_coors)
        ax2.set_xticklabels([])
        ax2.legend()
        ax2.grid()
        
        
        
        fig.tight_layout()
        
        save_fig = True
        if save_fig:
            timepoint = dates[fileleg[1]].strftime('H%HM%M')
            fig.savefig(f'../../Figures/PNG/301/{suite}/Series/{exp}/turbprops/{res}/{config}_{res}_turbprops_series_{exp}_{timepoint}_leg{fileleg[0]}.png')
        plt.show()
        #sys.exit()