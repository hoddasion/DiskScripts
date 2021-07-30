#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 17:08:43 2020

@author: kse18nru
"""

#%% module imports
import iris
from iris.analysis.cartography import rotate_pole
import data_management as dmna
import pandas as pd
import numpy as np
import persinterpolate as pert
import matplotlib.pyplot as plt
import xsection_workshop as workshop
#%%
resolutions = ['0p5km', '1p5km', '4p4km']
varname_u = 'm01s15i002'
varname_v = 'm01s15i003'
varname_wsp = 'windspeed'
varname_w = 'upward_air_velocity'
varname_theta = 'air_potential_temperature'
varname_q = 'specific_humidity'
leg = 8
plotleg = 8
lat_lon_ord = ('grid_latitude', 'grid_longitude')
lon_lat_ord = ('grid_longitude', 'grid_latitude')
leg_meta = ((1,26,1, lon_lat_ord, (0,-1)),(2,27,2,lon_lat_ord, (0,-1)),(3,28,3,lon_lat_ord, (0,-5)),(4,29,3,lon_lat_ord, (5,-1)),
            (5,30,3, lon_lat_ord, (0,-1)),(6,30,3, lon_lat_ord, (5,-1)),(8,31,8, lon_lat_ord, (0,-1)),(9,32,8, lon_lat_ord,(0,-1)))

flight = 301
config = 'RA1M'
exp = 'Control'
suite = 'u-bu807'
standard_plotting = False
if standard_plotting:
    for fileleg in leg_meta:
        for res in resolutions:
            #fileleg = leg_meta[0]
            #%% load and process core model diagnostics from vertical sliced files
            
            wsp_cube = dmna.load_model_wsp(res, varname_u, varname_v, 'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
            w_cube = dmna.load_model_w(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
            q_cube = dmna.load_model_specific_hum(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
            theta_cube = dmna.load_model_theta(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35, secondary_var = 'upward_air_velocity')
            #%% load and process orography
            
            orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
            orogmax = orogcube.collapsed(fileleg[3][1], iris.analysis.MAX)
            orogdata = orogmax.data
            orogcoor = orogmax.coord(fileleg[3][0]).points
            
            #%% load and process databse shortrun obs data
            # using pandas
            database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_sr_301.txt', delimiter = ' ')
            
            db_legno = database['legno']
            #legcond = np.array(db_legno == 3) + np.array(db_legno == 4) + np.array(db_legno == 5) + np.array(db_legno == 6)
            legcond = db_legno == fileleg[0]
            left = fileleg[4][0]; right = fileleg[4][1]
            db_wsp = np.array(database['wsp'][legcond])[left:right] 
            db_w = np.array(database['w'][legcond])[left:right]
            db_q = np.array(database['q_BUCK'][legcond])[left:right]
            db_theta = np.array(database['theta'][legcond])[left:right]
            
            db_lon = np.array(database['lon'][legcond])[left:right]
            db_lat = np.array(database['lat'][legcond])[left:right]
            db_alt = np.array(database['altgps'][legcond])[left:right]
            
            #%% rotate db coords
            compcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/RA1M_0p5km_umh1_flt{flight}.nc', varname_u)[0,0]
            polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
            polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
            rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
            rot_db_lon = rot_db_lon + 360
            
            #%% interpolate model data onto obs altitudes
            time_index = fileleg[1]
            if fileleg[3][0] == 'grid_longitude':
                matched_wsp, matched_coors = pert.interpolate_matched_series(wsp_cube, db_alt, rot_db_lon, time_index)
                simple_wsp, simple_coors_wsp = pert.interpolate_simple_series(wsp_cube, np.nanmean(db_alt), time_index)
                
                matched_w, matched_coors = pert.interpolate_matched_series(w_cube, db_alt, rot_db_lon, time_index)
                simple_w, simple_coors_w = pert.interpolate_simple_series(w_cube, np.nanmean(db_alt), time_index)
                
                matched_q, matched_coors = pert.interpolate_matched_series(q_cube, db_alt, rot_db_lon, time_index)
                simple_q, simple_coors_q = pert.interpolate_simple_series(q_cube, np.nanmean(db_alt), time_index)
                
                matched_theta, matched_coors = pert.interpolate_matched_series(theta_cube, db_alt, rot_db_lon, time_index)
                simple_theta, simple_coors_t = pert.interpolate_simple_series(theta_cube, np.nanmean(db_alt), time_index)
                db_coor = rot_db_lon
            elif fileleg[3][0] == 'grid_latitude':
                matched_wsp, matched_coors = pert.interpolate_matched_series(wsp_cube, db_alt, rot_db_lat, time_index, grid_axis = 'grid_latitude')
                simple_wsp, simple_coors_wsp = pert.interpolate_simple_series(wsp_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                
                matched_w, matched_coors = pert.interpolate_matched_series(w_cube, db_alt, rot_db_lat, time_index, grid_axis = 'grid_latitude')
                simple_w, simple_coors_w = pert.interpolate_simple_series(w_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                
                matched_q, matched_coors = pert.interpolate_matched_series(q_cube, db_alt, rot_db_lat, time_index, grid_axis = 'grid_latitude')
                simple_q, simple_coors_q = pert.interpolate_simple_series(q_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                
                matched_theta, matched_coors = pert.interpolate_matched_series(theta_cube, db_alt, rot_db_lat, time_index, grid_axis = 'grid_latitude')
                simple_theta, simple_coors_t = pert.interpolate_simple_series(theta_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                db_coor = rot_db_lat
            
            #%% dynamically adjust horizontal bounds
            xmin = np.nanmin(simple_coors_w)
            xmax = np.nanmax(simple_coors_w)
            
            #%% get scalar distance ticklist for x axis labelling
            tick_coors, tick_dists = workshop.dist_tick_generator(w_cube, fileleg[3][0], ticknumber = 6, decprecision = 2)
            #%% format date time for file name and figure title
            time = w_cube.coord('time')
            dates = time.units.num2date(time.points)
            print(dates[fileleg[1]])
            #%% plot some stuff
            fig, ((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2, figsize = (20,20))
            
            
            ## wsp
            ax0.scatter(db_wsp, matched_wsp, label = 'UM wsp')
            
            ax0.set_ylabel(f'UM, {w_cube.units}')
            ax0.set_ylim(bottom = 0, top = 16)
            #plt.title(f'Atmospheric state variables, leg {fileleg[0]} xsection, {dates[fileleg[1]]}, {(np.nanmean(db_alt)//1)}m altitude')
            ax0.set_xlim(left = 0, right = 16)
            ax0.set_xlabel(f'Obs, {w_cube.units}')
            ax0.set_title('windspeed')
            ax0.grid()
            ## w
            ax1.scatter(db_w, matched_w)
            ax1.set_ylim(bottom = -2, top = 2)
            ax1.set_xlim(left = -2, right = 2)
            ax1.set_ylabel(f'UM, {w_cube.units}')
            ax1.set_xlabel(f'Obs, {w_cube.units}')
            ax1.set_title('vertical windspeed')
            ax1.grid()
            ## theta
            ax2.scatter(db_theta, matched_theta)
            ax2.set_ylim(bottom = 270, top = 290)
            ax2.set_xlim(left = 270, right = 290)
            ax2.set_ylabel(f'UM, {theta_cube.units}')
            ax2.set_xlabel(f'Obs, {theta_cube.units}')
            ax2.set_title('potential temperature')
            ax2.grid()
            ## specfific humidity
            ax3.scatter(db_q, matched_q)
            ax3.set_ylim(bottom = 0, top = 3)
            ax3.set_xlim(left = 0, right = 3)
            ax3.set_ylabel(r'UM, $g kg^{-1}$')
            ax3.set_xlabel(r'Obs, $g kg^{-1}$')
            ax3.set_title('sepcific humidity')
            ax3.grid()
            
            
            fig.suptitle(f'Atmospheric state variables, leg {fileleg[0]} xsection, {dates[fileleg[1]]}, {(np.nanmean(db_alt)//1)}m mean altitude')
            fig.tight_layout()
            
            save_fig = True
            if save_fig:
                timepoint = dates[fileleg[1]].strftime('H%HM%M')
                fig.savefig(f'../../Figures/PNG/301/{suite}/Scatter/{exp}/atmostate/{res}/{config}_{res}_atmostate_scatter_{exp}_{timepoint}_leg{fileleg[0]}.png')
                
boundary_layer = True
firstleg = 1; secondleg = 7
if boundary_layer:
    fig, ((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2, figsize = (20,20))
    for res in resolutions:
        if res == '0p5km' or '4p4km':
            fileleg = leg_meta[firstleg]
            store_leg = fileleg
            #%% load and process core model diagnostics from vertical sliced files
            
            wsp_cube = dmna.load_model_wsp(res, varname_u, varname_v, 'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
            w_cube = dmna.load_model_w(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
            q_cube = dmna.load_model_specific_hum(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
            theta_cube = dmna.load_model_theta(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35, secondary_var = 'upward_air_velocity')
            #%% load and process orography
            
            orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
            orogmax = orogcube.collapsed(fileleg[3][1], iris.analysis.MAX)
            orogdata = orogmax.data
            orogcoor = orogmax.coord(fileleg[3][0]).points
            
            #%% load and process databse shortrun obs data
            # using pandas
            database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_sr_301.txt', delimiter = ' ')
            
            db_legno = database['legno']
            #legcond = np.array(db_legno == 3) + np.array(db_legno == 4) + np.array(db_legno == 5) + np.array(db_legno == 6)
            legcond = db_legno == fileleg[0]
            left = fileleg[4][0]; right = fileleg[4][1]
            db_wsp = np.array(database['wsp'][legcond])[left:right] 
            db_w = np.array(database['w'][legcond])[left:right]
            db_q = np.array(database['q_BUCK'][legcond])[left:right]
            db_theta = np.array(database['theta'][legcond])[left:right]
            
            db_lon = np.array(database['lon'][legcond])[left:right]
            db_lat = np.array(database['lat'][legcond])[left:right]
            db_alt = np.array(database['altgps'][legcond])[left:right]
            
            #%% rotate db coords
            compcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/RA1M_0p5km_umh1_flt{flight}.nc', varname_u)[0,0]
            polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
            polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
            rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
            rot_db_lon = rot_db_lon + 360
            
            #%% interpolate model data onto obs altitudes
            time_index = fileleg[1]
            if fileleg[3][0] == 'grid_longitude':
                matched_wsp, matched_coors = pert.interpolate_matched_series(wsp_cube, db_alt, rot_db_lon, time_index)
                simple_wsp, simple_coors_wsp = pert.interpolate_simple_series(wsp_cube, np.nanmean(db_alt), time_index)
                
                matched_w, matched_coors = pert.interpolate_matched_series(w_cube, db_alt, rot_db_lon, time_index)
                simple_w, simple_coors_w = pert.interpolate_simple_series(w_cube, np.nanmean(db_alt), time_index)
                
                matched_q, matched_coors = pert.interpolate_matched_series(q_cube, db_alt, rot_db_lon, time_index)
                simple_q, simple_coors_q = pert.interpolate_simple_series(q_cube, np.nanmean(db_alt), time_index)
                
                matched_theta, matched_coors = pert.interpolate_matched_series(theta_cube, db_alt, rot_db_lon, time_index)
                simple_theta, simple_coors_t = pert.interpolate_simple_series(theta_cube, np.nanmean(db_alt), time_index)
                db_coor = rot_db_lon
            elif fileleg[3][0] == 'grid_latitude':
                matched_wsp, matched_coors = pert.interpolate_matched_series(wsp_cube, db_alt, rot_db_lat, time_index, grid_axis = 'grid_latitude')
                simple_wsp, simple_coors_wsp = pert.interpolate_simple_series(wsp_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                
                matched_w, matched_coors = pert.interpolate_matched_series(w_cube, db_alt, rot_db_lat, time_index, grid_axis = 'grid_latitude')
                simple_w, simple_coors_w = pert.interpolate_simple_series(w_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                
                matched_q, matched_coors = pert.interpolate_matched_series(q_cube, db_alt, rot_db_lat, time_index, grid_axis = 'grid_latitude')
                simple_q, simple_coors_q = pert.interpolate_simple_series(q_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                
                matched_theta, matched_coors = pert.interpolate_matched_series(theta_cube, db_alt, rot_db_lat, time_index, grid_axis = 'grid_latitude')
                simple_theta, simple_coors_t = pert.interpolate_simple_series(theta_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                db_coor = rot_db_lat
            
            #%% dynamically adjust horizontal bounds
            xmin = np.nanmin(simple_coors_w)
            xmax = np.nanmax(simple_coors_w)
            mean_alt1 = np.nanmean(db_alt)//1
            #%% get scalar distance ticklist for x axis labelling
            tick_coors, tick_dists = workshop.dist_tick_generator(w_cube, fileleg[3][0], ticknumber = 6, decprecision = 2)
            #%% format date time for file name and figure title
            time = w_cube.coord('time')
            dates = time.units.num2date(time.points)
            print(dates[fileleg[1]])
            #%% plot some stuff
            round1_keys = {'marker' : 'x', 's' : 300}
            #fig, ((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2, figsize = (20,20))
            
            
            ## wsp
            ax0.scatter(db_wsp, matched_wsp, label = f'{res} leg {fileleg[0]}', **round1_keys)
            
            ax0.set_ylabel(f'UM, {w_cube.units}')
            ax0.set_ylim(bottom = 0, top = 16)
            #plt.title(f'Atmospheric state variables, leg {fileleg[0]} xsection, {dates[fileleg[1]]}, {(np.nanmean(db_alt)//1)}m altitude')
            ax0.set_xlim(left = 0, right = 16)
            ax0.set_xlabel(f'Obs, {w_cube.units}')
            ax0.set_title('windspeed')
            ax0.grid()
            ## w
            ax1.scatter(db_w, matched_w, label = f'{res} leg {fileleg[0]}', **round1_keys)
            ax1.set_ylim(bottom = -2, top = 2)
            ax1.set_xlim(left = -2, right = 2)
            ax1.set_ylabel(f'UM, {w_cube.units}')
            ax1.set_xlabel(f'Obs, {w_cube.units}')
            ax1.set_title('vertical windspeed')
            ax1.grid()
            ## theta
            ax2.scatter(db_theta, matched_theta, label = f'{res} leg {fileleg[0]}', **round1_keys)
            ax2.set_ylim(bottom = 270, top = 290)
            ax2.set_xlim(left = 270, right = 290)
            ax2.set_ylabel(f'UM, {theta_cube.units}')
            ax2.set_xlabel(f'Obs, {theta_cube.units}')
            ax2.set_title('potential temperature')
            ax2.grid()
            ## specfific humidity
            ax3.scatter(db_q, matched_q, label = f'{res} leg {fileleg[0]}', **round1_keys)
            ax3.set_ylim(bottom = 0, top = 3)
            ax3.set_xlim(left = 0, right = 3)
            ax3.set_ylabel(r'UM, $g kg^{-1}$')
            ax3.set_xlabel(r'Obs, $g kg^{-1}$')
            ax3.set_title('sepcific humidity')
            ax3.grid()
            
            
            
            ### now do it all again
            fileleg = leg_meta[secondleg]
            #%% load and process core model diagnostics from vertical sliced files
            
            wsp_cube = dmna.load_model_wsp(res, varname_u, varname_v, 'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
            w_cube = dmna.load_model_w(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
            q_cube = dmna.load_model_specific_hum(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
            theta_cube = dmna.load_model_theta(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35, secondary_var = 'upward_air_velocity')
            #%% load and process orography
            
            orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
            orogmax = orogcube.collapsed(fileleg[3][1], iris.analysis.MAX)
            orogdata = orogmax.data
            orogcoor = orogmax.coord(fileleg[3][0]).points
            
            #%% load and process databse shortrun obs data
            # using pandas
            database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_sr_301.txt', delimiter = ' ')
            
            db_legno = database['legno']
            #legcond = np.array(db_legno == 3) + np.array(db_legno == 4) + np.array(db_legno == 5) + np.array(db_legno == 6)
            legcond = db_legno == fileleg[0]
            left = fileleg[4][0]; right = fileleg[4][1]
            db_wsp = np.array(database['wsp'][legcond])[left:right] 
            db_w = np.array(database['w'][legcond])[left:right]
            db_q = np.array(database['q_BUCK'][legcond])[left:right]
            db_theta = np.array(database['theta'][legcond])[left:right]
            
            db_lon = np.array(database['lon'][legcond])[left:right]
            db_lat = np.array(database['lat'][legcond])[left:right]
            db_alt = np.array(database['altgps'][legcond])[left:right]
            
            #%% rotate db coords
            compcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/RA1M_0p5km_umh1_flt{flight}.nc', varname_u)[0,0]
            polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
            polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
            rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
            rot_db_lon = rot_db_lon + 360
            
            #%% interpolate model data onto obs altitudes
            time_index = fileleg[1]
            if fileleg[3][0] == 'grid_longitude':
                matched_wsp, matched_coors = pert.interpolate_matched_series(wsp_cube, db_alt, rot_db_lon, time_index)
                simple_wsp, simple_coors_wsp = pert.interpolate_simple_series(wsp_cube, np.nanmean(db_alt), time_index)
                
                matched_w, matched_coors = pert.interpolate_matched_series(w_cube, db_alt, rot_db_lon, time_index)
                simple_w, simple_coors_w = pert.interpolate_simple_series(w_cube, np.nanmean(db_alt), time_index)
                
                matched_q, matched_coors = pert.interpolate_matched_series(q_cube, db_alt, rot_db_lon, time_index)
                simple_q, simple_coors_q = pert.interpolate_simple_series(q_cube, np.nanmean(db_alt), time_index)
                
                matched_theta, matched_coors = pert.interpolate_matched_series(theta_cube, db_alt, rot_db_lon, time_index)
                simple_theta, simple_coors_t = pert.interpolate_simple_series(theta_cube, np.nanmean(db_alt), time_index)
                db_coor = rot_db_lon
            elif fileleg[3][0] == 'grid_latitude':
                matched_wsp, matched_coors = pert.interpolate_matched_series(wsp_cube, db_alt, rot_db_lat, time_index, grid_axis = 'grid_latitude')
                simple_wsp, simple_coors_wsp = pert.interpolate_simple_series(wsp_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                
                matched_w, matched_coors = pert.interpolate_matched_series(w_cube, db_alt, rot_db_lat, time_index, grid_axis = 'grid_latitude')
                simple_w, simple_coors_w = pert.interpolate_simple_series(w_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                
                matched_q, matched_coors = pert.interpolate_matched_series(q_cube, db_alt, rot_db_lat, time_index, grid_axis = 'grid_latitude')
                simple_q, simple_coors_q = pert.interpolate_simple_series(q_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                
                matched_theta, matched_coors = pert.interpolate_matched_series(theta_cube, db_alt, rot_db_lat, time_index, grid_axis = 'grid_latitude')
                simple_theta, simple_coors_t = pert.interpolate_simple_series(theta_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                db_coor = rot_db_lat
            
            #%% dynamically adjust horizontal bounds
            xmin = np.nanmin(simple_coors_w)
            xmax = np.nanmax(simple_coors_w)
            mean_alt2 = np.nanmean(db_alt)//1
            #%% get scalar distance ticklist for x axis labelling
            tick_coors, tick_dists = workshop.dist_tick_generator(w_cube, fileleg[3][0], ticknumber = 6, decprecision = 2)
            #%% format date time for file name and figure title
            time = w_cube.coord('time')
            dates = time.units.num2date(time.points)
            print(dates[fileleg[1]])
            
            #%%
            round2_keys = {'marker' : '+', 's' : 300}
            ax0.scatter(db_wsp, matched_wsp, label = f'{res} leg {fileleg[0]}', **round2_keys)
            ax1.scatter(db_w, matched_w, label = f'{res} leg {fileleg[0]}', **round2_keys)
            ax2.scatter(db_theta, matched_theta, label = f'{res} leg {fileleg[0]}', **round2_keys)
            ax3.scatter(db_q, matched_q, label = f'{res} leg {fileleg[0]}', **round2_keys)
            
            ax0.legend()
            ax1.legend()
            ax2.legend()
            ax3.legend()
            
            fig.suptitle(f'Atmospheric state variables, legs {store_leg[0]} & {fileleg[0]} xsection, {dates[store_leg[0]]} & {dates[fileleg[1]]}, {mean_alt1} & {mean_alt2}m mean altitude')
            fig.tight_layout() 
            
    plt.show()
    save_fig = True
    if save_fig:
        timepoint = dates[fileleg[1]].strftime('H%HM%M')
        fig.savefig(f'../../Figures/PNG/301/{suite}/Scatter/{exp}/atmostate/{config}_0p5km_4p4km_atmostate_scatter_{exp}_leg{store_leg[0]}n{fileleg[0]}.png')
            