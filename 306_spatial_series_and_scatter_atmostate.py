# -*- coding: utf-8 -*-
"""
Created on Tue May 18 22:44:56 2021

@author: kse18nru
"""

#%% module imports
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

resolutions = ['0p5km','1p5km','4p4km']
suite = 'u-cc134'
#expt = 'Control'
config = 'RA1M'
flight = 306
spatial_60s_atmostate = False
if spatial_60s_atmostate:
    # work for 60s averaged database obs
    print(0)
    for res in resolutions:
        print(f'\n\nRes {res}')
        #%% load model data
        ### 3D whole cubes
        alt_idx = 35
        u_3dcube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{config}_{res}_um_x_wind_24hrs_ph_306.nc', 'x_wind')[8:-4,:alt_idx]
        v_3dcube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{config}_{res}_um_y_wind_24hrs_ph_306.nc', 'y_wind')[8:-4,:alt_idx]
        theta_3dcube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{config}_{res}_um_air_potential_temperature_24hrs_pi_306.nc', 'air_potential_temperature')[8:-4,:alt_idx]
        w_3dcube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{config}_{res}_um_upward_air_velocity_24hrs_ph_306.nc', 'upward_air_velocity')[8:-4,:alt_idx]
        q_3dcube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{config}_{res}_um_specific_humidity_24hrs_pi_306.nc', 'specific_humidity')[8:-4,:alt_idx]
        
        #%% process u and v into wsp
        wsp_3dcube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{config}_{res}_um_x_wind_24hrs_ph_306.nc', 'x_wind')[8:-4,:alt_idx]
        wsp_3dcube.data = (u_3dcube.data**2 + v_3dcube.data[:,:,:-1]**2)**0.5
        
        #%%
        #%% load and process database shortrun obs data
        # using pandas
        database = pd.read_csv(f'D:/Project/Obvs_Data/Databases/IGP_flights_database_obs_60s_{flight}.txt', delimiter = ' ')
        
        db_legno = database['legno']
        slice_nos = [1,4,7,13,15]
        #legcond = np.array(db_legno == 3) + np.array(db_legno == 4) + np.array(db_legno == 5) + np.array(db_legno == 6)
        for j in range(len(slice_nos)):
            print(f'Slice {j}')
            ## go through the primary x-sections first, and we will cycle through various obs legs later (those that are stacked ontop of each other)
            ## load sliced model data
            ### sliced cubes saved previously
            cube_leg = slice_nos[j]
            alt_idx = 35
            u_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/x_wind/{config}_vertslice_{res}_x_wind_flt306_leg{cube_leg}.nc', 'x_wind')[8:-4,:alt_idx]
            v_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/y_wind/{config}_vertslice_{res}_y_wind_flt306_leg{cube_leg}.nc', 'y_wind')[8:-4,:alt_idx]
            theta_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/air_potential_temperature/{config}_vertslice_{res}_air_potential_temperature_flt306_leg{cube_leg}.nc', 'air_potential_temperature')[8:-4,:alt_idx]
            q_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/specific_humidity/{config}_vertslice_{res}_specific_humidity_flt306_leg{cube_leg}.nc', 'specific_humidity')[8:-4,:alt_idx]
            w_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/upward_air_velocity/{config}_vertslice_{res}_upward_air_velocity_flt306_leg{cube_leg}.nc','upward_air_velocity')[8:-4,:alt_idx]
            
            ## load u cube as wsp cube and replace data
            wsp_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/x_wind/{config}_vertslice_{res}_x_wind_flt306_leg{cube_leg}.nc', 'x_wind')[8:-4,:alt_idx]
            x_data = u_cube.data; y_data = v_cube.data
            wsp_data = (x_data**2 + y_data**2)**0.5
            wsp_cube.data = wsp_data
            #%%
            dblegs = [[1,2],[4,5,6],[7,8,9],[13],[15]]
            for leg in dblegs[j]:
                print(f'Leg {leg}')
                legcond = db_legno == leg
                left = 0; right = -1
                db_wsp = np.array(database['wsp'][legcond])[left:right] 
                db_w = np.array(database['w'][legcond])[left:right]
                db_q = np.array(database['q_BUCK'][legcond])[left:right]
                db_theta = np.array(database['theta'][legcond])[left:right]
                
                db_lon = np.array(database['lon'][legcond])[left:right]
                db_lat = np.array(database['lat'][legcond])[left:right]
                db_alt = np.array(database['altgps'][legcond])[left:right]
                
                #%% rotate db coords
            
                polelat = u_3dcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
                polelon = u_3dcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
                rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
                rot_db_lon = rot_db_lon + 360
                
                #%% format date time and comute mean leg time
                time = w_cube.coord('time')
                dates = time.units.num2date(time.points)
                db_time = np.array(database['meantime'][legcond])
                mid_point = np.nanmean(db_time)
                model_time_ssm = (time.points-time.points[0]+4.5)*60*60 # get seconds since midnight from model output timesteps
                model_time_hsm = (time.points-time.points[0]+4.5) # hours since midnight
                print(model_time_hsm)
                #%%
                for tidx in range(len(model_time_hsm)):
                    # go through each time step and check which is valid in next line
                    diff = mid_point/(60*60) - model_time_hsm[tidx]
                    if (diff < 0.25) and (diff >= -0.25):
                        # if meantime falls within 15min of model output
                        print('Time: obs', mid_point/(60*60), ', model', model_time_hsm[tidx]);
                        ## perform interpolations for current timestep
                        try:
                            time_index = tidx
                            matched_wsp, matched_coors = pert.interpolate_matched_series(wsp_3dcube, db_alt, rot_db_lon, time_index, obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                            #simple_wsp, simple_coors_wsp = pert.interpolate_simple_series(wsp_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                            
                            matched_w, matched_coors = pert.interpolate_matched_series(w_3dcube, db_alt, rot_db_lon, time_index, obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                            #simple_w, simple_coors_w = pert.interpolate_simple_series(w_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                            
                            matched_q, matched_coors = pert.interpolate_matched_series(q_3dcube, db_alt, rot_db_lon, time_index,obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                            #simple_q, simple_coors_q = pert.interpolate_simple_series(q_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                            
                            matched_theta, matched_coors = pert.interpolate_matched_series(theta_3dcube, db_alt, rot_db_lon, time_index, obs_lats = rot_db_lat,grid_axis = 'grid_latitude', mode = '3d')
                            #simple_theta, simple_coors_t = pert.interpolate_simple_series(theta_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                            db_coor = rot_db_lat
                            ## perform interpolations for prior timestep
                            time_index = tidx - 2
                            matched_wsp_m1h, matched_coors = pert.interpolate_matched_series(wsp_3dcube, db_alt, rot_db_lon, time_index, obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                            #simple_wsp_m1h, simple_coors_wsp = pert.interpolate_simple_series(wsp_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                            
                            matched_w_m1h, matched_coors = pert.interpolate_matched_series(w_3dcube, db_alt, rot_db_lon, time_index, obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                            #simple_w_m1h, simple_coors_w = pert.interpolate_simple_series(w_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                            
                            matched_q_m1h, matched_coors = pert.interpolate_matched_series(q_3dcube, db_alt, rot_db_lon, time_index,obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                            #simple_q_m1h, simple_coors_q = pert.interpolate_simple_series(q_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                            
                            matched_theta_m1h, matched_coors = pert.interpolate_matched_series(theta_3dcube, db_alt, rot_db_lon, time_index, obs_lats = rot_db_lat,grid_axis = 'grid_latitude', mode = '3d')
                            #simple_theta_m1h, simple_coors_t = pert.interpolate_simple_series(theta_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                            ## perform interpolations for prior timestep
                            time_index = tidx + 2
                            matched_wsp_p1h, matched_coors = pert.interpolate_matched_series(wsp_3dcube, db_alt, rot_db_lon, time_index, obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                            #simple_wsp_p1h, simple_coors_wsp = pert.interpolate_simple_series(wsp_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                            
                            matched_w_p1h, matched_coors = pert.interpolate_matched_series(w_3dcube, db_alt, rot_db_lon, time_index, obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                            #simple_w_p1h, simple_coors_w = pert.interpolate_simple_series(w_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                            
                            matched_q_p1h, matched_coors = pert.interpolate_matched_series(q_3dcube, db_alt, rot_db_lon, time_index,obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                            #simple_q_p1h, simple_coors_q = pert.interpolate_simple_series(q_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                            
                            matched_theta_p1h, matched_coors = pert.interpolate_matched_series(theta_3dcube, db_alt, rot_db_lon, time_index, obs_lats = rot_db_lat,grid_axis = 'grid_latitude', mode = '3d')
                            #simple_theta_p1h, simple_coors_t = pert.interpolate_simple_series(theta_cube, np.nanmean(db_alt), time_index, grid_axis = 'grid_latitude')
                            
                            #%% create uniform altitude array for indication on plot
                            #uni_alt = np.ones(len(cube_coor))*np.nanmean(db_alt)
                            #%% dynamically adjust horizontal bounds
                            xmin = np.nanmin(matched_coors)-0.1
                            xmax = np.nanmax(matched_coors)+0.1
                            
                            #%% get scalar distance ticklist for x axis labelling
                            tick_coors, tick_dists = workshop.obs_dist_tick_gen(rot_db_lon,rot_db_lat, 'grid_latitude', ticknumber = 6, decprecision = 2)
                            
                            #%% plot some stuff
                            fig, (ax0,ax1,ax2,ax3) = plt.subplots(4,1, figsize = (20,20))
                            
                            
                            legd_kw = {'fontsize' : 12, 'loc' : 'upper left'}
                            alpha = 0.4
                            ## orographic panel (bottom, ax4)
                            # =============================================================================
                            #                         ax4.plot(orogcoor,orogdata, label = 'Model background orography', color = 'k')
                            #                         ax4.plot(cube_coor, surf_alt, linestyle = '--', label = 'Model orography directly below', color = 'k')
                            #                         ax4.plot(cube_coor, uni_alt, color = 'g', label = 'Interpolation surface')
                            #                         ax4.plot(db_coor, db_alt, color = 'g', linestyle = '--', label = 'Obs gps')
                            #                         ax4.set_ylabel('m')
                            #                         ax4.set_ylim(bottom = 0, top = 1700)
                            #                         ax4.legend()
                            #                         ax4.set_xlim(left = xmin, right = xmax)
                            #                         ax4.set_xticks(tick_coors)
                            #                         ax4.set_xticklabels(tick_dists)
                            #                         ax4.set_xlabel('Scalar surface distance [km]')
                            #                         ax4.grid()
                            # =============================================================================
                            ## wsp
                            #ax0.plot(simple_coors_wsp, simple_wsp, color = 'blue', label = 'UM wsp')
                            #ax0.fill_between(simple_coors_wsp,simple_wsp,simple_wsp_m1h, alpha = alpha, color = 'deepskyblue', label = 'UM -1h')
                            #ax0.fill_between(simple_coors_wsp,simple_wsp,simple_wsp_p1h, alpha = alpha, color = 'coral', label = 'UM +1h')
                            
                            ax0.plot(matched_coors,matched_wsp, color = 'red', label = r'UM matched')
                            ax0.fill_between(matched_coors,matched_wsp,matched_wsp_m1h, alpha = alpha, color = 'peru', label = 'UM matched -1h')
                            ax0.fill_between(matched_coors,matched_wsp,matched_wsp_p1h, alpha = alpha, color = 'blueviolet', label = 'UM matched +1h')
                            ax0.plot(db_coor,db_wsp, label = 'obs', color = 'black')
                            ax0.set_ylabel(w_cube.units)
                            ax0.set_ylim(bottom = 0, top = 30)
                            #ax0.set_xlim(left = xmin, right = xmax)
                            ax0.set_xticks(tick_coors)
                            ax0.set_xticklabels([])
                            ax0.set_title(f'{res} Atmospheric state vars, leg {leg} xsection, {dates[tidx]}, {(np.nanmean(db_alt)//1)}m altitude')
                            ax0.legend(**legd_kw)
                            ax0.grid()
                            ## w
                            #ax1.plot(simple_coors_w, simple_w, color = 'blue', label = 'UM w')
                            
                            #ax1.fill_between(simple_coors_w,simple_w,simple_w_m1h, alpha = alpha, color = 'deepskyblue', label = 'UM -1h')
                            #ax1.fill_between(simple_coors_w,simple_w,simple_w_p1h, alpha = alpha, color = 'coral', label = 'UM +1h')
                            #ax1.plot(db_coor,db_wsp, label = 'obs', color = 'black')
                            ax1.plot(matched_coors,matched_w, color = 'red', label = r'UM $w$ matched')
                            
                            ax1.fill_between(matched_coors,matched_w,matched_w_m1h, alpha = alpha, color = 'peru', label = 'UM matched -1h')
                            ax1.fill_between(matched_coors,matched_w,matched_w_p1h, alpha = alpha, color = 'blueviolet', label = 'UM matched +1h')
                            ax1.plot(db_coor,db_w, label = 'obs', color = 'black')
                            ax1.set_ylabel(w_cube.units)
                            ax1.set_ylim(bottom = -5, top = 5)
                            #ax1.set_xlim(left = xmin, right = xmax)
                            ax1.set_xticks(tick_coors)
                            ax1.set_xticklabels([])
                            #ax1.legend(**legd_kw)
                            ax1.grid()
                            ## theta
                            #ax2.plot(simple_coors_t, simple_theta, color = 'blue', label = 'UM theta')
                            
                            #ax2.fill_between(simple_coors_t,simple_theta,simple_theta_m1h, alpha = alpha, color = 'deepskyblue', label = 'UM -1h')
                            #ax2.fill_between(simple_coors_t,simple_theta,simple_theta_p1h, alpha = alpha, color = 'coral', label = 'UM +1h')
                            
                            ax2.plot(matched_coors,matched_theta, color = 'red', label = r'UM $\theta$ matched')
                            
                            ax2.fill_between(matched_coors,matched_theta,matched_theta_m1h, alpha = alpha, color = 'peru', label = 'UM matched -1h')
                            ax2.fill_between(matched_coors,matched_theta,matched_theta_p1h, alpha = alpha, color = 'blueviolet', label = 'UM matched +1h')
                            ax2.plot(db_coor,db_theta, label = 'obs', color = 'black')
                            ax2.set_ylabel(theta_cube.units)
                            ax2.set_ylim(bottom = 270, top = 290)
                            #ax2.set_xlim(left = xmin, right = xmax)
                            ax2.set_xticks(tick_coors)
                            ax2.set_xticklabels([])
                            #ax2.legend(**legd_kw)
                            ax2.grid()
                            ## specfific humidity
                            #ax3.plot(simple_coors_q, simple_q*1000, color = 'blue', label = 'UM q')
                            
                            #ax3.fill_between(simple_coors_q,simple_q*1000,simple_q_m1h*1000, alpha = alpha, color = 'deepskyblue', label = 'UM -1h')
                            #ax3.fill_between(simple_coors_q,simple_q*1000,simple_q_p1h*1000, alpha = alpha, color = 'coral', label = 'UM +1h')
                            #ax3.plot(db_coor,db_q, label = 'obs', color = 'red')
                            ax3.plot(matched_coors,matched_q*1000, color = 'red', label = r'UM $q$ matched')
                            
                            ax3.fill_between(matched_coors,matched_q*1000,matched_q_m1h*1000, alpha = alpha, color = 'peru', label = 'UM matched -1h')
                            ax3.fill_between(matched_coors,matched_q*1000,matched_q_p1h*1000, alpha = alpha, color = 'blueviolet', label = 'UM matched +1h')
                            ax3.plot(db_coor,db_q, label = 'obs', color = 'black')
                            ax3.set_ylabel(r'$gkg^{-1}$')
                            ax3.set_ylim(bottom = 0, top = 5)
                            #ax3.set_xlim(left = xmin, right = xmax)
                            #ax3.legend(**legd_kw)
                            
                            ax3.set_xticks(tick_coors)
                            ax3.set_xticklabels(tick_dists)
                            ax3.set_xlabel('Scalar surface distance [km]')
                            ax3.grid()
                            
                            plt.subplots_adjust(wspace=0.06, hspace=0.07) 
                            
                            save_fig = True
                            if save_fig:
                                timepoint = dates[tidx].strftime('H%HM%M')
                                exp = 'control'
                                fig.savefig(f'D:/Project/Figures/PNG/306/{suite}/spatial_series/60s/{leg}/{config}_{exp}_{res}_atmostate_series_60s_matchonly_{timepoint}_leg{leg}_306.png')
                                fig.savefig(f'D:/Project/Figures/PDF/306/{suite}/spatial_series/60s/{leg}/{config}_{exp}_{res}_atmostate_series_60s_matchonly_{timepoint}_leg{leg}_306.pdf')
                            plt.show()
                            #%%
                        except:
                            print('Exception occured, passing.')
                            pass

                                                
spatial_60s_fluxes = True
if spatial_60s_fluxes:
    for res in resolutions:
        print(f'\n\nRes {res}')
        #%% load model data
        ### 3D whole cubes
        #res = '0p5km'
        alt_idx = 35
        ## windstress (two cumponents, need to be combined)
        # wse = windstrss eastward, wsn = windstress northward
        wse_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_{res}_um_atmosphere_downward_eastward_stress_24hrs_ph_306.nc', 'atmosphere_downward_eastward_stress')[8:-2,:alt_idx]
        wsn_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_{res}_um_atmosphere_downward_northward_stress_24hrs_ph_306.nc','atmosphere_downward_northward_stress')[8:-2,:alt_idx]
        # regrid wsn cube and combine cubes for magnitude
        #%%
        wsn_cube_reg = wsn_cube.regrid(wse_cube, iris.analysis.Linear())
        ws_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_{res}_um_atmosphere_downward_eastward_stress_24hrs_ph_306.nc', 'atmosphere_downward_eastward_stress')[8:-2,:alt_idx]
        ws_cube.data = (wse_cube.data**2 + wsn_cube_reg.data**2)**0.5
        #%%
        ## load heatfluxes
        # lh = latent heat, sh = sensible heat
        #lh_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_{res}_um_surface_upward_latent_heat_flux_24hrs_pg_306.nc', 'atmosphere_upward_latent_heat_flux')[:,:alt_idx]
        sh_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_{res}_um_upward_heat_flux_in_air_24hrs_pi_306.nc', 'upward_heat_flux_in_air')[8:-2,:alt_idx]
        # calculate lh
        #%%
        theta_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{config}_{res}_um_air_potential_temperature_24hrs_pi_306.nc', 'air_potential_temperature')[8:-2,:alt_idx]
        p_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{config}_{res}_um_air_pressure_24hrs_pi_306.nc', 'air_pressure')[8:-2,:alt_idx]
        q_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{config}_{res}_um_specific_humidity_24hrs_pi_306.nc', 'specific_humidity')[8:-2,:alt_idx]
        vf_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{config}_{res}_um_upward_water_vapor_flux_in_air_24hrs_pi_306.nc', 'upward_water_vapor_flux_in_air')[8:-2,:alt_idx] # vapor flux
        lh_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_{res}_um_upward_heat_flux_in_air_24hrs_pi_306.nc', 'upward_heat_flux_in_air')[8:-2,:alt_idx]
        lh_data = dmna.calc_lhf(q_cube.data,theta_cube.data,p_cube.data,vf_cube.data)
        lh_cube.data = lh_data
        #%% load and process database shortrun obs data
        # using pandas
        database = pd.read_csv('D:/Project/Obvs_Data/Databases/IGP_flights_database_obs_60s_306.txt', delimiter = ' ')
        
        db_legno = database['legno']
        
        dblegs = [1,2,4,5,6,7,8,9,13,15]
        for leg in dblegs:
            print(f'Leg {leg}')
            legcond = db_legno == leg
            left = 0; right = -1
            db_ws = np.array(database['windstress'][legcond])[left:right] 
            db_lh = np.array(database['lh'][legcond])[left:right]
            db_sh = np.array(database['sh'][legcond])[left:right]
            
            
            db_lon = np.array(database['lon'][legcond])[left:right]
            db_lat = np.array(database['lat'][legcond])[left:right]
            db_alt = np.array(database['altgps'][legcond])[left:right]
            #%%
            print(len(db_alt), len(db_lon), len(db_ws))
            #%% rotate db coords
        
            polelat = wse_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
            polelon = wse_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
            rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
            rot_db_lon = rot_db_lon + 360
            
            #%% format date time and comute mean leg time
            time = wse_cube.coord('time')
            dates = time.units.num2date(time.points)
            db_time = np.array(database['meantime'][legcond])
            mid_point = np.nanmean(db_time)
            model_time_ssm = (time.points-time.points[0]+4.5)*60*60 # get seconds since midnight from model output timesteps
            model_time_hsm = (time.points-time.points[0]+4.5) # hours since midnight
            print(model_time_hsm)
            #%%
            for tidx in range(len(model_time_hsm)):
                # go through each time step and check which is valid in next line
                diff = mid_point/(60*60) - model_time_hsm[tidx]
                if (diff < 0.25) and (diff >= -0.25):
                    # if meantime falls within 15min of model output
                    print('Time: obs', mid_point/(60*60), ', model', model_time_hsm[tidx]);
                    ## perform interpolations for current timestep
                    try:
                        #%%
                        
                        print('Interpolating index...')
                        time_index = tidx
                        print(time_index)
                        matched_ws, matched_coors = pert.interpolate_matched_series(ws_cube, db_alt, rot_db_lon, time_index, obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                        print(len(matched_coors), len(matched_ws))
                        matched_lh, matched_coors = pert.interpolate_matched_series(lh_cube, db_alt, rot_db_lon, time_index, obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                        
                        matched_sh, matched_coors = pert.interpolate_matched_series(sh_cube, db_alt, rot_db_lon, time_index,obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                        
                        db_coor = rot_db_lat
                        ## perform interpolations for prior timestep
                        print('Interpolating index -2...')
                        time_index = tidx - 2
                        print(time_index)
                        matched_ws_m1h, matched_coors = pert.interpolate_matched_series(ws_cube, db_alt, rot_db_lon, time_index, obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                        print(len(matched_coors), len(matched_ws_m1h))
                        matched_lh_m1h, matched_coors = pert.interpolate_matched_series(lh_cube, db_alt, rot_db_lon, time_index, obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                        
                        matched_sh_m1h, matched_coors = pert.interpolate_matched_series(sh_cube, db_alt, rot_db_lon, time_index,obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                        ## perform interpolations for prior timestep
                        print('Interpolating index +2...')
                        time_index = tidx + 2
                        print(time_index)
                        matched_ws_p1h, matched_coors = pert.interpolate_matched_series(ws_cube, db_alt, rot_db_lon, time_index, obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                        print(len(matched_coors), len(matched_ws_p1h))
                        matched_lh_p1h, matched_coors = pert.interpolate_matched_series(lh_cube, db_alt, rot_db_lon, time_index, obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                        
                        matched_sh_p1h, matched_coors = pert.interpolate_matched_series(sh_cube, db_alt, rot_db_lon, time_index,obs_lats = rot_db_lat, grid_axis = 'grid_latitude', mode = '3d')
                        
                        #%%
                        print(len(matched_coors), len(matched_ws), len(matched_ws_p1h))
                        #%% create uniform altitude array for indication on plot
                        #uni_alt = np.ones(len(cube_coor))*np.nanmean(db_alt)
                        #%% dynamically adjust horizontal bounds
                        xmin = np.nanmin(matched_coors)-0.1
                        xmax = np.nanmax(matched_coors)+0.1
                        
                        #%% get scalar distance ticklist for x axis labelling
                        tick_coors, tick_dists = workshop.obs_dist_tick_gen(rot_db_lon,rot_db_lat, 'grid_latitude', ticknumber = 6, decprecision = 2)
                        
                        #%% plot some stuff
                        print('Commence plotting')
                        fig, (ax0,ax1,ax2) = plt.subplots(3,1, figsize = (20,15))
                        
                        
                        legd_kw = {'fontsize' : 12, 'loc' : 'upper left'}
                        alpha = 0.4
                        
                        
                        ax0.plot(matched_coors,matched_ws, color = 'red', label = r'UM matched')
                        ax0.fill_between(matched_coors,matched_ws,matched_ws_m1h, alpha = alpha, color = 'peru', label = 'UM matched -1h')
                        ax0.fill_between(matched_coors,matched_ws,matched_ws_p1h, alpha = alpha, color = 'blueviolet', label = 'UM matched +1h')
                        ax0.plot(db_coor,db_ws, label = 'obs', color = 'black')
                        ax0.set_ylabel(f'WS {ws_cube.units}')
                        #ax0.set_ylim(bottom = 0, top = 30)
                        #ax0.set_xlim(left = xmin, right = xmax)
                        ax0.set_xticks(tick_coors)
                        ax0.set_xticklabels([])
                        
                        ax0.legend(**legd_kw)
                        ax0.grid()
                        ## w
                        #ax1.plot(simple_coors_w, simple_w, color = 'blue', label = 'UM w')
                        
                        #ax1.fill_between(simple_coors_w,simple_w,simple_w_m1h, alpha = alpha, color = 'deepskyblue', label = 'UM -1h')
                        #ax1.fill_between(simple_coors_w,simple_w,simple_w_p1h, alpha = alpha, color = 'coral', label = 'UM +1h')
                        #ax1.plot(db_coor,db_wsp, label = 'obs', color = 'black')
                        ax1.plot(matched_coors,matched_sh, color = 'red', label = r'UM $sh$ matched')
                        
                        ax1.fill_between(matched_coors,matched_sh,matched_sh_m1h, alpha = alpha, color = 'peru', label = 'UM matched -1h')
                        ax1.fill_between(matched_coors,matched_sh,matched_sh_p1h, alpha = alpha, color = 'blueviolet', label = 'UM matched +1h')
                        ax1.plot(db_coor,db_sh, label = 'obs', color = 'black')
                        ax1.set_ylabel(f'SH {sh_cube.units}')
                        #ax1.set_ylim(bottom = -5, top = 5)
                        #ax1.set_xlim(left = xmin, right = xmax)
                        ax1.set_xticks(tick_coors)
                        ax1.set_xticklabels([])
                        #ax1.legend(**legd_kw)
                        ax1.grid()
                        ## theta
                        #ax2.plot(simple_coors_t, simple_theta, color = 'blue', label = 'UM theta')
                        
                        #ax2.fill_between(simple_coors_t,simple_theta,simple_theta_m1h, alpha = alpha, color = 'deepskyblue', label = 'UM -1h')
                        #ax2.fill_between(simple_coors_t,simple_theta,simple_theta_p1h, alpha = alpha, color = 'coral', label = 'UM +1h')
                        
                        ax2.plot(matched_coors,matched_lh, color = 'red', label = r'UM $\theta$ matched')
                        
                        ax2.fill_between(matched_coors,matched_lh,matched_lh_m1h, alpha = alpha, color = 'peru', label = 'UM matched -1h')
                        ax2.fill_between(matched_coors,matched_lh,matched_lh_p1h, alpha = alpha, color = 'blueviolet', label = 'UM matched +1h')
                        ax2.plot(db_coor,db_lh, label = 'obs', color = 'black')
                        ax2.set_ylabel(f'LH {lh_cube.units}')
                        #ax2.set_ylim(bottom = 270, top = 290)
                        #ax2.set_xlim(left = xmin, right = xmax)
                        
                        #ax2.legend(**legd_kw)
                        ax2.set_xticks(tick_coors)
                        ax2.set_xticklabels(tick_dists)
                        ax2.set_xlabel('Scalar surface distance [km]')
                        ax2.grid()
                        
                        
                        plt.subplots_adjust(wspace=0.06, hspace=0.07) 
                        
                        save_fig = True
                        if save_fig:
                            timepoint = dates[tidx].strftime('H%HM%M')
                            exp = 'control'
                            fig.savefig(f'D:/Project/Figures/PDF/{flight}/{suite}/spatial_series/60s/{leg}/{config}_{exp}_{res}_fluxes_wsshlh_series_60s_matchonly_{timepoint}_leg{leg}_{flight}.pdf')
                            ax0.set_title(f'{res} Turbulent flux vars, leg {leg} xsection, {dates[tidx]}, {(np.nanmean(db_alt)//1)}m altitude')
                            fig.savefig(f'D:/Project/Figures/PNG/{flight}/{suite}/spatial_series/60s/{leg}/{config}_{exp}_{res}_fluxes_wsshlh_series_60s_matchonly_{timepoint}_leg{leg}_{flight}.png')
                            
                        plt.show()
                        #%%
                    except:
                        print('Exception occured, passing.')
                        pass