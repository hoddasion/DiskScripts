#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 16:53:40 2020

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
#%%
plt.rcParams['font.size'] = 26
resolutions = ['0p5km', '1p5km', '4p4km']
varname_u = 'm01s15i002'
varname_v = 'm01s15i003'
varname_wsp = 'windspeed'
varname_w = 'upward_air_velocity'
varname_theta = 'air_potential_temperature'
varname_q = 'specific_humidity'

lat_lon_ord = ('grid_latitude', 'grid_longitude')
lon_lat_ord = ('grid_longitude', 'grid_latitude')
leg_meta_obs = ((1,26,1, lon_lat_ord),(2,27,2,lon_lat_ord),(3,28,3,lon_lat_ord),(4,29,3,lon_lat_ord),
            (5,30,3, lon_lat_ord),(6,30,3, lon_lat_ord),(8,31,8, lat_lon_ord),(9,32,8, lat_lon_ord))
leg_meta_time = ((1,26,1, lon_lat_ord),(2,27,2,lon_lat_ord),(3,28,3,lon_lat_ord),(8,31,8, lat_lon_ord))
flight = 301
config = 'RA1M'
exp = 'Control'
suite = 'u-bu807'
 #%% load and process core model diagnostics from vertical sliced files
obs_plotting = True
if obs_plotting:
    for fileleg in leg_meta_obs:
        for res in resolutions:
            #%%
            #fileleg = (1,26,1, lon_lat_ord)
            plt.close()
            wsp_cube = dmna.load_model_wsp(res, varname_u, varname_v, 'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
            w_cube = dmna.load_model_w(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
            q_cube = dmna.load_model_specific_hum(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
            theta_cube = dmna.load_model_theta(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35, secondary_var = 'upward_air_velocity')
            
            #%% extract coordinates
            
            model_coors = w_cube.coord(fileleg[3][0])
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
            if leg == 3:
                print('3456 active')
                legcond = np.array(db_legno == 3) + np.array(db_legno == 4) + np.array(db_legno == 5) + np.array(db_legno == 6)
            elif leg == 4:
                legcond = np.array(db_legno == 3) + np.array(db_legno == 4) + np.array(db_legno == 5) + np.array(db_legno == 6)
            elif leg == 5:
                legcond = np.array(db_legno == 3) + np.array(db_legno == 4) + np.array(db_legno == 5) + np.array(db_legno == 6)
            elif leg == 6:
                legcond = np.array(db_legno == 3) + np.array(db_legno == 4) + np.array(db_legno == 5) + np.array(db_legno == 6)
            elif fileleg[0] == 8:
                print('89 active')
                legcond = np.array(db_legno == 8) + np.array(db_legno == 9)
            elif leg == 9:
                legcond = np.array(db_legno == 8) + np.array(db_legno == 9)
            else:
                print('12 active')
                legcond = db_legno == fileleg[0]
            db_wsp = np.array(database['wsp'][legcond]) 
            db_w = np.array(database['w'][legcond])
            db_q = np.array(database['q_BUCK'][legcond])
            db_theta = np.array(database['theta'][legcond])
            
            db_lon = np.array(database['lon'][legcond])
            db_lat = np.array(database['lat'][legcond])
            db_alt = np.array(database['altgps'][legcond])
            
            #%% rotate db coords
            compcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/RA1M_0p5km_umh1_flt{flight}.nc', varname_u)[0,0]
            polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
            polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
            rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
            rot_db_lon = rot_db_lon + 360
            if fileleg[3][0] == 'grid_longitude':
                db_coor = rot_db_lon
            else:
                db_coor = rot_db_lat
            #%% get scalar distance ticklist for x axis labelling
            tick_coors, tick_dists = workshop.dist_tick_generator(w_cube, fileleg[3][0], ticknumber = 6, decprecision = 2)
            #%% format date time for file name and figure title
            time = w_cube.coord('time')
            dates = time.units.num2date(time.points)
            print(dates[fileleg[1]])
            
            
            #%% set colorbar normalisations
            
            wspmin = 0; wspmax = 18
            wmin = -5; wmax = -wmin
            Tmin = 270; Tmax = 285
            qmin = 0; qmax = 3
            wsp_norm = matplotlib.colors.Normalize(vmin = wspmin, vmax = wspmax)
            w_norm = matplotlib.colors.Normalize(vmin = wmin, vmax = wmax)
            theta_norm = matplotlib.colors.Normalize(vmin = Tmin, vmax = Tmax)
            q_norm = matplotlib.colors.Normalize(vmin = qmin, vmax = qmax)
            
            wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = 'Oranges')
            w_map = matplotlib.cm.ScalarMappable(norm = w_norm, cmap = 'seismic')
            theta_map = matplotlib.cm.ScalarMappable(norm = theta_norm, cmap = 'plasma')
            q_map = matplotlib.cm.ScalarMappable(norm = q_norm, cmap = 'Blues')
            #%% dynamically adjust horizontal bounds
            xmin = np.nanmin(model_coors.points)
            xmax = np.nanmax(model_coors.points)
            #%% commence plotting
            tidx = fileleg[1]
            geo_axis = fileleg[3][0]
            common_keywargs = {'coords' : [geo_axis,'altitude'], 'levels' : 35}
            com_cbarkeys = {'pad' : 0.005}
            fig = plt.figure(figsize = (24,20))
            gs = gridspec.GridSpec(nrows = 4, ncols = 1)
            
            ## wsp
            
            ax0 = fig.add_subplot(gs[0,0])
            ax0.set_ylim(top = 2000)
            ax0.set_title('horizontal windspeed')
            ax0.set_xticks([])
            ax0.set_ylabel('m')
            ax0.set_xlim(left = xmin, right = xmax)
            pwsp = iplt.contourf(wsp_cube[tidx], **common_keywargs, cmap = 'Oranges', vmin = wspmin, vmax  = wspmax)
            plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
            wsp_cbar = plt.colorbar(wsp_map, ax = ax0, **com_cbarkeys)
            wsp_cbar.ax.set(ylabel = r'$ms^{-1}$')
            ## wsp obs
            ax0.scatter(db_coor, db_alt, c = db_wsp, cmap = 'Oranges', vmin = wspmin, vmax = wspmax,
                        s = 300, edgecolor = 'k')
            ## w
            ax1 = fig.add_subplot(gs[1,0])
            pw = iplt.contourf(w_cube[tidx], **common_keywargs, cmap = 'seismic', vmin = wmin, vmax = wmax)
            ax1.set_ylim(top = 2000)
            ax1.set_title('upward windspeed')
            ax1.set_xticks([])
            ax1.set_ylabel('m')
            ax1.set_xlim(left = xmin, right = xmax)
            plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
            w_cbar = plt.colorbar(w_map, ax = ax1, **com_cbarkeys)
            w_cbar.ax.set(ylabel = r'$ms^{-1}$')
            ## w obs
            ax1.scatter(db_coor, db_alt, c = db_w, cmap = 'seismic', vmin = wmin, vmax = wmax,
                        s = 300, edgecolor = 'k')
            ## theta
            
            ax2 = fig.add_subplot(gs[2,0])
            ptheta = iplt.contourf(theta_cube[tidx], **common_keywargs,cmap = 'plasma', vmin = Tmin, vmax = Tmax)
            ax2.set_ylim(top = 2000)
            ax2.set_title('potential temperature')
            ax2.set_xticks([])
            ax2.set_ylabel('m')
            ax2.set_xlim(left = xmin, right = xmax)
            plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
            theta_cbar = plt.colorbar(theta_map, ax = ax2, **com_cbarkeys)
            theta_cbar.ax.set(ylabel = 'K')
            ## theta obs
            ax2.scatter(db_coor, db_alt, c = db_theta, cmap = 'plasma', vmin = Tmin, vmax = Tmax,
                        s = 300, edgecolor = 'k')
            ## q
            ax3 = fig.add_subplot(gs[3,0])
            pq = iplt.contourf(q_cube[tidx], **common_keywargs, cmap = 'Blues', vmin = 0, vmax = 3)
            ax3.set_ylim(top = 2000)
            ax3.set_title('specific humidity')
            ax3.set_xticks(tick_coors)
            ax3.set_xticklabels(tick_dists)
            ax3.set_xlabel('Scalar distance, km')
            ax3.set_ylabel('m')
            ax3.set_xlim(left = xmin, right = xmax)
            plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
            q_cbar = plt.colorbar(q_map, ax = ax3, **com_cbarkeys)
            q_cbar.ax.set(ylabel = '$gkg^{-1}$')
            ## q obs
            ax3.scatter(db_coor, db_alt, c = db_q, cmap = 'Blues', vmin = qmin, vmax = qmax,
                        s = 300, edgecolor = 'k')
            fig.suptitle(f'Atmospheric state variables, {res}, leg {fileleg[0]} xsection, {dates[fileleg[1]]}')
            fig.tight_layout()
            save_fig = True
            if save_fig:
                timepoint = dates[tidx].strftime('H%HM%M')
                fig.savefig(f'../../Figures/PNG/301/{suite}/xsections_vert/{exp}/atmostate/{res}/obs/{config}_{res}_atmostate_vert_xsection_with_obs_{exp}_{timepoint}_leg{fileleg[0]}.png')

#%%
#leg_meta_time = ((8,31,8, lat_lon_ord),(0,0))



time_plotting = False
if time_plotting:
    for fileleg in leg_meta_time:
        for res in resolutions:
            for tidx in range(48):
                plt.close()
            
                wsp_cube = dmna.load_model_wsp(res, varname_u, varname_v, 'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
                w_cube = dmna.load_model_w(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
                q_cube = dmna.load_model_specific_hum(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
                theta_cube = dmna.load_model_theta(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35, secondary_var = 'upward_air_velocity')
                
                #%% extract coordinates
                
                model_coors = w_cube.coord(fileleg[3][0])
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
                db_wsp = np.array(database['wsp'][legcond]) 
                db_w = np.array(database['w'][legcond])
                db_q = np.array(database['q_BUCK'][legcond])
                db_theta = np.array(database['theta'][legcond])
                
                db_lon = np.array(database['lon'][legcond])
                db_lat = np.array(database['lat'][legcond])
                db_alt = np.array(database['altgps'][legcond])
                
                #%% rotate db coords
                compcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/RA1M_0p5km_umh1_flt{flight}.nc', varname_u)[0,0]
                polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
                polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
                rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
                rot_db_lon = rot_db_lon + 360
                
                #%% get scalar distance ticklist for x axis labelling
                tick_coors, tick_dists = workshop.dist_tick_generator(w_cube, fileleg[3][0], ticknumber = 6, decprecision = 2)
                #%% format date time for file name and figure title
                time = w_cube.coord('time')
                dates = time.units.num2date(time.points)
                print(dates[tidx])
                
                
                #%% set colorbar normalisations
                
                wspmin = 0; wspmax = 18
                wmin = -5; wmax = -wmin
                Tmin = 270; Tmax = 285
                qmin = 0; qmax = 3
                wsp_norm = matplotlib.colors.Normalize(vmin = wspmin, vmax = wspmax)
                w_norm = matplotlib.colors.Normalize(vmin = wmin, vmax = wmax)
                theta_norm = matplotlib.colors.Normalize(vmin = Tmin, vmax = Tmax)
                q_norm = matplotlib.colors.Normalize(vmin = qmin, vmax = qmax)
                
                wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = 'Oranges')
                w_map = matplotlib.cm.ScalarMappable(norm = w_norm, cmap = 'seismic')
                theta_map = matplotlib.cm.ScalarMappable(norm = theta_norm, cmap = 'plasma')
                q_map = matplotlib.cm.ScalarMappable(norm = q_norm, cmap = 'Blues')
                #%% dynamically adjust horizontal bounds
                xmin = np.nanmin(model_coors.points)
                xmax = np.nanmax(model_coors.points)
                #%% commence plotting
                
                geo_axis = fileleg[3][0]
                common_keywargs = {'coords' : [geo_axis,'altitude'], 'levels' : 35}
                com_cbarkeys = {'pad' : 0.005}
                fig = plt.figure(figsize = (24,20))
                gs = gridspec.GridSpec(nrows = 4, ncols = 1)
                
                ## wsp
                
                ax0 = fig.add_subplot(gs[0,0])
                ax0.set_ylim(top = 2000)
                ax0.set_title('horizontal windspeed')
                ax0.set_xticks([])
                ax0.set_ylabel('m')
                ax0.set_xlim(left = xmin, right = xmax)
                pwsp = iplt.contourf(wsp_cube[tidx], **common_keywargs, cmap = 'Oranges', vmin = wspmin, vmax  = wspmax)
                plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                wsp_cbar = plt.colorbar(wsp_map, ax = ax0, **com_cbarkeys)
                wsp_cbar.ax.set(ylabel = r'$ms^{-1}$')
                ## w
                ax1 = fig.add_subplot(gs[1,0])
                pw = iplt.contourf(w_cube[tidx], **common_keywargs, cmap = 'seismic', vmin = wmin, vmax = wmax)
                ax1.set_ylim(top = 2000)
                ax1.set_title('upward windspeed')
                ax1.set_xticks([])
                ax1.set_ylabel('m')
                ax1.set_xlim(left = xmin, right = xmax)
                plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                w_cbar = plt.colorbar(w_map, ax = ax1, **com_cbarkeys)
                w_cbar.ax.set(ylabel = r'$ms^{-1}$')
                ## theta
                
                ax2 = fig.add_subplot(gs[2,0])
                ptheta = iplt.contourf(theta_cube[tidx], **common_keywargs,cmap = 'plasma', vmin = Tmin, vmax = Tmax)
                ax2.set_ylim(top = 2000)
                ax2.set_title('potential temperature')
                ax2.set_xticks([])
                ax2.set_ylabel('m')
                ax2.set_xlim(left = xmin, right = xmax)
                plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                theta_cbar = plt.colorbar(theta_map, ax = ax2, **com_cbarkeys)
                theta_cbar.ax.set(ylabel = 'K')
                ## q
                ax3 = fig.add_subplot(gs[3,0])
                pq = iplt.contourf(q_cube[tidx], **common_keywargs, cmap = 'Blues', vmin = 0, vmax = 3)
                ax3.set_ylim(top = 2000)
                ax3.set_title('specific humidity')
                ax3.set_xticks(tick_coors)
                ax3.set_xticklabels(tick_dists)
                ax3.set_xlabel('Scalar distance, km')
                ax3.set_ylabel('m')
                ax3.set_xlim(left = xmin, right = xmax)
                plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                q_cbar = plt.colorbar(q_map, ax = ax3, **com_cbarkeys)
                q_cbar.ax.set(ylabel = '$gkg^{-1}$')
                
                fig.suptitle(f'Atmospheric state variables, {res}, leg {fileleg[0]} xsection, {dates[tidx]}')
                fig.tight_layout()
                save_fig = True
                if save_fig:
                    timepoint = dates[tidx].strftime('H%HM%M')
                    fig.savefig(f'../../Figures/PNG/301/{suite}/xsections_vert/{exp}/atmostate/{res}/leg{fileleg[0]}/{config}_{res}_atmostate_vert_xsection_{exp}_{timepoint}_leg{fileleg[0]}.png')