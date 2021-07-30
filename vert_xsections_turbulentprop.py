#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 17:38:28 2020

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
 #%% load and process core model diagnostics from vertical sliced files
obs_plotting = False
if obs_plotting:
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
            heat_cube = dmna.load_model_heat_flux(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
            lh_cube = dmna.load_model_heat_flux(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
            #tke_cube = dmna.load_model_tke(res,'Vertical_slices', fileleg[2],301, 'Control')
            print(vapor_cube)
            #print(tke_cube)
            
            
            #%% extract coordinates
            
            model_coors = p_cube.coord(fileleg[3][0])
            #%% load and process orography
            
            orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
            orogmax = orogcube.collapsed(fileleg[3][1], iris.analysis.MAX)
            orogdata = orogmax.data
            orogcoor = orogmax.coord(fileleg[3][0]).points
            
            #%%
            print(np.shape(theta_cube.data))
            print(np.shape(q_cube.data))
            print(np.shape(p_cube.data))
            print(np.shape(vapor_cube.data))
            print(np.shape(heat_cube.data))
            #%% calculate latent heat flux
            lh_cube.data = dmna.calc_lhf(q_cube.data, theta_cube.data, p_cube.data, vapor_cube.data)
            print(lh_cube)
            #%% give units
            #print(q_cube.units)
            # #%% test plotting
            # q = iplt.contourf(lh_cube[30])
            # plt.colorbar(q)
            # plt.show()

            #%% load and process databse shortrun obs data
            # using pandas
            database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_sr_301.txt', delimiter = ' ')
            leg = fileleg[0]
            print(database)
            #%%
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
            db_wst = np.array(database['windstress'][legcond]) 
            db_tke = np.array(database['tke'][legcond])
            db_sh = np.array(database['sh'][legcond])
            db_lh = np.array(database['lh'][legcond])
            
            db_lon = np.array(database['lon'][legcond])
            db_lat = np.array(database['lat'][legcond])
            db_alt = np.array(database['altgps'][legcond])
            
            #%% rotate db coords
            compcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/RA1M_0p5km_umi1_flt{flight}.nc', varname_p)[0,0]
            polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
            polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
            rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
            rot_db_lon = rot_db_lon + 360
            if fileleg[3][0] == 'grid_longitude':
                db_coor = rot_db_lon
            else:
                db_coor = rot_db_lat
            #%% get scalar distance ticklist for x axis labelling
            tick_coors, tick_dists = workshop.dist_tick_generator(p_cube, fileleg[3][0], ticknumber = 6, decprecision = 2)
            #%% format date time for file name and figure title
            time = p_cube.coord('time')
            dates = time.units.num2date(time.points)
            print(dates[fileleg[1]])
            
            
            #%% set colorbar normalisations
            
            wstmin = 0; wstmax = 2.0
            tkemin = 0; tkemax = 1.5
            shmin = -200; shmax = -shmin
            lhmin = -400; lhmax = - lhmin
            wst_norm = matplotlib.colors.Normalize(vmin = wstmin, vmax = wstmax)
            tke_norm = matplotlib.colors.Normalize(vmin = tkemin, vmax = tkemax)
            sh_norm = matplotlib.colors.Normalize(vmin = shmin, vmax = shmax)
            lh_norm = matplotlib.colors.Normalize(vmin = lhmin, vmax = lhmax)
            
            wst_map = matplotlib.cm.ScalarMappable(norm = wst_norm, cmap = 'Blues')
            tke_map = matplotlib.cm.ScalarMappable(norm = tke_norm, cmap = 'Greens')
            sh_map = matplotlib.cm.ScalarMappable(norm = sh_norm, cmap = plt.cm.get_cmap('BrBG').reversed())
            lh_map = matplotlib.cm.ScalarMappable(norm = lh_norm, cmap = 'seismic')
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
            ax0.set_title('windstress')
            ax0.set_xticks([])
            ax0.set_ylabel('m')
            ax0.set_xlim(left = xmin, right = xmax)
            pwst = iplt.contourf(wst_cube[tidx], **common_keywargs, cmap = 'Blues', vmin = wstmin, vmax  = wstmax)
            plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
            wst_cbar = plt.colorbar(wst_map, ax = ax0, **com_cbarkeys)
            wst_cbar.ax.set(ylabel = r'$Nm^{-1}$')
            ## wsp obs
            ax0.scatter(db_coor, db_alt, c = db_wst, cmap = 'Blues', vmin = wstmin, vmax = wstmax,
                        s = 300, edgecolor = 'k')
            ## w
            ax1 = fig.add_subplot(gs[1,0])
            psh = iplt.contourf(heat_cube[tidx], **common_keywargs, cmap = plt.cm.get_cmap('BrBG').reversed(), vmin = shmin, vmax = shmax)
            ax1.set_ylim(top = 2000)
            ax1.set_title('upward sensible heat flux')
            ax1.set_xticks([])
            ax1.set_ylabel('m')
            ax1.set_xlim(left = xmin, right = xmax)
            plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
            sh_cbar = plt.colorbar(sh_map, ax = ax1, **com_cbarkeys)
            sh_cbar.ax.set(ylabel = r'$Wm^{-2}$')
            ## w obs
            ax1.scatter(db_coor, db_alt, c = db_sh, cmap = plt.cm.get_cmap('BrBG').reversed(), vmin = shmin, vmax = shmax,
                        s = 300, edgecolor = 'k')
            ## theta
            
            ax2 = fig.add_subplot(gs[2,0])
            ptheta = iplt.contourf(lh_cube[tidx], **common_keywargs,cmap = 'seismic', vmin = lhmin, vmax = lhmax)
            ax2.set_ylim(top = 2000)
            ax2.set_title('upward latent heat flux')
            ax2.set_xticks([])
            ax2.set_ylabel('m')
            ax2.set_xlim(left = xmin, right = xmax)
            plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
            lh_cbar = plt.colorbar(lh_map, ax = ax2, **com_cbarkeys)
            lh_cbar.ax.set(ylabel = r'$Wm^{-2}$')
            ## theta obs
            ax2.scatter(db_coor, db_alt, c = db_lh, cmap = 'seismic', vmin = lhmin, vmax = lhmax,
                        s = 300, edgecolor = 'k')
            # ## q
            # ax3 = fig.add_subplot(gs[3,0])
            # pq = iplt.contourf(tke_cube[tidx], coords = [geo_axis, 'level_height'], levels =35, cmap = 'Greens', vmin = tkemin, vmax = tkemax)
            # ax3.set_ylim(top = 2000)
            # ax3.set_title('turbulent kinetic energy')
            ax2.set_xticks(tick_coors)
            ax2.set_xticklabels(tick_dists)
            ax2.set_xlabel('Scalar distance, km')
            # ax3.set_ylabel('m')
            # ax3.set_xlim(left = xmin, right = xmax)
            # plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
            # tke_cbar = plt.colorbar(tke_map, ax = ax3, **com_cbarkeys)
            # tke_cbar.ax.set(ylabel = '$m^{-2}s^{-2}$')
            # ## q obs
            # ax3.scatter(db_coor, db_alt, c = db_tke, cmap = 'Greens', vmin = tkemin, vmax = tkemax,
            #             s = 300, edgecolor = 'k')
            fig.suptitle(f'Turbulent quantities, {res}, leg {fileleg[0]} xsection, {dates[fileleg[1]]}')
            fig.tight_layout()
            plt.show()
           
            save_fig = True
            if save_fig:
                timepoint = dates[tidx].strftime('H%HM%M')
                fig.savefig(f'../../Figures/PNG/301/{suite}/xsections_vert/{exp}/turbprops/{res}/obs/{config}_{res}_turbprops_vert_xsection_with_obs_{exp}_{timepoint}_leg{fileleg[0]}.png')
            #sys.exit()

#%%
#leg_meta_time = ((8,31,8, lat_lon_ord),(0,0))



time_plotting = True
if time_plotting:
    for fileleg in leg_meta_time:
        for res in resolutions:
            for tidx in range(48):
                plt.close()
                wst_cube = dmna.load_model_windstress(res, varname_wse, varname_wsn, 'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
                p_cube = dmna.load_model_pressure(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
                q_cube = dmna.load_model_specific_hum(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
                theta_cube = dmna.load_model_theta(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35, secondary_var = 'upward_air_velocity')
                vapor_cube = dmna.load_model_vapor_flux(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
                heat_cube = dmna.load_model_heat_flux(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
                lh_cube = dmna.load_model_heat_flux(res,'Vertical_slices', fileleg[2],301, 'Control', level_range = 35)
                #tke_cube = dmna.load_model_tke(res,'Vertical_slices', fileleg[2],301, 'Control')
                #%%
                #%% extract coordinates
            
                model_coors = p_cube.coord(fileleg[3][0])
                #%% load and process orography
                
                orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
                orogmax = orogcube.collapsed(fileleg[3][1], iris.analysis.MAX)
                orogdata = orogmax.data
                orogcoor = orogmax.coord(fileleg[3][0]).points
                
                
                #%% calculate latent heat flux
                lh_cube.data = dmna.calc_lhf(q_cube.data, theta_cube.data, p_cube.data, vapor_cube.data)
                
                #%% get scalar distance ticklist for x axis labelling
                tick_coors, tick_dists = workshop.dist_tick_generator(p_cube, fileleg[3][0], ticknumber = 6, decprecision = 2)
                #%% format date time for file name and figure title
                time = p_cube.coord('time')
                dates = time.units.num2date(time.points)
                print(dates[tidx])
                
                
                #%% set colorbar normalisations
                
                wstmin = 0; wstmax = 2.0
                tkemin = 0; tkemax = 1.5
                shmin = -200; shmax = -shmin
                lhmin = -400; lhmax = - lhmin
                wst_norm = matplotlib.colors.Normalize(vmin = wstmin, vmax = wstmax)
                tke_norm = matplotlib.colors.Normalize(vmin = tkemin, vmax = tkemax)
                sh_norm = matplotlib.colors.Normalize(vmin = shmin, vmax = shmax)
                lh_norm = matplotlib.colors.Normalize(vmin = lhmin, vmax = lhmax)
                
                wst_map = matplotlib.cm.ScalarMappable(norm = wst_norm, cmap = 'Blues')
                tke_map = matplotlib.cm.ScalarMappable(norm = tke_norm, cmap = 'Greens')
                sh_map = matplotlib.cm.ScalarMappable(norm = sh_norm, cmap = plt.cm.get_cmap('BrBG').reversed())
                lh_map = matplotlib.cm.ScalarMappable(norm = lh_norm, cmap = 'seismic')
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
                ax0.set_title('windstress')
                ax0.set_xticks([])
                ax0.set_ylabel('m')
                ax0.set_xlim(left = xmin, right = xmax)
                pwst = iplt.contourf(wst_cube[tidx], **common_keywargs, cmap = 'Blues', vmin = wstmin, vmax  = wstmax)
                plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                wst_cbar = plt.colorbar(wst_map, ax = ax0, **com_cbarkeys)
                wst_cbar.ax.set(ylabel = r'$Nm^{-1}$')
                ## wsp obs
                
                ## w
                ax1 = fig.add_subplot(gs[1,0])
                psh = iplt.contourf(heat_cube[tidx], **common_keywargs, cmap = plt.cm.get_cmap('BrBG').reversed(), vmin = shmin, vmax = shmax)
                ax1.set_ylim(top = 2000)
                ax1.set_title('upward sensible heat flux')
                ax1.set_xticks([])
                ax1.set_ylabel('m')
                ax1.set_xlim(left = xmin, right = xmax)
                plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                sh_cbar = plt.colorbar(sh_map, ax = ax1, **com_cbarkeys)
                sh_cbar.ax.set(ylabel = r'$Wm^{-2}$')
                ## w obs
                
                ## theta
                
                ax2 = fig.add_subplot(gs[2,0])
                ptheta = iplt.contourf(lh_cube[tidx], **common_keywargs,cmap = 'seismic', vmin = lhmin, vmax = lhmax)
                ax2.set_ylim(top = 2000)
                ax2.set_title('upward latent heat flux')
                ax2.set_xticks([])
                ax2.set_ylabel('m')
                ax2.set_xlim(left = xmin, right = xmax)
                plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                lh_cbar = plt.colorbar(lh_map, ax = ax2, **com_cbarkeys)
                lh_cbar.ax.set(ylabel = r'$Wm^{-2}$')
                ## theta obs
                
                # ## q
                # ax3 = fig.add_subplot(gs[3,0])
                # pq = iplt.contourf(tke_cube[tidx], coords = [geo_axis, 'level_height'], levels =35, cmap = 'Greens', vmin = tkemin, vmax = tkemax)
                # ax3.set_ylim(top = 2000)
                # ax3.set_title('turbulent kinetic energy')
                ax2.set_xticks(tick_coors)
                ax2.set_xticklabels(tick_dists)
                ax2.set_xlabel('Scalar distance, km')
                # ax3.set_ylabel('m')
                # ax3.set_xlim(left = xmin, right = xmax)
                # plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                # tke_cbar = plt.colorbar(tke_map, ax = ax3, **com_cbarkeys)
                # tke_cbar.ax.set(ylabel = '$m^{-2}s^{-2}$')
                # ## q obs
                # ax3.scatter(db_coor, db_alt, c = db_tke, cmap = 'Greens', vmin = tkemin, vmax = tkemax,
                #             s = 300, edgecolor = 'k')
                fig.suptitle(f'Turbulent quantities, {res}, leg {fileleg[0]} xsection, {dates[tidx]}')
                fig.tight_layout()
                plt.show()
               
                save_fig = True
                if save_fig:
                    timepoint = dates[tidx].strftime('H%HM%M')
                    fig.savefig(f'../../Figures/PNG/301/{suite}/xsections_vert/{exp}/turbprops/{res}/{fileleg[0]}/{config}_{res}_turbprops_vert_xsection_with_obs_{exp}_{timepoint}_leg{fileleg[0]}.png')
                #sys.exit()
                