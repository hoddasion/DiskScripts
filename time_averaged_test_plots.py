# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 14:52:04 2022

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
import matplotlib.colors as colors
import warnings
from datetime import datetime


#warnings.filterwarnings("ignore")
#%%
plt.rcParams['font.size'] = 28
resolutions = ['0p5km']#,'1p5km', '4p4km']


lat_lon_ord = ('grid_latitude', 'grid_longitude')
lon_lat_ord = ('grid_longitude', 'grid_latitude')
config = 'LONGTAIL'
suite = 'u-cf117'
alt_idx = 35

#%% initiate resolution and flight leg  loops
for res in resolutions:
    for leg in [1,4,7,13,15]:
        ## load experiment data
        experiment = 'LONGTAIL'
        suite = 'u-cf117'
        config = 'LONGTAIL'
        factor_suffix = ''
        
        test_xwind_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/x_wind/{config}_vertslice_{res}_x_wind_flt306_leg{leg}{factor_suffix}.nc', 'x_wind')[:,:alt_idx]
        test_ywind_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/y_wind/{config}_vertslice_{res}_y_wind_flt306_leg{leg}{factor_suffix}.nc', 'y_wind')[:,:alt_idx]
        test_theta_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/air_potential_temperature/{config}_vertslice_{res}_air_potential_temperature_flt306_leg{leg}{factor_suffix}.nc', 'air_potential_temperature')[:,:alt_idx]
        test_q_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/specific_humidity/{config}_vertslice_{res}_specific_humidity_flt306_leg{leg}{factor_suffix}.nc', 'specific_humidity')[:,:alt_idx]
        test_w_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/upward_air_velocity/{config}_vertslice_{res}_upward_air_velocity_flt306_leg{leg}{factor_suffix}.nc', 'upward_air_velocity')[:,:alt_idx]
        ## compute wind magnitude
        wind_data = (test_xwind_cube.data**2 + test_ywind_cube.data**2)**0.5
        test_wind_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/x_wind/{config}_vertslice_{res}_x_wind_flt306_leg{leg}{factor_suffix}.nc', 'x_wind')[:,:alt_idx]
        test_wind_cube.data = wind_data
        ## load available flux vertslice cubes
        test_sh_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/upward_heat_flux_in_air/{config}_vertslice_{res}_upward_heat_flux_in_air_flt306_leg{leg}{factor_suffix}.nc', 'upward_heat_flux_in_air')[:,:alt_idx]
        test_lh_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/upward_latent_heat_flux_in_air/{config}_vertslice_{res}_upward_latent_heat_flux_in_air_flt306_leg{leg}{factor_suffix}.nc', 'upward_latent_heat_flux_in_air')[:,:alt_idx]
        test_wse_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/atmosphere_downward_eastward_stress/{config}_vertslice_{res}_atmosphere_downward_eastward_stress_flt306_leg{leg}{factor_suffix}.nc', 'atmosphere_downward_eastward_stress')[:,:alt_idx]
        test_wsn_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/atmosphere_downward_northward_stress/{config}_vertslice_{res}_atmosphere_downward_northward_stress_flt306_leg{leg}{factor_suffix}.nc', 'atmosphere_downward_northward_stress')[:,:alt_idx]
        ## compute stress magnitude
        test_ws_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/atmosphere_downward_eastward_stress/{config}_vertslice_{res}_atmosphere_downward_eastward_stress_flt306_leg{leg}{factor_suffix}.nc', 'atmosphere_downward_eastward_stress')[:,:alt_idx]
        ws_data = (test_wse_cube.data**2 + test_wsn_cube.data**2)**0.5
        test_ws_cube.data = ws_data
        print('Test data loaded')
        
        #%% format date time for file name and figure title
        time = test_w_cube.coord('time')
        dates = time.units.num2date(time.points)
        
        #%% get scalar distance ticklist for x axis labelling
        tick_coors, tick_dists = workshop.dist_tick_generator(test_w_cube, 'grid_latitude', ticknumber = 6, decprecision = 2)
        
        
        #%% load observations
        # using pandas
        database = pd.read_csv('D:/Project/Obvs_Data/Databases/IGP_flights_database_obs_60s_306.txt', delimiter = ' ')
        print(database.columns)
        #%%
        db_legno = database['legno']
        if leg == 1:
            legcond = (db_legno == 1) + (db_legno == 2)
        if leg == 4:
            legcond = (db_legno == 4) + (db_legno == 5) + (db_legno == 6) + (db_legno == 10)
        if leg == 7:
            legcond = (db_legno == 7) + (db_legno == 8) + (db_legno == 9)
        if leg == 13:
            legcond = db_legno == 13
        if leg == 15:
            legcond = db_legno == 15
        #legcond = db_legno == leg
        left = 0; right = -1
        db_wsp = np.array(database['wsp'][legcond])[left:right] 
        db_w = np.array(database['w'][legcond])[left:right]
        db_q = np.array(database['q_BUCK'][legcond])[left:right]
        db_theta = np.array(database['theta'][legcond])[left:right]
        db_ws = np.array(database['windstress'][legcond])[left:right] 
        db_sh = np.array(database['sh'][legcond])[left:right] 
        db_lh = np.array(database['lh'][legcond])[left:right]
        
        db_lon = np.array(database['lon'][legcond])[left:right]
        db_lat = np.array(database['lat'][legcond])[left:right]
        db_alt = np.array(database['altgps'][legcond])[left:right]
        db_start = np.array(database['starttime'][legcond])[0]/60/60
        
        print(db_start)
        #%% rotate db coords
        # first load 3d cube
        CUBE3D = iris.load_cube(f'D:/Project/Model_Data/{suite}/{config}_0p5km_um_air_potential_temperature_24hrs_pi_306.nc', 'air_potential_temperature')
        polelat = CUBE3D.coord('grid_latitude').coord_system.grid_north_pole_latitude
        polelon = CUBE3D.coord('grid_longitude').coord_system.grid_north_pole_longitude
        rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
        rot_db_lon = rot_db_lon + 360
        #%% start aggregate plots section
        aggregate_3hour_plots = True
        if aggregate_3hour_plots:
            ## aggregate cube data in 6 hour means
            for chunk in range(int(24/3)):
                ss = 6 # ss : sample size; for half hourly outputs, a sample of 6 represents 3 hours
                mean_wind_cube = test_wind_cube[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
                mean_theta_cube = test_theta_cube[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
                mean_q_cube = test_q_cube[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
                mean_ws_cube = test_ws_cube[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
                mean_sh_cube = test_sh_cube[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
                mean_lh_cube = test_lh_cube[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
                mean_w_cube = test_w_cube[ss*chunk:ss+chunk*ss].collapsed('time', iris.analysis.MEAN)
                print(mean_theta_cube)
                
                #%%
                print(mean_wind_cube.coords('surface_altitude'))
                surface = mean_wind_cube.coords('surface_altitude')[0].points
                x_coors = mean_wind_cube.coords('grid_latitude')[0].points
                print(surface)
                #%% set colorbar normalisations
                uni_cmap = 'bwr'             
                wsp_cmap = 'Oranges'
                w_cmap = 'bwr'
                theta_cmap = 'viridis'
                q_cmap = 'Blues'
                ws_cmap = 'Blues'
                sh_cmap = 'bwr'
                lh_cmap = 'bwr'
                ## atmostate
                wspmin = 0; wspmax = 25
                wmin = -3; wmax = -wmin
                Tmin = 270; Tmax = 295
                qmin = 0; qmax = 5
                q_levels = np.arange(0,50)*0.5
                wsp_norm = matplotlib.colors.Normalize(vmin = wspmin, vmax = wspmax)
                w_norm = matplotlib.colors.Normalize(vmin = wmin, vmax = wmax)
                theta_norm = matplotlib.colors.Normalize(vmin = Tmin, vmax = Tmax)
                q_norm = matplotlib.colors.Normalize(vmin = qmin, vmax = qmax)
                
                wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = wsp_cmap)
                w_map = matplotlib.cm.ScalarMappable(norm = w_norm, cmap = w_cmap)
                theta_map = matplotlib.cm.ScalarMappable(norm = theta_norm, cmap = theta_cmap)
                q_map = matplotlib.cm.ScalarMappable(norm = q_norm, cmap = q_cmap)
                
                ## fluxes
                wsmin = 0; wsmax = 2
                shmin = -250; shmax = -shmin
                lhmin = -250; lhmax = -lhmin
                
                ws_norm = matplotlib.colors.Normalize(vmin = wsmin, vmax = wsmax)
                sh_norm = matplotlib.colors.Normalize(vmin = shmin, vmax = shmax)
                lh_norm = matplotlib.colors.Normalize(vmin = lhmin, vmax = lhmax)
                
                ws_map = matplotlib.cm.ScalarMappable(norm = ws_norm, cmap = ws_cmap)
                sh_map = matplotlib.cm.ScalarMappable(norm = sh_norm, cmap = sh_cmap)
                lh_map = matplotlib.cm.ScalarMappable(norm = lh_norm, cmap = lh_cmap)
                
                #%%
                model_coors = mean_w_cube.coords('grid_latitude')[0]
                xmin = np.nanmin(model_coors.points)
                xmax = np.nanmax(model_coors.points)
                    
                #%% commence plotting
                        
                geo_axis = 'grid_latitude'
                common_keywargs = {'coords' : [geo_axis,'altitude'], 'levels' : 30}#,'norm':colors.CenteredNorm()}
                com_cbarkeys = {'pad' : 0.005}
                fig = plt.figure(figsize = (20,28))
                gs = gridspec.GridSpec(nrows = 7, ncols = 1)
                top_alt = 4000
                ## wsp
                
                ax0 = fig.add_subplot(gs[0,0])
                ax0.set_ylim(top = top_alt)
                ax0.set_title('horizontal windspeed')
                ax0.set_xticks([])
                ax0.set_ylabel('m')
                ax0.set_xlim(left = xmin, right = xmax)
                pwsp = iplt.contourf(mean_wind_cube, **common_keywargs, norm = wsp_norm, cmap = wsp_cmap)
                wsp_cbar = plt.colorbar(wsp_map, ax = ax0, **com_cbarkeys)
                wsp_cbar.ax.set(ylabel = r'$ms^{-1}$') 
                
                
                ax1 = fig.add_subplot(gs[1,0])
                ax1.set_ylim(top = top_alt)
                ax1.set_title('vertical windspeed')
                ax1.set_xticks([])
                ax1.set_ylabel('m')
                ax1.set_xlim(left = xmin, right = xmax)
                pw = iplt.contourf(mean_w_cube, **common_keywargs, cmap = w_cmap, norm = w_norm)
                w_cbar = plt.colorbar(w_map, ax = ax1, **com_cbarkeys)
                w_cbar.ax.set(ylabel = r'$ms^{-1}$')
                
                ax2 = fig.add_subplot(gs[2,0])
                ax2.set_ylim(top = top_alt)
                ax2.set_title('potential temperature')
                ax2.set_xticks([])
                ax2.set_ylabel('m')
                ax2.set_xlim(left = xmin, right = xmax)
                pw = iplt.contourf(mean_theta_cube, **common_keywargs, cmap = theta_cmap, norm = theta_norm)
                theta_cbar = plt.colorbar(theta_map, ax = ax2, **com_cbarkeys)
                theta_cbar.ax.set(ylabel = 'K')
                
                ax3 = fig.add_subplot(gs[3,0])
                ax3.set_ylim(top = top_alt)
                ax3.set_title('specific humidity')
                ax3.set_xticks([])
                ax3.set_ylabel('m')
                ax3.set_xlim(left = xmin, right = xmax)
                pw = iplt.contourf(mean_q_cube*1000, **common_keywargs, cmap = q_cmap, norm = q_norm)
                q_cbar = plt.colorbar(q_map, ax = ax3, **com_cbarkeys)
                q_cbar.ax.set(ylabel = '$gkg^{-1}$')
                
                ax4 = fig.add_subplot(gs[4,0])
                ax4.set_ylim(top = top_alt)
                ax4.set_title('windstress')
                ax4.set_xticks([])
                ax4.set_ylabel('m')
                ax4.set_xlim(left = xmin, right = xmax)
                pws = iplt.contourf(mean_ws_cube, **common_keywargs, cmap = ws_cmap, norm = ws_norm)
                ws_cbar = plt.colorbar(ws_map, ax = ax4, **com_cbarkeys)
                ws_cbar.ax.set(ylabel = '$Nm^{-2}$')
                
                ax5 = fig.add_subplot(gs[5,0])
                ax5.set_ylim(top = top_alt)
                ax5.set_title('sensible heat flux')
                ax5.set_xticks([])
                ax5.set_ylabel('m')
                ax5.set_xlim(left = xmin, right = xmax)
                psh = iplt.contourf(mean_sh_cube, **common_keywargs, cmap = sh_cmap, norm = sh_norm)
                sh_cbar = plt.colorbar(sh_map, ax = ax5, **com_cbarkeys)
                sh_cbar.ax.set(ylabel = '$Wm^{-2}$')
                
                ax6 = fig.add_subplot(gs[6,0])
                ax6.set_ylim(top = top_alt)
                ax6.set_title('latent heat flux')
                ax6.set_xticks([])
                ax6.set_ylabel('m')
                ax6.set_xlim(left = xmin, right = xmax)
                plh = iplt.contourf(mean_lh_cube, **common_keywargs, cmap = lh_cmap, norm = lh_norm)
                lh_cbar = plt.colorbar(lh_map, ax = ax6, **com_cbarkeys)
                lh_cbar.ax.set(ylabel = '$Wm^{-2}$')
                
                ax6.set_xticks(tick_coors)
                ax6.set_xticklabels(tick_dists)
                ax6.set_xlabel('Scalar distance, km')
                
                timepoint = mean_wind_cube.coord('time').bounds[0]
                print(dates[chunk*ss:ss+ss*chunk])
                boundstart = dates[chunk*ss:ss+ss*chunk][0].strftime('H%HM%M')
                print(boundstart)
                boundend = dates[chunk*ss:ss+ss*chunk][-1].strftime('H%HM%M')
                
                ## add in surface altitude outline
                for k in [ax0,ax1,ax2,ax3,ax4,ax5,ax6]:
                    k.plot(x_coors,surface, color = 'k', linewidth = 2)
                
                fig.suptitle(f'{config} - leg {leg} - 3 hour mean {boundstart}-{boundend}')
                
                plt.tight_layout()
                if True:
                    
                    plt.savefig(f'D:/Project/Figures/PNG/306/{experiment}/{suite}/vertical/averages/{res}/three_hour/leg{leg}/{config}_3houravr_vert_xsections_{boundstart}_306_leg{leg}.png')
                
                runstart = chunk*ss/48*24
                runend = runstart + 3
                if db_start > runstart and db_start < runend:
                    size = 280
                    ax0.scatter(rot_db_lat, db_alt, c = db_wsp, cmap = wsp_cmap, norm = wsp_norm, s = size,edgecolor = 'k')
                    ax1.scatter(rot_db_lat, db_alt, c = db_w, cmap = w_cmap, norm = w_norm, s = size,edgecolor = 'k')
                    ax2.scatter(rot_db_lat, db_alt, c = db_theta, cmap = theta_cmap, norm = theta_norm, s = size,edgecolor = 'k')
                    ax3.scatter(rot_db_lat, db_alt, c = db_q, cmap = q_cmap, norm = q_norm, s = size,edgecolor = 'k')
                    ax4.scatter(rot_db_lat, db_alt, c = db_ws, cmap = ws_cmap, norm = ws_norm, s = size,edgecolor = 'k')
                    ax5.scatter(rot_db_lat, db_alt, c = db_sh, cmap = sh_cmap, norm = sh_norm, s = size,edgecolor = 'k')
                    ax6.scatter(rot_db_lat, db_alt, c = db_lh, cmap = lh_cmap, norm = lh_norm, s = size,edgecolor = 'k')
                    if True:
                        
                        plt.savefig(f'D:/Project/Figures/PNG/306/{experiment}/{suite}/vertical/averages/{res}/three_hour/obs/{config}_3houravr_vert_xsections_{boundstart}_306_leg{leg}_with_obs.png')
                #plt.show()    
                plt.close()
                
        #%% start 