# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 11:37:50 2021

@author: kse18nru
"""

#%% module imports
import sys
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
import model_foundry_toolkit as toolkit
import data_management as dmna
import time
import model_foundry as foundry
import matplotlib.gridspec as gridspec
#%% Global settings
plt.rcParams.update({'font.size': 20})


time_series_plot_atmostate = False
if time_series_plot_atmostate:
    testmeta = [(5,0,'A'),(5,-1, 'B'),(8,0,'A'), (8,-1,'B')]
    test_coordinates = []
    for test in testmeta:
        testnum = test[0]
        print('Initiate test', test)
        #%% load and process databse shortrun obs data
        # using pandas
        database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_sr_301.txt', delimiter = ' ')
        leg = testnum
        db_legno = database['legno']; condition = db_legno == testnum
        db_wsp = np.array(database['wsp'][condition])
        db_w = np.array(database['w'][condition])
        db_q = np.array(database['q_BUCK'][condition])
        db_theta = np.array(database['theta'][condition])
        db_lon = np.array(database['lon'][condition])
        db_lat = np.array(database['lat'][condition])
        db_alt = np.array(database['altgps'][condition])
        db_time = np.array(database['meantime'][condition])
        #%% rotate db coords
        compcube = iris.load_cube('../../Model_Data/u-bu807/nc/Control/RA1M_0p5km_umh1_flt301.nc', 'm01s15i002')[0,0]
        polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
        polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
        rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
        rot_db_lon = rot_db_lon + 360
        
        db_coor = rot_db_lat
        compcube = None
        #%% select coordniate for time series point
        tpoint_idx = test[1]
        tpoint_lat = rot_db_lat[tpoint_idx]
        tpoint_lon = rot_db_lon[tpoint_idx]
        tpoint_alt = db_alt[tpoint_idx]
        #print(db_lat[tpoint_idx], db_lon[tpoint_idx], db_alt[tpoint_idx])
        #print(tpoint_lat, tpoint_lon, tpoint_alt)
        
        #%% define resolutions and initialise figure
        resolutions = ['0p5km','1p5km', '4p4km']
        common_keywargs = {'s':300, 'color' : 'r'}
        fig, axes = plt.subplots(4,1, figsize = (20,10))
        #%% reprocess obs time into needed format
        for res in resolutions:
            #print(db_time[0]/60/60)
            time_tpoint = time.strftime('%H:%M',time.gmtime(db_time[0]))
            #print(time_tpoint)
            ## load data variables from vertical slice cubes ( simplifies picking the point for now )
            slice8_path = 'C:/Users/kse18nru/University/Project/Year_2/Model_Data/u-bu807/nc/Control/Vertical_slices/'
            #u_slice8 = iris.load_cube(f'{slice8_path}m01s15i002/RA1M_0p5km_m01s15i002_flt301_leg8.nc', 'm01s15i002')
            #v_slice8 = iris.load_cube(f'{slice8_path}m01s15i003/RA1M_0p5km_m01s15i003_flt301_leg8.nc', 'm01s15i003')
            wsp_slice8_0p5 = dmna.load_model_wsp(res, 'm01s15i002', 'm01s15i003', 'Vertical_slices', 8,301, 'Control', level_range = 35)
            w_slice8_0p5 =  dmna.load_model_w(res,'Vertical_slices', 8,301, 'Control', level_range = 35)
            
            u_3dcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_{res}_um_m01s15i002_24hrs_h_301.nc', 'm01s15i002')
            v_3dcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_{res}_um_m01s15i003_24hrs_h_301.nc', 'm01s15i003')
            wsp_3dcube = (u_3dcube**2 + v_3dcube**2)**0.5
            w_3dcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_{res}_um_upward_air_velocity_24hrs_h_301.nc', 'upward_air_velocity')
            th_3dcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_{res}_um_air_potential_temperature_24hrs_i_301.nc', 'air_potential_temperature')
            q_3dcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_{res}_um_specific_humidity_24hrs_i_301.nc', 'specific_humidity')
            
            
            
            
            #%% interpolate time series
            
            wsp_series, wsptcoor = pert.interpolate_timeseries(wsp_3dcube, tpoint_lon, tpoint_lat, tpoint_alt)
            w_series, wtcoor = pert.interpolate_timeseries(w_3dcube, tpoint_lon, tpoint_lat, tpoint_alt)
            th_series, thtcoor = pert.interpolate_timeseries(th_3dcube, tpoint_lon, tpoint_lat, tpoint_alt)
            q_series, qtcoor = pert.interpolate_timeseries(q_3dcube, tpoint_lon, tpoint_lat, tpoint_alt)
            
            
            #%% reformat time coordinate
            cube_time = w_3dcube.coord('time')
            dates = cube_time.units.num2date(cube_time.points)
            timepoints = []
            for t in range(len(dates)):
                timepoint = dates[t].strftime('%H:%M')
                timepoints.append(timepoint)
            times = np.array(timepoints)
            #print(times)
            
            #%%
            #print((cube_time.points - cube_time.points[0])*60*60)
            model_time_ssm = (cube_time.points - cube_time.points[0])*60*60 + (30*60) # convert to seconds since last midnight (1st model point is 00:30)
            #print(type(time_tpoint))
            time_tpoint = db_time[0]
            #%% plot time series
            
            ## windspeed plot
            axes[0].plot(model_time_ssm, wsp_series, label = f'{res}')
            
            ## vert windspeed plot
            axes[1].plot(model_time_ssm, w_series, label = f'{res}')
            
            ## theta plot
            axes[2].plot(model_time_ssm, th_series, label = f'{res}')
            
            
            
            ## specific humidity plot
            axes[3].plot(model_time_ssm, q_series*1000, label = f'{res}')
            print(f'{res} cycle complete')
        #%%    
        axes[0].set_ylabel(r'$wsp$, ms$^{-1}$')
        axes[0].set_xticks(model_time_ssm[1::4])
        axes[0].set_xticklabels([])
        axes[0].grid(True) 
        axes[0].legend()
        axes[0].scatter([time_tpoint], [db_wsp[tpoint_idx]], **common_keywargs)
        axes[1].set_ylim(top = 1.2, bottom = -1.2)
        axes[1].set_ylabel(r'$w$, ms$^{-1}$')
        axes[1].set_xticks(model_time_ssm[1::4])
        axes[1].set_xticklabels([])
        axes[1].grid(True)
        #axes[1].legend()
        axes[1].scatter([time_tpoint], [db_w[tpoint_idx]], **common_keywargs)
        axes[2].set_ylabel(r'$\theta$, K')
        axes[2].set_xticks(model_time_ssm[1::4])
        axes[2].set_xticklabels([])
        axes[2].grid(True)
        #axes[2].legend()
        axes[2].scatter([time_tpoint], [db_theta[tpoint_idx]], **common_keywargs)    
        axes[3].set_ylabel(r'$q$, gkg$^{-1}$')
        axes[3].set_xticks(model_time_ssm[1::4])
        axes[3].set_xticklabels(times[1::4])
        axes[3].grid(True)
        axes[3].set_xlabel('Model time, UTC')
        #axes[3].legend()
        axes[3].scatter([time_tpoint], [db_q[tpoint_idx]], **common_keywargs)
        axes[0].set_title(f'Control atmospheric state variable time series at test point {testnum}{test[2]} ({db_lat[tpoint_idx]:.3f},{db_lon[tpoint_idx]:.3f},{tpoint_alt:.0f}m)')
        plt.subplots_adjust(wspace=0, hspace=0.1)
        #plt.tight_layout()
            
        plt.savefig(f'../../Figures/PNG/301/u-bu807/SlowBubble/atmostate__control_tseries_testpoint_{testnum}{test[2]}.png')
        plt.show()
        #%%
        test_coordinates.append([f'{testnum}{test[2]}',tpoint_lon, tpoint_lat, db_lon[test[1]], db_lat[test[1]], tpoint_alt])
    df = pd.DataFrame(test_coordinates, columns = ['label','grid_lon','grid_lat','lon','lat','alt_m'])
    df.to_csv('../../Figures/PNG/301/u-bu807/SlowBubble/tpoint_coors_301.csv')

    
#%%
time_series_contour_atmostate = True
if time_series_contour_atmostate:
    testmeta = [(5,0,'A'),(5,-1, 'B'),(8,0,'A'), (8,-1,'B'), (99,99,99), (99,99,99), (99,99,99), (99,99,99), (99,99,99)]
    test_coordinates = []
    i = 0
    for i in range(9):
        df = pd.read_csv('../../Figures/PNG/301/u-bu807/SlowBubble/tpoint_coors_301.csv')
        label = df['label'][i]
        
        test = testmeta[i]
            
        testnum = test[0]
        if testnum < 9:
            print('Initiate test', test)
       
            database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_sr_301.txt', delimiter = ' ')
            leg = testnum
            db_legno = database['legno']; condition = db_legno == testnum
            db_wsp = np.array(database['wsp'][condition])
            db_w = np.array(database['w'][condition])
            db_q = np.array(database['q_BUCK'][condition])
            db_theta = np.array(database['theta'][condition])
            db_lon = np.array(database['lon'][condition])
            db_lat = np.array(database['lat'][condition])
            db_alt = np.array(database['altgps'][condition])
            db_time = np.array(database['meantime'][condition])
        
            #%% rotate db coords

            compcube = iris.load_cube('../../Model_Data/u-bu807/nc/Control/RA1M_0p5km_umh1_flt301.nc', 'm01s15i002')[0,0]
            polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
            polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
            rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
            rot_db_lon = rot_db_lon + 360
            
            db_coor = rot_db_lat
            compcube = None
            #%% select coordniate for time series point
            tpoint_idx = test[1]
            tpoint_lat = rot_db_lat[tpoint_idx]
            tpoint_lon = rot_db_lon[tpoint_idx]
            tpoint_alt = db_alt[tpoint_idx]
            #print(db_lat[tpoint_idx], db_lon[tpoint_idx], db_alt[tpoint_idx])
            #print(tpoint_lat, tpoint_lon, tpoint_alt)
        
            
        #%% define resolutions and initialise figure
        resolutions = ['0p5km','1p5km', '4p4km']
        common_keywargs = {'s':300, 'color' : 'r'}
        
        #%% reprocess obs time into needed format
        for res in resolutions:
            #print(db_time[0]/60/60)
            if testnum < 9:
                time_tpoint = time.strftime('%H:%M',time.gmtime(db_time[0]))
            
            
    
            u_3dcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_{res}_um_m01s15i002_24hrs_h_301.nc', 'm01s15i002')
            v_3dcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_{res}_um_m01s15i003_24hrs_h_301.nc', 'm01s15i003')
            wsp_3dcube = (u_3dcube**2 + v_3dcube**2)**0.5
            w_3dcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_{res}_um_upward_air_velocity_24hrs_h_301.nc', 'upward_air_velocity')
            th_3dcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_{res}_um_air_potential_temperature_24hrs_i_301.nc', 'air_potential_temperature')
            q_3dcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_{res}_um_specific_humidity_24hrs_i_301.nc', 'specific_humidity')
            p_3dcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/RA1M_{res}_um_air_pressure_24hrs_i_301.nc', 'air_pressure')
            
            
            
            #%% interpolate time series
            if testnum < 9:
                wsp_series = pert.interpolate_timeseries_coloumns(wsp_3dcube, tpoint_lon, tpoint_lat)
                w_series = pert.interpolate_timeseries_coloumns(w_3dcube, tpoint_lon, tpoint_lat)
                th_series = pert.interpolate_timeseries_coloumns(th_3dcube, tpoint_lon, tpoint_lat)
                q_series = pert.interpolate_timeseries_coloumns(q_3dcube, tpoint_lon, tpoint_lat)
                p_series = pert.interpolate_timeseries_coloumns(p_3dcube, tpoint_lon, tpoint_lat)
            else:
                
                glon = df['grid_lon'][i]
                glat = df['grid_lat'][i]
                alt = df['alt_m'][i]
                wsp_series = pert.interpolate_timeseries_coloumns(wsp_3dcube, glon,glat)
                w_series = pert.interpolate_timeseries_coloumns(w_3dcube, glon, glat)
                th_series = pert.interpolate_timeseries_coloumns(th_3dcube, glon, glat)
                q_series = pert.interpolate_timeseries_coloumns(q_3dcube, glon, glat)
                p_series = pert.interpolate_timeseries_coloumns(p_3dcube, glon, glat)
                
            #%%
            wspmin = 0; wspmax = 22
            wmin = -2; wmax = -wmin
            Tmin = 268; Tmax = 282
            qmin = 0; qmax = 2.6
            
            wsp_norm = matplotlib.colors.Normalize(vmin = wspmin, vmax = wspmax)
            w_norm = matplotlib.colors.Normalize(vmin = wmin, vmax = wmax)
            theta_norm = matplotlib.colors.Normalize(vmin = Tmin, vmax = Tmax)
            q_norm = matplotlib.colors.Normalize(vmin = qmin, vmax = qmax)
            
            wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = 'Oranges')
            w_map = matplotlib.cm.ScalarMappable(norm = w_norm, cmap = 'seismic')
            theta_map = matplotlib.cm.ScalarMappable(norm = theta_norm, cmap = 'plasma')
            q_map = matplotlib.cm.ScalarMappable(norm = q_norm, cmap = 'Blues')
            #%% reformat time coordinate
            cube_time = w_series.coord('time')
            dates = cube_time.units.num2date(cube_time.points)
            timepoints = []
            for t in range(len(dates)):
                timepoint = dates[t].strftime('%H:%M')
                timepoints.append(timepoint)
            times = np.array(timepoints)
            print(times)
            print(cube_time.points)
            print(w_series.coord('forecast_period').points)
            #%%
            #print((cube_time.points - cube_time.points[0])*60*60)
            model_time_ssm = (cube_time.points - cube_time.points[0])*60*60 + (30*60) # convert to seconds since last midnight (1st model point is 00:30)
            #print(type(time_tpoint))
            time_tpoint = db_time[0]
            #%% plot time series
            print(wsp_series[:,:,0])
            cbar_keywargs = {'orientation':'vertical', 'pad' : 0.01}
            uni_kwargs = {'coords' : ['forecast_period', 'altitude'],'levels' : 30}
            ylim_kw = {'top' : 4000, 'bottom' : 0}
            top_level = 35
            yticks = [1000,2000,3000,4000]
            fig = plt.figure(figsize = (20,20))
            
            gs = gridspec.GridSpec(nrows = 4, ncols = 1)
            ax0 = fig.add_subplot(gs[0,0])
            iplt.contourf(wsp_series[:,:top_level,0], cmap = 'Oranges', vmin = wspmin, vmax = wspmax, **uni_kwargs)
            cbar0 = plt.colorbar(wsp_map, ax = ax0, **cbar_keywargs)
            ax0.set_title(f'{res} {label}')
            cbar0.ax.set(ylabel = r'wsp $ms^{-1}$')
            ax0.set_xticks(w_series.coord('forecast_period').points[1::4])
            ax0.set_xticklabels([])
            ax0.set_ylim(**ylim_kw)
            ax0.set_yticks(yticks)
            plt.grid(axis = 'x')
            
            ax1 = fig.add_subplot(gs[1,0])
            iplt.contourf(w_series[:,:top_level,0], cmap = 'seismic', vmin = wmin, vmax = wmax, **uni_kwargs)
            cbar1 = plt.colorbar(w_map, ax = ax1, **cbar_keywargs)
            cbar1.ax.set(ylabel = r'w $ms^{-1}$')
            ax1.set_xticks(w_series.coord('forecast_period').points[1::4])
            ax1.set_xticklabels([])
            ax1.set_ylim(**ylim_kw)
            ax1.set_yticks(yticks)
            plt.grid(axis = 'x')
            
            ax2 = fig.add_subplot(gs[2,0])
            iplt.contourf(th_series[:,:top_level,0], cmap = 'viridis', vmin = Tmin, vmax = Tmax, **uni_kwargs)
            cbar2 = plt.colorbar(theta_map, ax = ax2, **cbar_keywargs)
            cbar2.ax.set(ylabel = r'Theta $K$')
            ax2.set_xticks(w_series.coord('forecast_period').points[1::4])
            ax2.set_xticklabels([])
            ax2.set_ylim(**ylim_kw)
            ax2.set_yticks(yticks)
            plt.grid(axis = 'x')
            
            ax3 = fig.add_subplot(gs[3,0])
            iplt.contourf(q_series[:,:top_level,0]*1000, cmap = 'Blues', vmin = qmin, vmax = qmax, **uni_kwargs)
            cbar3 = plt.colorbar(q_map, ax = ax3, **cbar_keywargs)
            cbar3.ax.set(ylabel = r'Sp. Hum. $gkg^{-1}$')
            ax3.set_yticks(yticks)
            
            ## make tick list
            ax3.set_xticks(w_series.coord('forecast_period').points[1::4])
            ax3.set_xticklabels(times[1::4])
            ax3.set_ylim(**ylim_kw)
            plt.grid(axis = 'x')
            plt.subplots_adjust(wspace=0.06, hspace=0.05)
            png_path = f'../../Figures/PNG/301/u-bu807/Vertical/timeseries/Control/testpoints/{label}/{res}_control_atmostate_timexsections_tpoint_{label}.png'
            plt.savefig(png_path)
        
            plt.show()
        
            
  #%%
test_point_map_plot = False
if test_point_map_plot:
    #%% load test point coordinates from file
    df = pd.read_csv('../../Figures/PNG/301/u-bu807/SlowBubble/tpoint_coors_301.csv')
    
    print(df)
    tpoint_labels = np.array(df['label'])
    tpoint_lon = np.array(df['grid_lon'])
    tpoint_lat = np.array(df['grid_lat'])
    tpoint_alt = np.array(df['alt_m'])
    concat_labels = []
    print(f'{tpoint_alt[0]:.0f}')
    cl_n = 0
    for cl in tpoint_labels:
        newstring = f'{cl}, {tpoint_alt[cl_n]:.0f}m'
        cl_n +=1
        concat_labels.append(newstring)
    #%% rotate db coords
    compcube = iris.load_cube('../../Model_Data/u-bu807/nc/Control/RA1M_0p5km_umh1_flt301.nc', 'm01s15i002')[0,0]
    polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    
    compcube = None
    
    #%% load orography
    orogcube = iris.load_cube('../../Model_Data/u-bu807/nc/Control/surface_altitude/RA1M_0p5km_surface_altitude_fulldomain_flt301.nc', 'surface_altitude')
    print(orogcube)
    gridlon = orogcube.coords('grid_longitude')[0].points
    gridlat = orogcube.coords('grid_latitude')[0].points
    lons, lats = foundry.modf.unrotate_coords(gridlon, gridlat, polelon, polelat)
    land_cube = iris.load_cube('../../Model_Data/u-bu807/nc/Control/surface_altitude/RA1M_0p5km_surface_altitude_fulldomain_flt301.nc', 'land_binary_mask')
    lm_data = land_cube.data
    or_data = orogcube.data
    #%% mask out sea points from surface altitude
    #print(lm_data)
    or_data[lm_data == 0] = np.nan
    #print(or_data)
    sea_points = (lm_data -1)*(-1)
    #%% get scalar distance ticklist for x axis labelling
    tick_lon_coors, tick_lon_dists = workshop.dist_tick_generator(orogcube, 'grid_longitude', ticknumber = 6, decprecision = 2, flat = True)
    tick_lat_coors, tick_lat_dists = workshop.dist_tick_generator(orogcube, 'grid_latitude', ticknumber = 6, decprecision = 2, flat = True)
    print(tick_lon_coors, tick_lon_dists)
    #%%
    fig,ax = plt.subplots(1,1, figsize = (15,20))
    ax.pcolormesh(gridlon, gridlat, sea_points, cmap = 'Spectral')
    ax.pcolormesh(gridlon, gridlat, or_data, cmap = 'summer')
    #ax.contourf(orogcube.coords('grid_longitude')[0].points, orogcube.coords('grid_latitude')[0].points,or_data)
    ax.scatter(tpoint_lon,tpoint_lat, s = 1000, c = 'k', marker = '+', linewidths = 8)
    q4 = ax.contour(gridlon, gridlat, lons, colors = 'k')
    q5 = ax.contour(gridlon, gridlat, lats, colors = 'k')
    label_keys = {'fontsize' : 10, 'inline' : 1}
        
    q6 = ax.clabel(q4, **label_keys)
    q7 = ax.clabel(q5, **label_keys)
    ax.set_xticks(tick_lon_coors)
    ax.set_xticklabels(tick_lon_dists)
    ax.set_yticks(tick_lat_coors)
    ax.set_yticklabels(tick_lat_dists)
    ax.set_xlabel('km'); ax.set_ylabel('km')
    ax.set_title('Locations of test points for in situ time series, IGP case study 301')
    for i, txt in enumerate(concat_labels):
        ax.annotate(tpoint_labels[i], (tpoint_lon[i]-0.05, tpoint_lat[i]-0.1))
    plt.savefig('../../Figures/PNG/301/u-bu807/SlowBubble/tpoint_locations_301.png')
    plt.show()
#%%
time_series_plot_turbprop = False
if time_series_plot_turbprop:
    placeholder = 0
    
#%%
resolutions = ['0p5km', '1p5km', '4p4km']
horz_pressure_plot = False
if horz_pressure_plot:
    for res in resolutions:
        
        subset_name = 'fulldomain'
        th_cube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/interp4d/RA1M_{res}_air_potential_temperature_interp_{subset_name}_flt301_alt500to1200_tidx23to36.nc', 'air_potential_temperature')
        p_cube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/interp4d/RA1M_{res}_air_pressure_interp_{subset_name}_flt301_alt500to1200_tidx23to36.nc', 'air_pressure')
        u_cube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/interp4d/RA1M_{res}_m01s15i002_interp_{subset_name}_flt301_alt500to1200_tidx23to36.nc', 'm01s15i002')
        v_cube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/interp4d/RA1M_{res}_m01s15i003_interp_{subset_name}_flt301_alt500to1200_tidx23to36.nc', 'm01s15i003')
        wsp_cube = (u_cube**2 + v_cube**2)**0.5
        w_cube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/interp4d/RA1M_{res}_upward_air_velocity_interp_{subset_name}_flt301_alt500to1200_tidx23to36.nc', 'upward_air_velocity')
        q_cube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/interp4d/RA1M_{res}_specific_humidity_interp_{subset_name}_flt301_alt500to1200_tidx23to36.nc', 'specific_humidity')
        print(np.shape(p_cube.data))
        land_cube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/surface_altitude/RA1M_{res}_surface_altitude_{subset_name}_flt301.nc', 'land_binary_mask')
        orog_cube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/surface_altitude/RA1M_{res}_surface_altitude_{subset_name}_flt301.nc', 'surface_altitude')
        #%% extract data and coordinate arrays
        th_data = th_cube.data
        p_data = p_cube.data
        wsp_data = wsp_cube.data
        w_data = w_cube.data
        u_data = u_cube.data
        v_data = v_cube.data
        wsp_lat = u_cube.coords('grid_latitude')[0].points
        wsp_lon = u_cube.coords('grid_longitude')[0].points
        lat_coor = w_cube.coords('grid_latitude')[0].points
        lon_coor = w_cube.coords('grid_longitude')[0].points
        
        q_data = q_cube.data
        or_data = orog_cube.data
        lm_data = land_cube.data
        
        #%% time handling
        midnight_hour = orog_cube.coords('time')[0][0]
        print(midnight_hour)
        dimtimes = q_cube.coords('time')[0].points -midnight_hour.points
        
        
        
        altitude = q_cube.coords('altitude')[0].points
        print(altitude)
        #%% mask out sea points from surface altitude
        #print(lm_data)
        or_data[lm_data == 0] = np.nan
        #print(or_data)
        sea_points = (lm_data -1)*(-1)
        print(type(sea_points))
        #sea_points[lm_data == 1] = np.nan
        #%% create scalar distanc axis labels and posotions
        
        tick_lon_coors, tick_lon_dists = workshop.dist_tick_generator(p_cube, 'grid_longitude', ticknumber = 5, decprecision = 0, flat = True)
        tick_lat_coors, tick_lat_dists = workshop.dist_tick_generator(p_cube, 'grid_latitude', ticknumber = 8, decprecision = 0, flat = True)
        print(tick_lon_coors, tick_lon_dists)
        #%%
        wspmin = 0; wspmax = 22
        wmin = -7; wmax = -wmin
        Tmin = 268; Tmax = 282
        pmin = 888; pmax = 891
        qmin = 0; qmax = 2.6
        wsp_norm = matplotlib.colors.Normalize(vmin = wspmin, vmax = wspmax)
        w_norm = matplotlib.colors.Normalize(vmin = wmin, vmax = wmax)
        theta_norm = matplotlib.colors.Normalize(vmin = Tmin, vmax = Tmax)
        p_norm = matplotlib.colors.Normalize(vmin = pmin, vmax = pmax)
        q_norm = matplotlib.colors.Normalize(vmin = qmin, vmax = qmax)
        
        wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = 'Oranges')
        w_map = matplotlib.cm.ScalarMappable(norm = w_norm, cmap = 'seismic')
        theta_map = matplotlib.cm.ScalarMappable(norm = theta_norm, cmap = 'plasma')
        p_map = matplotlib.cm.ScalarMappable(norm = p_norm, cmap = 'viridis')
        q_map = matplotlib.cm.ScalarMappable(norm = q_norm, cmap = 'Blues')
        #%%
        for t in range(14):
            for level in range(3):
                #%% plot a horizontal x-section
                
                if res == '0p5km':
                    figsize = (18,28)
                if res == '1p5km':
                    figsize = (20,20)
                if res == '4p4km':
                    figsize = (20,20)
                cbar_keywargs = {'orientation':'vertical'}
                #fig, ax = plt.subplots(1,1, figsize = (20,20))
                fig = plt.figure(figsize = figsize)
                gs = gridspec.GridSpec(nrows = 3, ncols = 2)
                ax0 = fig.add_subplot(gs[0,0])
                q0 = ax0.pcolormesh(lon_coor, lat_coor,p_data[t][level]/100, cmap = 'viridis')#, vmin = pmin, vmax = pmax)
                cbar0 = plt.colorbar(q0, ax = ax0, **cbar_keywargs)
                cbar0.ax.set(ylabel = r'Pressure $hPa$')
                ax0.set_xticks([])
                ax0.set_yticks(tick_lat_coors)
                ax0.set_yticklabels(tick_lat_dists)
                
                ax1 = fig.add_subplot(gs[0,1])
                q1 = ax1.pcolormesh(lon_coor, lat_coor,th_data[t][level], cmap = 'plasma', vmin = Tmin, vmax = Tmax) 
                cbar1 = plt.colorbar(theta_map, ax = ax1, **cbar_keywargs)
                cbar1.ax.set(ylabel = r'Theta $K$')
                ax1.set_xticks([])
                ax1.set_yticks([])
                if subset_name == 'fulldomain':
                    if res == '0p5km': steps = 17; 
                    if res == '1p5km' : steps = 15; 
                    if res == '4p4km' : steps = 15;
                if subset_name == 'southbay':
                    if res == '0p5km': steps = 12; 
                    if res == '1p5km' : steps = 4; 
                    if res == '4p4km' : steps = 2;
                
                ax2 = fig.add_subplot(gs[1,0])
                q2 = ax2.pcolormesh(wsp_lon, wsp_lat,wsp_data[t][level], cmap = 'Oranges', vmin = wspmin, vmax = wspmax)
                cbar2 = plt.colorbar(wsp_map, ax = ax2, **cbar_keywargs)
                cbar2.ax.set(ylabel = r'Hoz. Windspeed $ms^{-1}$')
                ax2.quiver(wsp_lon[::steps],wsp_lat[::steps],u_data[t,level,::steps,::steps],v_data[t,level,::steps,::steps])
                ax2.set_xticks([])
                
                ax2.set_yticks(tick_lat_coors)
                ax2.set_yticklabels(tick_lat_dists)
                
                ax3 = fig.add_subplot(gs[1,1])
                q3 = ax3.pcolormesh(lon_coor, lat_coor,w_data[t][level], cmap = 'seismic', vmin = wmin, vmax = wmax)
                cbar3 = plt.colorbar(w_map, ax = ax3, **cbar_keywargs)
                cbar3.ax.set(ylabel = r'Vert. Windsp. $ms^{-1}$')
                ax3.set_yticks([])
                ax3.set_xticks([])
                
                
                ax4 = fig.add_subplot(gs[2,0])
                q4 = ax4.pcolormesh(lon_coor, lat_coor, q_data[t][level]*1000, cmap = 'Blues')
                cbar4 = plt.colorbar(q_map, ax = ax4, **cbar_keywargs)
                cbar4.ax.set(ylabel = r'Spec. Humidity $gkg^{-1}$')
                ax4.set_xticks(tick_lon_coors)
                ax4.set_xticklabels(tick_lon_dists)
                ax4.set_yticks(tick_lat_coors)
                ax4.set_yticklabels(tick_lat_dists)
                
                ax5 = fig.add_subplot(gs[2,1])
                ax5.pcolormesh(lon_coor, lat_coor, sea_points, cmap = 'Spectral')
                q5 = ax5.pcolormesh(lon_coor, lat_coor, or_data, cmap = 'summer')
                cbar5 = plt.colorbar(q5, ax = ax5, cmap = 'summer', **cbar_keywargs)
                cbar5.ax.set(ylabel = r'Surface Alt. $m$')
                ax5.set_yticks([])
                ax5.set_xticks(tick_lon_coors)
                ax5.set_xticklabels(tick_lon_dists)
                
                plt.subplots_adjust(wspace=0.06, hspace=0.05)
                
                #%% title
                dectime = dimtimes[t]
                hour = int(dectime)
                minute = int((dectime*60)%60)
                if minute == 0:
                    minute = '00'
                title = f'Res {res} Alt {altitude[level]}m Time {hour}:{minute} UTC'
                fig.suptitle(title, y = 0.91)
                #%%
                pngpath = f'../../Figures/PNG/301/u-bu807/M_levels/Control/Interpolated/{subset_name}/{res}/{altitude[level]}m/RA1M_{res}_atmostate_interpolated_{subset_name}_alt{altitude[level]}m_H{hour}M{minute}.png'
                plt.savefig(pngpath)
                plt.show()