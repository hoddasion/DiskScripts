# -*- coding: utf-8 -*-
"""
Created on Wed May 12 16:06:19 2021

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
resolutions = ['0p5km','1p5km', '4p4km']
varname_u = 'x_wind'
varname_v = 'y_wind'
varname_wsp = 'windspeed'
varname_w = 'upward_air_velocity'
varname_theta = 'air_potential_temperature'
varname_q = 'specific_humidity'

lat_lon_ord = ('grid_latitude', 'grid_longitude')
lon_lat_ord = ('grid_longitude', 'grid_latitude')

config = 'RA1M'
exp = 'Control'
suite = 'u-cc134'
alt_idx = 35
long = 'long'
#%% load cubes from files
mass_plots = True
if mass_plots:
    for res in resolutions:
        for leg in [7]:#[1,4,7,13,15]:  
            
            
                
            #%%
            try:
               
                u_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/x_wind/{config}_vertslice_{res}_x_wind_flt306_leg{leg}{long}.nc', 'x_wind')[:,:alt_idx]
                v_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/y_wind/{config}_vertslice_{res}_y_wind_flt306_leg{leg}{long}.nc', 'y_wind')[:,:alt_idx]
                theta_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/air_potential_temperature/{config}_vertslice_{res}_air_potential_temperature_flt306_leg{leg}{long}.nc', 'air_potential_temperature')[:,:alt_idx]
                q_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/specific_humidity/{config}_vertslice_{res}_specific_humidity_flt306_leg{leg}{long}.nc', 'specific_humidity')[:,:alt_idx]
                w_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/upward_air_velocity/{config}_vertslice_{res}_upward_air_velocity_flt306_leg{leg}{long}.nc','upward_air_velocity')[:,:alt_idx]
           
                #%% process u_cube and v_cube into wsp_cube by completely circumventing and ignoring stupid metadata issue
                wsp_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/x_wind/{config}_vertslice_{res}_x_wind_flt306_leg{leg}{long}.nc', 'x_wind')[:,:alt_idx]
                x_data = u_cube.data; y_data = v_cube.data
                wsp_data = (x_data**2 + y_data**2)**0.5
                wsp_cube.data = wsp_data
                #wsp_cube.standard_name = 'windspeed'
                #%% get scalar distance ticklist for x axis labelling
                tick_coors, tick_dists = workshop.dist_tick_generator(w_cube, 'grid_latitude', ticknumber = 6, decprecision = 2)
                for tidx in range(48):
                    #%% format date time for file name and figure title
                    time = w_cube.coord('time')
                    dates = time.units.num2date(time.points)
                    print(dates[tidx])
                    
                    #%% set colorbar normalisations
                            
                    wspmin = 0; wspmax = np.nanpercentile(wsp_data, 99.99)
                    wmin = -7; wmax = -wmin
                    Tmin = np.nanpercentile(theta_cube.data, 0.01); Tmax = np.nanpercentile(theta_cube.data, 99.99)
                    qmin = 0; qmax = 6
                    q_levels = np.arange(0,50)*0.5
                    wsp_norm = matplotlib.colors.Normalize(vmin = wspmin, vmax = wspmax)
                    w_norm = matplotlib.colors.Normalize(vmin = wmin, vmax = wmax)
                    theta_norm = matplotlib.colors.Normalize(vmin = Tmin, vmax = Tmax)
                    q_norm = matplotlib.colors.Normalize(vmin = qmin, vmax = qmax)
                    
                    wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = 'Oranges')
                    w_map = matplotlib.cm.ScalarMappable(norm = w_norm, cmap = 'seismic')
                    theta_map = matplotlib.cm.ScalarMappable(norm = theta_norm, cmap = 'plasma')
                    q_map = matplotlib.cm.ScalarMappable(norm = q_norm, cmap = 'Blues')
                    #%% dynamically adjust horizontal bounds
                    model_coors = w_cube.coords('grid_latitude')[0]
                    xmin = np.nanmin(model_coors.points)
                    xmax = np.nanmax(model_coors.points)
                    
                    
                    #%% commence plotting
                    #tidx = 24
                    geo_axis = 'grid_latitude'
                    common_keywargs = {'coords' : [geo_axis,'altitude'], 'levels' : 40}
                    com_cbarkeys = {'pad' : 0.005}
                    fig = plt.figure(figsize = (24,20))
                    gs = gridspec.GridSpec(nrows = 4, ncols = 1)
                    top_alt = 4000
                    ## wsp
                    
                    ax0 = fig.add_subplot(gs[0,0])
                    ax0.set_ylim(top = top_alt)
                    ax0.set_title('horizontal windspeed')
                    ax0.set_xticks([])
                    ax0.set_ylabel('m')
                    ax0.set_xlim(left = xmin, right = xmax)
                    pwsp = iplt.contourf(wsp_cube[tidx], **common_keywargs, cmap = 'Oranges', vmin = wspmin, vmax  = wspmax)
                    #plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                    wsp_cbar = plt.colorbar(wsp_map, ax = ax0, **com_cbarkeys)
                    wsp_cbar.ax.set(ylabel = r'$ms^{-1}$')
                    ## wsp obs
                    #ax0.scatter(db_coor, db_alt, c = db_wsp, cmap = 'Oranges', vmin = wspmin, vmax = wspmax,
                    #            s = 300, edgecolor = 'k')
                    ## w
                    ax1 = fig.add_subplot(gs[1,0])
                    pw = iplt.contourf(w_cube[tidx], **common_keywargs, cmap = 'seismic', vmin = wmin, vmax = wmax)
                    ax1.set_ylim(top = top_alt)
                    ax1.set_title('upward windspeed')
                    ax1.set_xticks([])
                    ax1.set_ylabel('m')
                    ax1.set_xlim(left = xmin, right = xmax)
                    #plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                    w_cbar = plt.colorbar(w_map, ax = ax1, **com_cbarkeys)
                    w_cbar.ax.set(ylabel = r'$ms^{-1}$')
                    ## w obs
                    #ax1.scatter(db_coor, db_alt, c = db_w, cmap = 'seismic', vmin = wmin, vmax = wmax,
                    #            s = 300, edgecolor = 'k')
                    ## theta
                    
                    ax2 = fig.add_subplot(gs[2,0])
                    ptheta = iplt.contourf(theta_cube[tidx], **common_keywargs,cmap = 'plasma', vmin = Tmin, vmax = Tmax)
                    ax2.set_ylim(top = top_alt)
                    ax2.set_title('potential temperature')
                    ax2.set_xticks([])
                    ax2.set_ylabel('m')
                    ax2.set_xlim(left = xmin, right = xmax)
                    #plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                    theta_cbar = plt.colorbar(theta_map, ax = ax2, **com_cbarkeys)
                    theta_cbar.ax.set(ylabel = 'K')
                    ## theta obs
                    #ax2.scatter(db_coor, db_alt, c = db_theta, cmap = 'plasma', vmin = Tmin, vmax = Tmax,
                    #            s = 300, edgecolor = 'k')
                    ## q
                    ax3 = fig.add_subplot(gs[3,0])
                    pq = iplt.contourf(q_cube[tidx]*1000, **common_keywargs, cmap = 'Blues', vmin = qmin, vmax = qmax)
                    ax3.set_ylim(top = top_alt)
                    ax3.set_title('specific humidity')
                    ax3.set_xticks(tick_coors)
                    ax3.set_xticklabels(tick_dists)
                    ax3.set_xlabel('Scalar distance, km')
                    ax3.set_ylabel('m')
                    ax3.set_xlim(left = xmin, right = xmax)
                    #plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                    q_cbar = plt.colorbar(q_map, ax = ax3, **com_cbarkeys)
                    q_cbar.ax.set(ylabel = '$gkg^{-1}$')
                    ## q obs
                    #ax3.scatter(db_coor, db_alt, c = db_q, cmap = 'Blues', vmin = qmin, vmax = qmax,
                    #            s = 300, edgecolor = 'k')
                    fig.suptitle(f'Atmospheric state variables, {res}, leg {leg}{long}, {dates[tidx]}')
                    fig.tight_layout()
                    save_figure = True
                    if save_figure:
                        timepoint = dates[tidx].strftime('H%HM%M')
                        plt.savefig(f'D:/Project/Figures/PNG/306/{suite}/Vertical/atmostate/{res}/leg{leg}{long}/{config}_{res}_atmostate_vert_xsection_{exp}_{timepoint}_leg{leg}{long}.png')
                    
                    #plt.show()
                    #%% add observations
                    # using pandas
                    database = pd.read_csv('D:/Project/Obvs_Data/Databases/IGP_flights_database_obs_sr_306.txt', delimiter = ' ')
                    
                    db_legno = database['legno']
                    legcond = db_legno == leg
                    db_time = np.array(database['meantime'][legcond])
                    mid_point = np.mean(db_time)
                    model_time_ssm = (time.points-time.points[0])*60*60 # get seconds since midnight from model output timesteps
                    model_time_hsm = (time.points-time.points[0]) # hours since midnight
                    
                    
                        
                    if mid_point/(60*60) <= model_time_hsm[tidx] + 0.25 and mid_point/(60*60) > model_time_hsm[tidx] - 0.25:
                        print('Yay found the closest time! It is', mid_point/(60*60), 'and', model_time_hsm[tidx])
                        db_wsp = np.array(database['wsp'][legcond]) 
                        db_w = np.array(database['w'][legcond])
                        db_q = np.array(database['q_BUCK'][legcond])
                        db_theta = np.array(database['theta'][legcond])
                        
                        db_lon = np.array(database['lon'][legcond])
                        db_lat = np.array(database['lat'][legcond])
                        db_alt = np.array(database['altgps'][legcond])
                        
                        #%% rotate db coords
                        compcube = iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_0p5km_um_upward_air_velocity_24hrs_ph_306.nc', 'upward_air_velocity')[0,0]
                        polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
                        polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
                        rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
                        rot_db_lon = rot_db_lon + 360
                        db_coor = rot_db_lat
                        
                        #%% set colorbar normalisations
                            
                        wspmin = 0; wspmax = np.nanpercentile(wsp_data[tidx,:25], 99.999)
                        wmin = -7; wmax = -wmin
                        Tmin = np.nanpercentile(theta_cube.data[tidx,:25], 0.001); Tmax = np.nanpercentile(theta_cube.data[tidx,:25], 99.999)
                        qmin = 0; qmax = 6
                        q_levels = np.arange(0,50)*0.5
                        wsp_norm = matplotlib.colors.Normalize(vmin = wspmin, vmax = wspmax)
                        w_norm = matplotlib.colors.Normalize(vmin = wmin, vmax = wmax)
                        theta_norm = matplotlib.colors.Normalize(vmin = Tmin, vmax = Tmax)
                        q_norm = matplotlib.colors.Normalize(vmin = qmin, vmax = qmax)
                        
                        wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = 'Oranges')
                        w_map = matplotlib.cm.ScalarMappable(norm = w_norm, cmap = 'seismic')
                        theta_map = matplotlib.cm.ScalarMappable(norm = theta_norm, cmap = 'plasma')
                        q_map = matplotlib.cm.ScalarMappable(norm = q_norm, cmap = 'Blues')
                        
                        #%% commence plotting
                        #tidx = 24
                        geo_axis = 'grid_latitude'
                        common_keywargs = {'coords' : [geo_axis,'altitude'], 'levels' : 35}
                        com_cbarkeys = {'pad' : 0.005}
                        fig = plt.figure(figsize = (24,20))
                        gs = gridspec.GridSpec(nrows = 4, ncols = 1)
                        top_alt = 2500
                        ## wsp
                        
                        ax0 = fig.add_subplot(gs[0,0])
                        ax0.set_ylim(top = top_alt)
                        ax0.set_title('horizontal windspeed')
                        ax0.set_xticks([])
                        ax0.set_ylabel('m')
                        ax0.set_xlim(left = xmin, right = xmax)
                        pwsp = iplt.contourf(wsp_cube[tidx], **common_keywargs, cmap = 'Oranges', vmin = wspmin, vmax  = wspmax)
                        #plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                        wsp_cbar = plt.colorbar(wsp_map, ax = ax0, **com_cbarkeys)
                        wsp_cbar.ax.set(ylabel = r'$ms^{-1}$')
                        ## wsp obs
                        ax0.scatter(db_coor, db_alt, c = db_wsp, cmap = 'Oranges', vmin = wspmin, vmax = wspmax,
                                    s = 300, edgecolor = 'k')
                        ## w
                        ax1 = fig.add_subplot(gs[1,0])
                        pw = iplt.contourf(w_cube[tidx], **common_keywargs, cmap = 'seismic', vmin = wmin, vmax = wmax)
                        ax1.set_ylim(top = top_alt)
                        ax1.set_title('upward windspeed')
                        ax1.set_xticks([])
                        ax1.set_ylabel('m')
                        ax1.set_xlim(left = xmin, right = xmax)
                        #plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                        w_cbar = plt.colorbar(w_map, ax = ax1, **com_cbarkeys)
                        w_cbar.ax.set(ylabel = r'$ms^{-1}$')
                        ## w obs
                        ax1.scatter(db_coor, db_alt, c = db_w, cmap = 'seismic', vmin = wmin, vmax = wmax,
                                    s = 300, edgecolor = 'k')
                        ## theta
                        
                        ax2 = fig.add_subplot(gs[2,0])
                        ptheta = iplt.contourf(theta_cube[tidx], **common_keywargs,cmap = 'plasma', vmin = Tmin, vmax = Tmax)
                        ax2.set_ylim(top = top_alt)
                        ax2.set_title('potential temperature')
                        ax2.set_xticks([])
                        ax2.set_ylabel('m')
                        ax2.set_xlim(left = xmin, right = xmax)
                        #plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                        theta_cbar = plt.colorbar(theta_map, ax = ax2, **com_cbarkeys)
                        theta_cbar.ax.set(ylabel = 'K')
                        ## theta obs
                        ax2.scatter(db_coor, db_alt, c = db_theta, cmap = 'plasma', vmin = Tmin, vmax = Tmax,
                                    s = 300, edgecolor = 'k')
                        ## q
                        ax3 = fig.add_subplot(gs[3,0])
                        pq = iplt.contourf(q_cube[tidx]*1000, **common_keywargs, cmap = 'Blues', vmin = qmin, vmax = qmax)
                        ax3.set_ylim(top = top_alt)
                        ax3.set_title('specific humidity')
                        ax3.set_xticks(tick_coors)
                        ax3.set_xticklabels(tick_dists)
                        ax3.set_xlabel('Scalar distance, km')
                        ax3.set_ylabel('m')
                        ax3.set_xlim(left = xmin, right = xmax)
                        #plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                        q_cbar = plt.colorbar(q_map, ax = ax3, **com_cbarkeys)
                        q_cbar.ax.set(ylabel = '$gkg^{-1}$')
                        ## q obs
                        ax3.scatter(db_coor, db_alt, c = db_q, cmap = 'Blues', vmin = qmin, vmax = qmax,
                                    s = 300, edgecolor = 'k')
                        fig.suptitle(f'Atmospheric state variables, {res}, leg {leg}, {dates[tidx]}')
                        fig.tight_layout()
                        save_figure = True
                        if save_figure:
                            timepoint = dates[tidx].strftime('H%HM%M')
                            plt.savefig(f'D:/Project/Figures/PNG/306/{suite}/Vertical/atmostate/obs/{config}_{res}_atmostate_vert_xsection_obs_{exp}_{timepoint}_leg{leg}{long}.png')
                            plt.savefig(f'D:/Project/Figures/PDF/306/{suite}/Vertical/atmostate/obs/{config}_{res}_atmostate_vert_xsection_obs_{exp}_{timepoint}_leg{leg}{long}.pdf')
                        
                        #plt.show()
                        plt.close()
            except:
                print(f'Exception occured; pass on leg {leg}...\n')
                pass  
                    
#%% plot aggregated obs xsections

for res in resolutions:
    for i in range(1):
        #try:
            exit()
            cube_legs = [7]#[1,4,7,13,15]
            
            u_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/x_wind/{config}_vertslice_{res}_x_wind_flt306_leg{cube_legs[i]}{long}.nc', 'x_wind')[:,:alt_idx]
            v_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/y_wind/{config}_vertslice_{res}_y_wind_flt306_leg{cube_legs[i]}{long}.nc', 'y_wind')[:,:alt_idx]
            theta_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/air_potential_temperature/{config}_vertslice_{res}_air_potential_temperature_flt306_leg{cube_legs[i]}{long}.nc', 'air_potential_temperature')[:,:alt_idx]
            q_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/specific_humidity/{config}_vertslice_{res}_specific_humidity_flt306_leg{cube_legs[i]}{long}.nc', 'specific_humidity')[:,:alt_idx]
            w_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/upward_air_velocity/{config}_vertslice_{res}_upward_air_velocity_flt306_leg{cube_legs[i]}{long}.nc','upward_air_velocity')[:,:alt_idx]
            #%% process u_cube and v_cube into wsp_cube by completely circumventing and ignoring stupid metadata issue
            wsp_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/x_wind/{config}_vertslice_{res}_x_wind_flt306_leg{cube_legs[i]}{long}.nc', 'x_wind')[:,:alt_idx]
            x_data = u_cube.data; y_data = v_cube.data
            wsp_data = (x_data**2 + y_data**2)**0.5
            wsp_cube.data = wsp_data
            #%% get scalar distance ticklist for x axis labelling
            tick_coors, tick_dists = workshop.dist_tick_generator(w_cube, 'grid_latitude', ticknumber = 6, decprecision = 2)
            #%% format date time for file name and figure title
            time = w_cube.coord('time')
            dates = time.units.num2date(time.points)
           # print(dates)
            #%% add observations
            # using pandas
            database = pd.read_csv('D:/Project/Obvs_Data/Databases/IGP_flights_database_obs_60s_306.txt', delimiter = ' ')
            
            db_legno = database['legno']
            dbl = db_legno
            #print(dbl)
            dbl_list = [np.array(dbl == 1) + np.array(dbl == 2),
                         np.array(dbl == 4) + np.array(dbl == 5) + np.array(dbl == 6),
                         np.array(dbl == 7) + np.array(dbl == 8) + np.array(dbl == 9),
                         dbl == 13,
                         dbl == 15]
            legcond = dbl_list[i]
            db_time = np.array(database['meantime'][legcond])
            mid_point = np.nanmean(db_time)
            model_time_ssm = (time.points-time.points[0]+0.5)*60*60 # get seconds since midnight from model output timesteps
            model_time_hsm = (time.points-time.points[0]+0.5) # hours since midnight
            print(mid_point/(60*60))
            
            for tidx in range(48): 
                # go through each time step and check which is valid in next line
                diff = mid_point/(60*60) - model_time_hsm[tidx]
                if diff < 0.5 and diff > -0.5:
                    #specfically checking for which leg the midpoint falls within 30 minutes of the model output.
                    print('Yay found the closest time! It is', mid_point/(60*60), 'and', model_time_hsm[tidx])
                    ## now load remaining db data and coord arrays
                    db_wsp = np.array(database['wsp'][legcond]) 
                    db_w = np.array(database['w'][legcond])
                    db_q = np.array(database['q_BUCK'][legcond])
                    db_theta = np.array(database['theta'][legcond])
                    
                    db_lon = np.array(database['lon'][legcond])
                    db_lat = np.array(database['lat'][legcond])
                    db_alt = np.array(database['altgps'][legcond])
                            
                    #%% rotate db coords
                    compcube = iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_0p5km_um_upward_air_velocity_24hrs_ph_306.nc', 'upward_air_velocity')[0,0]
                    polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
                    polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
                    rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
                    rot_db_lon = rot_db_lon + 360
                    db_coor = rot_db_lat
                    
                    #%% now proceed with the plotting!
                    #%% set colorbar normalisations
                                
                    wspmin = 0; wspmax = np.nanpercentile(wsp_data[tidx,:25], 99.999)
                    wmin = -7; wmax = -wmin
                    Tmin = np.nanpercentile(theta_cube.data[tidx,:25], 0.001); Tmax = np.nanpercentile(theta_cube.data[tidx,:25], 99.999)
                    qmin = 0; qmax = 6
                    q_levels = np.arange(0,50)*0.5
                    wsp_norm = matplotlib.colors.Normalize(vmin = wspmin, vmax = wspmax)
                    w_norm = matplotlib.colors.Normalize(vmin = wmin, vmax = wmax)
                    theta_norm = matplotlib.colors.Normalize(vmin = Tmin, vmax = Tmax)
                    q_norm = matplotlib.colors.Normalize(vmin = qmin, vmax = qmax)
                    
                    wsp_map = matplotlib.cm.ScalarMappable(norm = wsp_norm, cmap = 'Oranges')
                    w_map = matplotlib.cm.ScalarMappable(norm = w_norm, cmap = 'seismic')
                    theta_map = matplotlib.cm.ScalarMappable(norm = theta_norm, cmap = 'plasma')
                    q_map = matplotlib.cm.ScalarMappable(norm = q_norm, cmap = 'Blues')
                    
                    #%% dynamically adjust horizontal bounds
                    model_coors = w_cube.coords('grid_latitude')[0]
                    xmin = np.nanmin(model_coors.points)
                    xmax = np.nanmax(model_coors.points)
                    
                    #%% commence plotting
                    #tidx = 24
                    geo_axis = 'grid_latitude'
                    common_keywargs = {'coords' : [geo_axis,'altitude'], 'levels' : 35}
                    com_cbarkeys = {'pad' : 0.005}
                    fig = plt.figure(figsize = (24,20))
                    gs = gridspec.GridSpec(nrows = 4, ncols = 1)
                    top_alt = 2500
                    ## wsp
                    
                    ax0 = fig.add_subplot(gs[0,0])
                    ax0.set_ylim(top = top_alt)
                    ax0.set_title('horizontal windspeed')
                    ax0.set_xticks([])
                    ax0.set_ylabel('m')
                    ax0.set_xlim(left = xmin, right = xmax)
                    pwsp = iplt.contourf(wsp_cube[tidx], **common_keywargs, cmap = 'Oranges', vmin = wspmin, vmax  = wspmax)
                    #plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                    wsp_cbar = plt.colorbar(wsp_map, ax = ax0, **com_cbarkeys)
                    wsp_cbar.ax.set(ylabel = r'$ms^{-1}$')
                    ## wsp obs
                    ax0.scatter(db_coor, db_alt, c = db_wsp, cmap = 'Oranges', vmin = wspmin, vmax = wspmax,
                                s = 300, edgecolor = 'k')
                    ## w
                    ax1 = fig.add_subplot(gs[1,0])
                    pw = iplt.contourf(w_cube[tidx], **common_keywargs, cmap = 'seismic', vmin = wmin, vmax = wmax)
                    ax1.set_ylim(top = top_alt)
                    ax1.set_title('upward windspeed')
                    ax1.set_xticks([])
                    ax1.set_ylabel('m')
                    ax1.set_xlim(left = xmin, right = xmax)
                    #plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                    w_cbar = plt.colorbar(w_map, ax = ax1, **com_cbarkeys)
                    w_cbar.ax.set(ylabel = r'$ms^{-1}$')
                    ## w obs
                    ax1.scatter(db_coor, db_alt, c = db_w, cmap = 'seismic', vmin = wmin, vmax = wmax,
                                s = 300, edgecolor = 'k')
                    ## theta
                    
                    ax2 = fig.add_subplot(gs[2,0])
                    ptheta = iplt.contourf(theta_cube[tidx], **common_keywargs,cmap = 'plasma', vmin = Tmin, vmax = Tmax)
                    ax2.set_ylim(top = top_alt)
                    ax2.set_title('potential temperature')
                    ax2.set_xticks([])
                    ax2.set_ylabel('m')
                    ax2.set_xlim(left = xmin, right = xmax)
                    #plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                    theta_cbar = plt.colorbar(theta_map, ax = ax2, **com_cbarkeys)
                    theta_cbar.ax.set(ylabel = 'K')
                    ## theta obs
                    ax2.scatter(db_coor, db_alt, c = db_theta, cmap = 'plasma', vmin = Tmin, vmax = Tmax,
                                s = 300, edgecolor = 'k')
                    ## q
                    ax3 = fig.add_subplot(gs[3,0])
                    pq = iplt.contourf(q_cube[tidx]*1000, **common_keywargs, cmap = 'Blues', vmin = qmin, vmax = qmax)
                    ax3.set_ylim(top = top_alt)
                    ax3.set_title('specific humidity')
                    ax3.set_xticks(tick_coors)
                    ax3.set_xticklabels(tick_dists)
                    ax3.set_xlabel('Scalar distance, km')
                    ax3.set_ylabel('m')
                    ax3.set_xlim(left = xmin, right = xmax)
                    #plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                    q_cbar = plt.colorbar(q_map, ax = ax3, **com_cbarkeys)
                    q_cbar.ax.set(ylabel = '$gkg^{-1}$')
                    ## q obs
                    ax3.scatter(db_coor, db_alt, c = db_q, cmap = 'Blues', vmin = qmin, vmax = qmax,
                                s = 300, edgecolor = 'k')
                    legnames = ['1_2', '4_5_6','7_8_9','13','15']
                    fig.suptitle(f'Atmospheric state variables, {res}, leg {legnames[i]}{long}, {dates[tidx]}')
                    fig.tight_layout()
                    save_figure = True
                    if save_figure:
                        timepoint = dates[tidx].strftime('H%HM%M')
                        print(dates[tidx])
                        plt.savefig(f'D:/Project/Figures/PNG/306/{suite}/Vertical/atmostate/obs/{legnames[i]}/{config}_{res}_atmostate_vert_xsection_obs_{exp}_{timepoint}_leg{legnames[i]}{long}.png')
                        plt.savefig(f'D:/Project/Figures/PDF/306/{suite}/Vertical/atmostate/obs/{legnames[i]}/{config}_{res}_atmostate_vert_xsection_obs_{exp}_{timepoint}_leg{legnames[i]}{long}.pdf')
                    
                    plt.show()
        #except:
         #   print(f'Exception occured for {res} {i}')