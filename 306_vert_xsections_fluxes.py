# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 17:45:06 2021

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
resolutions = ['0p5km']#,'1p5km', '4p4km']
var_we = 'atmosphere_downward_eastward_stress'
var_wn = 'atmosphere_downward_northward_stress'
var_sh = 'upward_heat_flux_in_air'
var_q = 'specific_humidity'
var_theta = 'air_potential_temperature'
var_p = 'air_pressure'
var_vap = 'upward_water_vapor_flux_in_air'

lat_lon_ord = ('grid_latitude', 'grid_longitude')
lon_lat_ord = ('grid_longitude', 'grid_latitude')

config = 'LONGTAIL'
exp = 'LONGTAIL'
suite = 'u-cf117'
alt_idx = 35

for res in resolutions:
    for i in range(5):
        #try:
            cube_legs = [1,4,7,13,15]
            
            
            theta_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/air_potential_temperature/{config}_vertslice_{res}_air_potential_temperature_flt306_leg{cube_legs[i]}.nc', 'air_potential_temperature')[:,:alt_idx]
            q_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/specific_humidity/{config}_vertslice_{res}_specific_humidity_flt306_leg{cube_legs[i]}.nc', 'specific_humidity')[:,:alt_idx]
           
            #%% load model data
            
            ## windstress (two cumponents, need to be combined)
            # wse = windstrss eastward, wsn = windstress northward
            wse_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/{var_we}/{config}_vertslice_{res}_{var_we}_flt306_leg{cube_legs[i]}.nc', var_we)[:,:alt_idx]
            wsn_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/{var_wn}/{config}_vertslice_{res}_{var_wn}_flt306_leg{cube_legs[i]}.nc',var_wn)[:,:alt_idx]
            # regrid wsn cube and combine cubes for magnitude
            #%%
            
            ws_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/{var_we}/{config}_vertslice_{res}_{var_we}_flt306_leg{cube_legs[i]}.nc', var_we)[:,:alt_idx]
            ws_cube.data = (wse_cube.data**2 + wsn_cube.data**2)**0.5
            #%%
            ## load heatfluxes
            # lh = latent heat, sh = sensible heat
            #lh_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_{res}_um_surface_upward_latent_heat_flux_24hrs_pg_306.nc', 'atmosphere_upward_latent_heat_flux')[:,:alt_idx]
            sh_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/{var_sh}/{config}_vertslice_{res}_{var_sh}_flt306_leg{cube_legs[i]}.nc', 'upward_heat_flux_in_air')[:,:alt_idx]
            # calculate lh
            #%%
            
            lh_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/Vertical/upward_latent_heat_flux_in_air/{config}_vertslice_{res}_upward_latent_heat_flux_in_air_flt306_leg{cube_legs[i]}.nc', 'upward_latent_heat_flux_in_air')[:,:alt_idx]
            
            
            net_cube = lh_cube + sh_cube
            #%% get scalar distance ticklist for x axis labelling
            tick_coors, tick_dists = workshop.dist_tick_generator(wse_cube, 'grid_latitude', ticknumber = 6, decprecision = 2)
            #%% format date time for file name and figure title
            time = wse_cube.coord('time')
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
                    db_ws = np.array(database['windstress'][legcond]) 
                    db_lh = np.array(database['lh'][legcond])
                    db_sh = np.array(database['sh'][legcond])
                    #db_theta = np.array(database['theta'][legcond])
                    
                    db_lon = np.array(database['lon'][legcond])
                    db_lat = np.array(database['lat'][legcond])
                    db_alt = np.array(database['altgps'][legcond])
                            
                    #%% rotate db coords
                    compcube = iris.load_cube(f'D:/Project/Model_Data/{suite}/{config}_0p5km_um_upward_air_velocity_24hrs_ph_306.nc', 'upward_air_velocity')[0,0]
                    polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
                    polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
                    rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
                    rot_db_lon = rot_db_lon + 360
                    db_coor = rot_db_lat
                    
                    #%% now proceed with the plotting!
                    #%% set colorbar normalisations
                                
                    #%% set colorbar normalisations
                                   
                    wsmin = 0; wsmax = np.nanpercentile(ws_cube.data[tidx], 99.99)
                    lhmin = -np.nanpercentile(np.positive(lh_cube.data[tidx]), 99.99); lhmax = -lhmin
                    shmin = -np.nanpercentile(np.positive(sh_cube.data[tidx]), 99.99); shmax = -shmin
                    netmin = -np.nanpercentile(np.positive(sh_cube.data[tidx] + lh_cube.data[tidx]), 99.99); netmax = -netmin
                    #lhmin = -200; lhmax = -lhmin
                    #shmin = -200; shmax = -shmin
                    #netmin = -200; netmax = -netmin
                    
                    ws_norm = matplotlib.colors.Normalize(vmin = wsmin, vmax = wsmax)
                    sh_norm = matplotlib.colors.Normalize(vmin = shmin, vmax = shmax)
                    lh_norm = matplotlib.colors.Normalize(vmin = lhmin, vmax = lhmax)
                    net_norm = matplotlib.colors.Normalize(vmin = netmin, vmax = netmax)
                    
                    ws_map = matplotlib.cm.ScalarMappable(norm = ws_norm, cmap = 'Blues')
                    sh_map = matplotlib.cm.ScalarMappable(norm = sh_norm, cmap = 'bwr')
                    lh_map = matplotlib.cm.ScalarMappable(norm = lh_norm, cmap = 'PuOr_r')
                    net_map = matplotlib.cm.ScalarMappable(norm = net_norm, cmap = 'PiYG')
                    
                    #%% dynamically adjust horizontal bounds
                    model_coors = wse_cube.coords('grid_latitude')[0]
                    xmin = np.nanmin(model_coors.points)
                    xmax = np.nanmax(model_coors.points)
                    
                    #%% commence plotting
                    #tidx = 24
                    geo_axis = 'grid_latitude'
                    common_keywargs = {'coords' : [geo_axis,'altitude']}
                    com_cbarkeys = {'pad' : 0.005}
                    fig = plt.figure(figsize = (24,20))
                    gs = gridspec.GridSpec(nrows = 4, ncols = 1)
                    top_alt = 2500
                    heat_levels = np.arange(-400,400,10)
                    print(heat_levels)
                    ## wsp
                    
                    
                    #%%
                    ax0 = fig.add_subplot(gs[0,0])
                    ax0.set_ylim(top = top_alt)
                    ax0.set_title('windstress')
                    ax0.set_xticks([])
                    ax0.set_ylabel('m')
                    ax0.set_xlim(left = xmin, right = xmax)
                    pwsp = iplt.contourf(ws_cube[tidx], **common_keywargs, cmap = 'Blues', vmin = wsmin, vmax  = wsmax, levels = 35)
                    #plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                    wsp_cbar = plt.colorbar(ws_map, ax = ax0, **com_cbarkeys)
                    wsp_cbar.ax.set(ylabel = r'$Pa$')
                    ## wsp obs
                    ax0.scatter(db_coor, db_alt, c = db_ws, cmap = 'Blues', vmin = wsmin, vmax = wsmax,
                                s = 300, edgecolor = 'k')
                    ## w
                    ax1 = fig.add_subplot(gs[1,0])
                    pw = iplt.contourf(sh_cube[tidx], **common_keywargs, cmap = 'bwr', vmin = shmin, vmax = shmax, levels = heat_levels)
                    ax1.set_ylim(top = top_alt)
                    ax1.set_title('upward sensible heat flux')
                    ax1.set_xticks([])
                    ax1.set_ylabel('m')
                    ax1.set_xlim(left = xmin, right = xmax)
                    #plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                    w_cbar = plt.colorbar(sh_map, ax = ax1, **com_cbarkeys)
                    w_cbar.ax.set(ylabel = r'$Wm^{-2}$')
                    ## w obs
                    ax1.scatter(db_coor, db_alt, c = db_sh, cmap = 'bwr', vmin = shmin, vmax = shmax,
                                s = 300, edgecolor = 'k')
                    ## theta
                    
                    ax2 = fig.add_subplot(gs[2,0])
                    ptheta = iplt.contourf(lh_cube[tidx], **common_keywargs,cmap = 'PuOr_r', vmin = lhmin, vmax = lhmax, levels = heat_levels)
                    ax2.set_ylim(top = top_alt)
                    ax2.set_title('upward latent heat flux')
                    ax2.set_xticks([])
                    ax2.set_ylabel('m')
                    ax2.set_xlim(left = xmin, right = xmax)
                    #plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                    theta_cbar = plt.colorbar(lh_map, ax = ax2, **com_cbarkeys)
                    theta_cbar.ax.set(ylabel = r'$Wm^{-2}$')
                    ## theta obs
                    ax2.scatter(db_coor, db_alt, c = db_lh, cmap = 'PuOr_r', vmin = lhmin, vmax = lhmax,
                                s = 300, edgecolor = 'k')
                    ## q
                    ax3 = fig.add_subplot(gs[3,0])
                    pq = iplt.contourf(net_cube[tidx], **common_keywargs, cmap = 'PiYG', vmin = netmin, vmax = netmax, levels = heat_levels)
                    ax3.set_ylim(top = top_alt)
                    ax3.set_title('upward latent + sensible heat flux')
                    ax3.set_xticks(tick_coors)
                    ax3.set_xticklabels(tick_dists)
                    ax3.set_xlabel('Scalar distance, km')
                    ax3.set_ylabel('m')
                    ax3.set_xlim(left = xmin, right = xmax)
                    #plt.plot(orogcoor, orogdata, color = 'k', linestyle = '--')
                    q_cbar = plt.colorbar(net_map, ax = ax3, **com_cbarkeys)
                    q_cbar.ax.set(ylabel = r'$Wm^{-2}$')
                    ## q obs
                    ax3.scatter(db_coor, db_alt, c = db_lh + db_sh, cmap = 'PiYG', vmin = netmin, vmax = netmax,
                                s = 300, edgecolor = 'k')
                    legnames = ['1_2', '4_5_6','7_8_9','13','15']
                    
                    fig.tight_layout()
                    save_figure = True
                    if save_figure:
                        timepoint = dates[tidx].strftime('H%HM%M')
                        print(dates[tidx])
                        #plt.savefig(f'D:/Project/Figures/PDF/306/{suite}/Vertical/fluxes/obs/{legnames[i]}/{config}_{res}_turbflux_vert_xsection_obs_{exp}_{timepoint}_leg{legnames[i]}.pdf')
                        fig.suptitle(f'{config} {res} Turbulent flux variables, leg {legnames[i]}, {dates[tidx]}')
                        plt.savefig(f'D:/Project/Figures/PNG/306/{exp}/{suite}/vertical/fluxes/obs/{legnames[i]}/{config}_{res}_turbflux_vert_xsection_obs_{exp}_{timepoint}_leg{legnames[i]}_dynamicscale.png')
                        
                    
                    plt.show()
        #except:
         #   print(f'Exception occured for {res} {i}')