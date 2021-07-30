# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 13:04:25 2021

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

for res in ['0p5km','1p5km','4p4km']:
    #%% load model cubes
    suite = 'u-cc134'
    ## windstress (two cumponents, need to be combined)
    # wse = windstrss eastward, wsn = windstress northward
    wse_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_{res}_um_surface_downward_eastward_stress_24hrs_pg_306.nc', 'surface_downward_eastward_stress')
    wsn_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_{res}_um_surface_downward_northward_stress_24hrs_pg_306.nc','surface_downward_northward_stress')
    # regrid wsn cube and combine cubes for magnitude
    wsn_cube_reg = wsn_cube.regrid(wse_cube, iris.analysis.Linear())
    
    ws_cube = (wse_cube**2 + wsn_cube_reg**2)**0.5
    
    ## load heatfluxes
    # lh = latent heat, sh = sensible heat
    lh_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_{res}_um_surface_upward_latent_heat_flux_24hrs_pg_306.nc', 'surface_upward_latent_heat_flux')
    sh_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_{res}_um_surface_upward_sensible_heat_flux_24hrs_pg_306.nc', 'surface_upward_sensible_heat_flux')
    
    ## load surface altitude and land mask
    orog_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_{res}_um_surface_altitude_0_24hrs_pi_306.nc', 'surface_altitude')
    lbm_cube = iris.load_cube(f'D:/Project/Model_Data/{suite}/RA1M_{res}_um_land_binary_mask_24hrs_pi_306.nc','land_binary_mask')
    #%% subset cubes by 0p5km domain
    if res == '0p5km':
        gridlats = wse_cube.coord('grid_latitude').points
        gridlons = wse_cube.coord('grid_longitude').points
        W0p5km = gridlons[0]
        E0p5km = gridlons[-1]
        S0p5km = gridlats[200]
        N0p5km = gridlats[-1]
        print(S0p5km, N0p5km, W0p5km,E0p5km)
        interskwargs = {'grid_longitude':(W0p5km,E0p5km),'grid_latitude':(S0p5km,N0p5km)}
        
    
    ws_cube = ws_cube.intersection(**interskwargs)
    lh_cube = lh_cube.intersection(**interskwargs)
    sh_cube = sh_cube.intersection(**interskwargs)
    orog_cube = orog_cube.intersection(**interskwargs)
    lbm_cube = lbm_cube.intersection(**interskwargs)
    
    #%% extract data and coords from data structures
    wslat = ws_cube.coord('grid_latitude').points
    wslon = ws_cube.coord('grid_longitude').points
    lhlat = lh_cube.coord('grid_latitude').points
    lhlon = lh_cube.coord('grid_longitude').points
    shlat = sh_cube.coord('grid_latitude').points
    shlon = sh_cube.coord('grid_longitude').points
    oroglat = orog_cube.coord('grid_latitude').points
    oroglon = orog_cube.coord('grid_longitude').points
    wsdata = ws_cube.data
    lhdata = lh_cube.data
    shdata = sh_cube.data
    lbmdata = lbm_cube.data
    orogdata = orog_cube.data; orogdata[lbmdata == 0] = np.nan
    
    #%% get scalar distance ticklist for x axis labelling
    tick_coors, tick_dists = workshop.dist_tick_generator(ws_cube, 'grid_latitude', ticknumber = 6, decprecision = 2)
    #%% format date time for file name and figure title
    time = ws_cube.coord('time')
    dates = time.units.num2date(time.points)
    # print(dates)
    #%% add observations
    # using pandas
    obs_on = False
    if obs_on == True:
        i = 0
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
    
    
    for tidx in range(len(wsdata[:,0,0])):
        #%% set colorbar normalisations
                                   
        wsmin = 0; wsmax = np.nanpercentile(wsdata, 99.99)
        lhmin = -np.nanpercentile(np.positive(lhdata), 99.999); lhmax = -lhmin
        shmin = -np.nanpercentile(np.positive(shdata), 99.999); shmax = -shmin
        ormin = 0; ormax = np.nanpercentile(orogdata, 99.99)
        
        ws_norm = matplotlib.colors.Normalize(vmin = wsmin, vmax = wsmax)
        sh_norm = matplotlib.colors.Normalize(vmin = shmin, vmax = shmax)
        lh_norm = matplotlib.colors.Normalize(vmin = lhmin, vmax = lhmax)
        or_norm = matplotlib.colors.Normalize(vmin = ormin, vmax = ormax)
        
        ws_map = matplotlib.cm.ScalarMappable(norm = ws_norm, cmap = 'Blues')
        sh_map = matplotlib.cm.ScalarMappable(norm = sh_norm, cmap = 'bwr')
        lh_map = matplotlib.cm.ScalarMappable(norm = lh_norm, cmap = 'PuOr_r')
        or_map = matplotlib.cm.ScalarMappable(norm = or_norm, cmap = 'gist_earth')
        
        #%% commence plotting
        
        fig, ((ax0,ax1),(ax2,ax3)) = plt.subplots(2,2, figsize = (20,16))
        com_cbarkeys = {'pad' : 0.005}
        
        
        ## windstress
        ax0.pcolormesh(wslon,wslat,wsdata[tidx], cmap = 'Blues', vmin = wsmin, vmax = wsmax)
        ax0.contour(oroglon, oroglat, lbmdata, colors = 'black')
        ws_cbar = plt.colorbar(ws_map, ax = ax0, **com_cbarkeys)
        ws_cbar.ax.set(ylabel = r'Windstress $Pa$')
        ax0.set_xticks([])
        ax0.set_yticks([])
        
        ax1.pcolormesh(oroglon,oroglat,orogdata, cmap = 'gist_earth', vmin = ormin, vmax = ormax)
        or_cbar = plt.colorbar(or_map, ax = ax1, **com_cbarkeys)
        or_cbar.ax.set(ylabel = r'Surf. Alt. $m$')
        ax1.set_xticks([])
        ax1.set_yticks([])
        
        ax2.pcolormesh(shlon,shlat,shdata[tidx], cmap = 'bwr', vmin = shmin, vmax = shmax)
        sh_cbar = plt.colorbar(sh_map, ax = ax2, **com_cbarkeys)
        sh_cbar.ax.set(ylabel = r'SH $Wm^{-2}$')
        ax2.set_xticks([])
        ax2.set_yticks([])
        ax2.contour(oroglon, oroglat, lbmdata, colors = 'black')
        
        ax3.pcolormesh(lhlon,lhlat,lhdata[tidx], cmap = 'PuOr_r', vmin = lhmin, vmax = lhmax)
        lh_cbar = plt.colorbar(lh_map, ax = ax3, **com_cbarkeys)
        lh_cbar.ax.set(ylabel = r'LH $Wm^{-2}$')
        ax3.set_xticks([])
        ax3.set_yticks([])
        ax3.contour(oroglon, oroglat, lbmdata, colors = 'black')
        
        
        fig.tight_layout()
        
        savefigure = True
        if savefigure:
            timepoint = dates[tidx].strftime('H%HM%M')
            print(dates[tidx])
            config = 'RA1M'; exp = 'Control'
            plt.savefig(f'D:/Project/Figures/PDF/306/{suite}/Surface/fluxes/{res}/{config}_{res}_turbflux_surface_{exp}_{timepoint}_306.pdf')
            fig.suptitle(f'Surface turbulent flux variables, {res}, {dates[tidx]}')
            plt.savefig(f'D:/Project/Figures/PNG/306/{suite}/Surface/fluxes/{res}/{config}_{res}_turbflux_surface_{exp}_{timepoint}_306.png')
            
            
        plt.show()
    