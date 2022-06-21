# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 18:45:51 2022

@author: kse18nru
"""

#%%
import iris
import iris.plot as iplt
from iris.analysis.cartography import rotate_pole
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from iris.analysis import trajectory
from iris.coords import DimCoord, AuxCoord
from iris.analysis.cartography import get_xy_grids
import timeit

#%% globals
flight = 306

#%%
if flight == 306:
    suite = 'u-cc134'
    config = 'RA1M'
    
    stream = 'pi'
    res = '0p5km'
    cube_path = f'D:/Project/Model_Data/{suite}/'
    theta_filename = f'{config}_{res}_um_air_potential_temperature_24hrs_pi_{flight}.nc'
    w_filename = f'{config}_{res}_um_upward_air_velocity_24hrs_ph_{flight}.nc'
    
    theta_cube = iris.load_cube(f'{cube_path}{theta_filename}','air_potential_temperature')
    w_cube = iris.load_cube(f'{cube_path}{w_filename}','upward_air_velocity')
    lbm_cube = iris.load_cube(f'{cube_path}{config}_{res}_um_land_binary_mask_24hrs_pi_{flight}.nc', 'land_binary_mask')
   

    #%% extract current cube grid coords
    grid_lon = theta_cube.coords('grid_longitude')[0].points
    grid_lat =theta_cube.coords('grid_latitude')[0].points
    print(grid_lon, grid_lat)
    
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
    wmin = -2; wmax = -wmin
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
    
    #%% plot something for reference
    tidx = 30
    hidx = 5
    
    fig, (ax0,ax1)  = plt.subplots(1,2, figsize = (12,10))
    
    ax0.contourf(grid_lon, grid_lat, theta_cube[tidx, hidx].data, levels = 20)
    ax0.contour(grid_lon, grid_lat,lbm_cube.data, colors = 'k')
    ax0.grid()
    
    ax1.contourf(grid_lon, grid_lat, w_cube[tidx, hidx].data, levels = 20, cmap = 'seismic',norm = w_norm)
    ax1.contour(grid_lon, grid_lat,lbm_cube.data, colors = 'k')
    ax1.grid()
    
    #%%
    #sea_w = np.where(lbm_cube.data == 0, w_cube.data[:,:])
    #print(sea_w)
    for tidx in [30]:
        for hidx in [5]:
            slice_w = w_cube.data[tidx, hidx]
            slice_theta = theta_cube.data[tidx, hidx]
            #%% mask land area out
            slice_w[lbm_cube.data == 1] = np.nan
            slice_theta[lbm_cube.data == 1] = np.nan
    
    
            analysis_subdomain = True
            if analysis_subdomain:
                print('Commencing analysis of subdomain mean fields')
                #%% subset cubes
                # define vertices
                
                gridSW = (359.5, 0)
                gridNE = (361, 1.5)
                # xgrid, ygrid = get_xy_grids(w_cube)
                # print(np.shape(xgrid), np.shape(ygrid))
                # xconds = (xgrid < gridNE[0]) * (xgrid > gridSW[0])
                # yconds = (ygrid < gridNE[1]) * (ygrid > gridSW[1])
                # xyconds = xconds * yconds
                xconds = (grid_lon < gridNE[0]) * (grid_lon > gridSW[0])
                yconds = (grid_lat < gridNE[1]) * (grid_lat > gridSW[1])
                print(np.shape(xconds), np.shape(yconds))
                print(np.shape(w_cube.data[0,0]))
                #local_gridlon = 
                local_w = slice_w[yconds,:]
                local_w = local_w[:, xconds]
                local_theta = slice_theta[yconds,:]
                local_theta = local_theta[:,xconds]
                local_lon = grid_lon[xconds]
                local_lat = grid_lat[yconds]
                #%%
                tidx = 30
                hidx = 5
                pcol_kwargs = {'shading':'auto'}
                fig, (ax0,ax1)  = plt.subplots(1,2, figsize = (12,10))
                
                ax0.pcolormesh(local_lon, local_lat, local_theta, **pcol_kwargs)
                ax0.contour(grid_lon, grid_lat,lbm_cube.data, colors = 'k')
                ax0.grid()
                
                ax1.pcolormesh(local_lon, local_lat, local_w, cmap = 'seismic',norm = w_norm,**pcol_kwargs)
                ax1.contour(grid_lon, grid_lat,lbm_cube.data, colors = 'k')
                ax1.grid()
                
                #%% take mean of fields
                mean_theta = np.nanmean(local_theta)
                std_theta = np.nanstd(local_theta)
                print(r'mean theta =',mean_theta, '+/-', std_theta)
                mean_w = np.nanmean(local_w)
                std_w = np.nanstd(local_w)
                print(r'mean w =',mean_w, '+/-', std_w)
                pert_theta = local_theta - mean_theta # perturbation field of theta
                pert_w = local_w - mean_w # perturbation field of w
                mean_pert_theta = np.nanmean(pert_theta)
                mean_pert_w = np.nanmean(pert_w)
                print(r'mean pert theta =',mean_pert_theta)
                print(r'mean pert w =',mean_pert_w)
                co_field = pert_theta*pert_w
                print(r'mean co field =', np.nanmean(co_field))
                #%% normalise fields for plotting
                pT_min =  -np.nanmax(np.abs(pert_theta)); pT_max = -pT_min
                pw_min = -np.nanpercentile(np.abs(pert_w),99.5); pw_max = -pw_min
                co_min = -np.nanpercentile(np.abs(co_field),99.5); co_max = -co_min
                #%% plot results this far
                fig, (ax0,ax1,ax2) = plt.subplots(1,3, figsize = (16,5))
                btm_bound = gridSW[1]; top_bound = gridNE[1]
                left = gridSW[0]; right = gridNE[0]
                
                ax0.pcolormesh(local_lon, local_lat, pert_theta, cmap = 'seismic', vmin = pT_min, vmax = pT_max, **pcol_kwargs)
                ax0.contour(grid_lon, grid_lat,lbm_cube.data, colors = 'k')
                ax0.set_ylim(bottom = btm_bound, top = top_bound)
                ax0.set_xlim(left = left, right = right)
                
                ax1.pcolormesh(local_lon, local_lat, pert_w, cmap = 'seismic', vmin = pw_min, vmax = pw_max,**pcol_kwargs)
                ax1.contour(grid_lon, grid_lat,lbm_cube.data, colors = 'k')
                ax1.set_ylim(bottom = btm_bound, top = top_bound)
                ax1.set_xlim(left = left, right = right)
                
                ax2.pcolormesh(local_lon,local_lat, co_field, cmap = 'seismic', vmin = co_min, vmax = co_max, **pcol_kwargs)
                ax2.contour(grid_lon, grid_lat,lbm_cube.data, colors = 'k')
                ax2.set_ylim(bottom = btm_bound, top = top_bound)
                ax2.set_xlim(left = left, right = right)