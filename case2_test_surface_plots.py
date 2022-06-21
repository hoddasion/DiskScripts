# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 15:31:20 2022

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
experiment = 'LONGTAIL'
config = 'LONGTAIL'
suite = 'u-cf117'

#%% thermodynamic and scalar flux analysis
thermo_scalar = True
if thermo_scalar:
    ### load control data
    control_path = 'D:/Project/Model_Data/u-cc134/'
    ## atmospheric variables
    CONTR_TEMP = iris.load_cube(f'{control_path}RA1M_0p5km_um_air_temperature_24hrs_pg_306.nc', 'air_temperature')
    CONTR_REL = iris.load_cube(f'{control_path}RA1M_0p5km_um_relative_humidity_24hrs_pg_306.nc', 'relative_humidity')
    CONTR_MSP = iris.load_cube(f'{control_path}RA1M_0p5km_um_air_pressure_at_sea_level_24hrs_pg_306.nc', 'air_pressure_at_sea_level')
    ## scalar flux variables
    CONTR_SH = iris.load_cube(f'{control_path}RA1M_0p5km_um_surface_upward_sensible_heat_flux_24hrs_pg_306.nc', 'surface_upward_sensible_heat_flux')
    CONTR_LH = iris.load_cube(f'{control_path}RA1M_0p5km_um_surface_upward_latent_heat_flux_24hrs_pg_306.nc', 'surface_upward_latent_heat_flux')
    ## geographic variables
    CONTR_LBM = iris.load_cube(f'{control_path}RA1M_0p5km_um_land_binary_mask_24hrs_pi_306.nc', 'land_binary_mask')
    ## windspeed
    CONTR_XWIND = iris.load_cube(f'{control_path}RA1M_0p5km_um_x_wind_24hrs_pg_306.nc', 'x_wind')
    CONTR_YWIND = iris.load_cube(f'{control_path}RA1M_0p5km_um_y_wind_24hrs_pg_306.nc', 'y_wind')
    CONTR_WSP = iris.load_cube(f'{control_path}RA1M_0p5km_um_x_wind_24hrs_pg_306.nc', 'x_wind')
    CONTR_WSP.data = (CONTR_XWIND.data**2 + CONTR_YWIND.data**2)**0.5
    
    ### load test data
    test_path = f'D:/Project/Model_Data/{suite}/'
    TEST_TEMP = iris.load_cube(f'{test_path}{config}_0p5km_um_air_temperature_24hrs_pg_306.nc', 'air_temperature')
    TEST_REL = iris.load_cube(f'{test_path}{config}_0p5km_um_relative_humidity_24hrs_pg_306.nc', 'relative_humidity')
    TEST_MSP = iris.load_cube(f'{test_path}{config}_0p5km_um_air_pressure_at_sea_level_24hrs_pg_306.nc', 'air_pressure_at_sea_level')
    ## scalar flux variables
    TEST_SH = iris.load_cube(f'{test_path}{config}_0p5km_um_surface_upward_sensible_heat_flux_24hrs_pg_306.nc', 'surface_upward_sensible_heat_flux')
    TEST_LH = iris.load_cube(f'{test_path}{config}_0p5km_um_surface_upward_latent_heat_flux_24hrs_pg_306.nc', 'surface_upward_latent_heat_flux')
    ## windspeed
    TEST_XWIND = iris.load_cube(f'{test_path}{config}_0p5km_um_x_wind_24hrs_pg_306.nc', 'x_wind')
    TEST_YWIND = iris.load_cube(f'{test_path}{config}_0p5km_um_y_wind_24hrs_pg_306.nc', 'y_wind')
    TEST_WSP = iris.load_cube(f'{test_path}{config}_0p5km_um_x_wind_24hrs_pg_306.nc', 'x_wind')
    TEST_WSP.data = (TEST_XWIND.data**2 + TEST_YWIND.data**2)**0.5
    
    ### extract grid lat/lon
    grid_lat = CONTR_TEMP.coords('grid_latitude')[0].points
    grid_lon = CONTR_TEMP.coords('grid_longitude')[0].points
    print(TEST_SH)
    
    #%% optional masking of land-area
    masked_land = True
    if masked_land:
        CONTR_TEMP.data[:,CONTR_LBM.data == 1] = np.nan
        CONTR_REL.data[:,CONTR_LBM.data == 1] = np.nan
        CONTR_MSP.data[:,CONTR_LBM.data == 1] = np.nan
        CONTR_SH.data[:,CONTR_LBM.data == 1] = np.nan
        CONTR_LH.data[:,CONTR_LBM.data == 1] = np.nan
        CONTR_WSP.data[:,CONTR_LBM.data == 1] = np.nan
        
        
    #%% format date time for file name and figure title
    time = CONTR_REL.coord('time')
    dates = time.units.num2date(time.points)
    
    #%%compute anomalies (as extracted ndarrays)
    
    ANO_TEMP = TEST_TEMP.data - CONTR_TEMP.data
    ANO_REL = TEST_REL.data - CONTR_REL.data
    ANO_SH = TEST_SH.data - CONTR_SH.data
    ANO_LH = TEST_LH.data - CONTR_LH.data
    ANO_WSP = TEST_WSP.data - CONTR_WSP.data
    #%% normalisation
    ## symmetrics first
    ANO_Tnorm, ANO_Tmap = workshop.mplt_symmetric_normalisation(ANO_TEMP, cmap = 'seismic')
    ANO_RELnorm, ANO_RELmap = workshop.mplt_symmetric_normalisation(ANO_TEMP, cmap = 'seismic')
    
    ## normal quantities now
    Tnorm = matplotlib.colors.Normalize(vmin = np.nanmin(CONTR_TEMP.data),vmax = np.nanmax(CONTR_TEMP.data))
    Tmap = matplotlib.cm.ScalarMappable(norm = Tnorm)
    RELnorm = matplotlib.colors.Normalize(vmin =0,vmax = np.nanmax(CONTR_REL.data))
    RELmap = matplotlib.cm.ScalarMappable(norm = RELnorm, cmap = 'Blues')
    #%%  plotting
    
    
    temp_rh_anomaly = True
    if temp_rh_anomaly:
        for tidx in range(48):
            
            com_cbarkeys = {'pad' : 0.005, 'orientation':'horizontal'}
            fig, (ax0,ax1,ax2,ax3) = plt.subplots(1,4, figsize = (18,8))
            
            ax0.pcolormesh(grid_lon, grid_lat, CONTR_TEMP.data[tidx], norm = Tnorm)
            T_cbar = plt.colorbar(Tmap, ax = ax0, **com_cbarkeys)
            T_cbar.ax.set(xlabel = 'K')
            ax0.set_title('1.5m Temp - Control')
            
            ax1.pcolormesh(grid_lon,grid_lat, ANO_TEMP[tidx], norm = ANO_Tnorm, cmap = 'seismic')
            ANO_T_cbar = plt.colorbar(ANO_Tmap, ax = ax1, **com_cbarkeys)
            ANO_T_cbar.ax.set(xlabel = 'K')
            ax1.set_title('1.5m Temp - Anomaly')
            
            ax2.pcolormesh(grid_lon,grid_lat, CONTR_REL.data[tidx], cmap = 'Blues', norm = RELnorm)
            REL_cbar = plt.colorbar(RELmap, ax = ax2, **com_cbarkeys)
            REL_cbar.ax.set(xlabel = '%')
            ax2.set_title('1.5m Rel. hum. - Control')
            
            ax3.pcolormesh(grid_lon, grid_lat, ANO_REL.data[tidx], cmap = 'seismic', norm = ANO_RELnorm)
            ANO_REL_cbar = plt.colorbar(ANO_RELmap, ax = ax3, **com_cbarkeys)
            ANO_REL_cbar.ax.set(xlabel = '%')
            ax3.set_title('1.5m Rel. hum. - Anomaly')
            
            for ax in (ax0,ax1,ax2,ax3):
                ax.contour(grid_lon, grid_lat, CONTR_LBM.data, colors = 'k')
                ax.set_ylim(bottom = -1.1, top = 1.5)
                ax.set_xticks([])
                ax.set_yticks([])
                
            timepoint = dates[tidx].strftime('H%HM%M')
            fig.suptitle(f'Case 2 - 0p5km CONTROL/{experiment} - {timepoint}')
            plt.tight_layout()
            if True:
                if masked_land:
                    plt.savefig(f'D:/Project/Figures/PNG/306/{experiment}/{suite}/surface/temp_rh_anomaly/masked/0p5km/{config}_0p5km_temp_rh_anomaly_surface_{timepoint}_306_masked.png')
                else:
                    plt.savefig(f'D:/Project/Figures/PNG/306/{experiment}/{suite}/surface/temp_rh_anomaly/0p5km/{config}_0p5km_temp_rh_anomaly_surface_{timepoint}_306.png')
                plt.close()
    
    
    #%% normalisation
    ## symmetrics first
    ANO_SHnorm, ANO_SHmap = workshop.mplt_symmetric_normalisation(ANO_SH, cmap = 'seismic')
    ANO_LHnorm, ANO_LHmap = workshop.mplt_symmetric_normalisation(ANO_LH, cmap = 'seismic')
    
    ## normal quantities now
    SHnorm, SHmap = workshop.mplt_symmetric_normalisation(CONTR_SH.data, cmap = 'seismic')
    LHnorm, LHmap = workshop.mplt_symmetric_normalisation(CONTR_LH.data, cmap = 'seismic')
    
    #%%  plotting
    
    
    sh_lh_anomaly = True
    if sh_lh_anomaly:
        for tidx in range(48):
            com_cbarkeys = {'pad' : 0.005, 'orientation':'horizontal'}
            fig, (ax0,ax1,ax2,ax3) = plt.subplots(1,4, figsize = (18,8))
            
            ax0.pcolormesh(grid_lon, grid_lat, CONTR_SH.data[tidx], norm = SHnorm, cmap = 'seismic')
            SH_cbar = plt.colorbar(SHmap, ax = ax0, **com_cbarkeys)
            SH_cbar.ax.set(xlabel = r'Wm$^{-2}$')
            ax0.set_title('SH - Control')
            
            ax1.pcolormesh(grid_lon,grid_lat, ANO_SH[tidx], norm = ANO_SHnorm, cmap = 'seismic')
            ANO_SH_cbar = plt.colorbar(ANO_SHmap, ax = ax1, **com_cbarkeys)
            ANO_SH_cbar.ax.set(xlabel =  r'Wm$^{-2}$')
            ax1.set_title('SH - Anomaly')
            
            ax2.pcolormesh(grid_lon,grid_lat, CONTR_LH.data[tidx], cmap = 'seismic', norm = LHnorm)
            LH_cbar = plt.colorbar(LHmap, ax = ax2, **com_cbarkeys)
            LH_cbar.ax.set(xlabel =  r'Wm$^{-2}$')
            ax2.set_title('LH - Control')
            
            ax3.pcolormesh(grid_lon, grid_lat, ANO_REL.data[tidx], cmap = 'seismic', norm = ANO_LHnorm)
            ANO_LH_cbar = plt.colorbar(ANO_LHmap, ax = ax3, **com_cbarkeys)
            ANO_LH_cbar.ax.set(xlabel =  r'Wm$^{-2}$')
            ax3.set_title('LH - Anomaly')
            
            for ax in (ax0,ax1,ax2,ax3):
                ax.contour(grid_lon, grid_lat, CONTR_LBM.data, colors = 'k')
                ax.set_ylim(bottom = -1.1, top = 1.5)
                ax.set_xticks([])
                ax.set_yticks([])
                
            timepoint = dates[tidx].strftime('H%HM%M')
            fig.suptitle(f'Case 2 - 0p5km CONTROL/{experiment} - {timepoint}')
            plt.tight_layout()
            if True:
                if masked_land:
                    plt.savefig(f'D:/Project/Figures/PNG/306/{experiment}/{suite}/surface/sh_lh_anomaly/masked/0p5km/{config}_0p5km_sh_lh_anomaly_surface_{timepoint}_306_masked.png')
                else:
                    plt.savefig(f'D:/Project/Figures/PNG/306/{experiment}/{suite}/surface/sh_lh_anomaly/0p5km/{config}_0p5km_sh_lh_anomaly_surface_{timepoint}_306.png')
                plt.close()
    