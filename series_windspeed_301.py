#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 12:39:00 2020

@author: kse18nru
"""

#%% Import modules and scripts


import matplotlib.pyplot as plt
from iris.analysis.cartography import rotate_pole
import numpy as np
import iris
import iris.coord_categorisation
import pandas as pd
import xsection_workshop as workshop
import persinterpolate as pert

#%%
resolutions = ['0p5km', '1p5km', '4p4km']
varname1 = 'm01s15i002'
varname2 = 'm01s15i003'
varname = 'windspeed'
vmin = 0; vmax = 16
cmap = 'Blues'
leg = 8
plotleg = 8
leg_meta = ((1,26,1),(2,27,2),(3,28,3),(4,29,3),(5,30,3),(6,30,3),(8,31,8),(9,32,8))
for D in resolutions:
    for plotleg in leg_meta:
        time_index = plotleg[1]
        Case = {'res' : D, 'flight' : 301, 'varname' : varname, 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
        res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
        fileleg = Case['fileleg']
        #%% load model data
        ucube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname1}/{config}_{res}_{varname1}_flt{flight}_leg{fileleg}.nc', varname1)
        vcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname2}/{config}_{res}_{varname2}_flt{flight}_leg{fileleg}.nc', varname2)
        maincube = (ucube**2 + vcube**2)**0.5
        orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
        print(maincube, '\n', orogcube)
        print(maincube.units)
        print('maximum',np.nanmax(maincube.data))
        compcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/RA1M_0p5km_umh1_flt{flight}.nc', varname1)
        polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
        polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
        print(polelat)
        
        orogmax = orogcube.collapsed('grid_latitude', iris.analysis.MAX)
        orogdata = orogmax.data
        oroglon = orogmax.coord('grid_longitude').points
        
        #%% load quality controlled and interval-meaned observational data from Database
        # using pandas
        database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_sr_301.txt', delimiter = ' ')
        
        db_legno = database['legno']
        #legcond = np.array(db_legno == 3) + np.array(db_legno == 4) + np.array(db_legno == 5) + np.array(db_legno == 6)
        legcond = db_legno == plotleg[0]
        db_w = np.array(database['wsp'][legcond]) # load and convert from g/kg to kg/kg
       
       
        db_lon = np.array(database['lon'][legcond])
        db_lat = np.array(database['lat'][legcond])
        db_alt = np.array(database['altgps'][legcond])
        
        #%% rotate db coords
    
        rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
        rot_db_lon = rot_db_lon + 360
        
        #%% interpolate model data onto obs altitudes
        
        matched_series, matched_coors = pert.interpolate_matched_series(maincube, db_alt, rot_db_lon, time_index)
        simple_series, simple_coors = pert.interpolate_simple_series(maincube, np.nanmean(db_alt), time_index)
        
        #%%
        
        
        
        timecoord = maincube.coord('time')
        datetimes = timecoord.units.num2date(timecoord.points)
        timepoint = datetimes[time_index].strftime('H%HM%M')
        
        from mpl_toolkits.axes_grid1 import host_subplot
        import mpl_toolkits.axisartist as AA
        
        fig = plt.figure(figsize = (20,8))
        host = host_subplot(111, axes_class=AA.Axes)
        plt.subplots_adjust(right=0.75)
        
        
        par2 = host.twinx()
        
        offset = 0
        new_fixed_axis = par2.get_grid_helper().new_fixed_axis
        par2.axis["right"] = new_fixed_axis(loc="right",
                                            axes=par2,
                                            offset=(offset, 0))
        
        par2.axis["right"].toggle(all=True)
        
        host.plot(simple_coors, simple_series, label = 'UM - Control')
        host.plot(rot_db_lon,db_w, label = 'Obs')
        host.set_ylim(bottom = vmin, top = vmax)
        #host.plot(rot_db_lon,db_q_licor, label = 'LICOR')
        plt.legend()
        par2.plot(oroglon,orogdata, color = 'k', linestyle = '--',alpha = 0.5)
        plt.xticks([])
        host.set_ylabel(r'ms$^{-1}$')
        par2.set_ylabel('Orographic altitude [m]')
        plt.title(f'Spatial series for {varname} T{timepoint} model res {res} leg{plotleg[0]} flt301')
        plt.grid(True)
        plt.draw()
        
        plt.savefig(f'../../Figures/PNG/301/{suite}/Series/{exp}/{varname}/{res}/{exp}_{res}_{varname}_spatial_series_obsl{plotleg[0]}_{timepoint}_flt301.png')
        
        
        #%% plot scatter
        plot_scatter = True
        if plot_scatter:
            fig, ax = plt.subplots(1,1, figsize = (16,16))
            ax.scatter(db_w, matched_series)
            ax.set_ylabel(r'Model windspeed, ms$^{-1}$')
            ax.set_xlabel(r'Observed windspeed, ms$^{-1}$')
            ax.grid(True)
            ax.set_ylim(bottom = vmin, top = vmax)
            ax.set_xlim(left = vmin, right = vmax)
            plt.title(f'Matched point amplitutdes for {varname} T{timepoint} model res {res} leg{plotleg[0]} flt301')
            
            plt.savefig(f'../../Figures/PNG/301/{suite}/Scatter/{exp}/{varname}/{res}/{exp}_{res}_{varname}_amplitude_scatter_obsl{plotleg[0]}_{timepoint}_flt301.png')