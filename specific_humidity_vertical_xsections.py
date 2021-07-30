# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 14:45:54 2020

@author: kse18nru
"""

#%% Import modules and scripts
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import iris
import iris.coord_categorisation
import iris.plot as iplt
import model_foundry as foundry
import pandas as pd
from iris.analysis.cartography import rotate_pole
from iris.analysis import trajectory
from model_foundry import toolkit
import datetime
import xsection_workshop as workshop

#%% execute function
Test = False
if Test:
    legs = [1,2,3]
    resolutions = ['0p5km', '1p5km', '4p4km']
    varname = 'specific_humidity'
    #%%
    D = '0p5km'; leg = 1
    Case = {'res' : D, 'flight' : 301, 'varname' : varname, 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
    res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
    fileleg = Case['fileleg']
    
    maincube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', varname)
    thetacube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
    orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
    print(maincube, '\n', thetacube, '\n', orogcube)
    print(maincube.units)
    print('maximum',np.nanmax(maincube.data))
    #%%
    fig = workshop.plot_vertical_xsection(maincube, thetacube, orogcube,index = 6, Case = Case, cmap = 'cividis', vmin = 0,vmax = 0.004, levels = 20)
    plt.show()
#%%
run_loop = False
if run_loop == True:
    for D in resolutions:
        for leg in legs:
            try:
                #%%
            
                Case = {'res' : D, 'flight' : 301, 'varname' : varname, 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
                res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
                fileleg = Case['fileleg']
                print(Case.values())
                #%%
                
                maincube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', varname)
                thetacube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
                orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
                #print(maincube, '\n', thetacube, '\n', orogcube)
                
            except:
                print('File error')
                continue
            path = f'../../Figures/PNG/301/u-bu807/Vertical/{varname}/leg{fileleg}/{res}'
            for index in range(47):
                
        
                
                fig = workshop.plot_vertical_xsection(maincube, thetacube, orogcube,index, Case, cmap = 'cividis', vmin = 0,vmax = 0.004, levels = 20)
                workshop.save_vert_xsection(fig, maincube,index, path, Case, fileleg = leg)
                #plt.close()
                plt.show()
                
obs_overlay = False
if obs_overlay:
    #%%
    varname = 'specific_humidity'
    vmin = 0;vmax = 0.004
    D = '0p5km'; leg = 2
    Case = {'res' : D, 'flight' : 301, 'varname' : varname, 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
    res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
    fileleg = Case['fileleg']
    
    
    maincube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc',varname)[:,:30]
    
    thetacube = maincube#iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
    orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
    print(maincube, '\n', orogcube)
    print(maincube.units)
    print('maximum',np.nanmax(maincube.data))
    compcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/RA1M_0p5km_umi1_flt{flight}.nc', varname)
    polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    print(polelat)
    #%%
    ## load quality controlled and interval-meaned observational data from Database
    # using pandas
    database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_60s_301.txt', delimiter = ' ')
    
    db_legno = database['legno']
    #legcond = np.array(db_legno == 3) + np.array(db_legno == 4) + np.array(db_legno == 5) + np.array(db_legno == 6)
    legcond = db_legno == leg
    db_q = np.array(database['q_LICOR'][legcond])/1000 # load and convert from g/kg to kg/kg
    
   
    db_lon = np.array(database['lon'][legcond])
    db_lat = np.array(database['lat'][legcond])
    db_alt = np.array(database['altgps'][legcond])
    
    db_mtime = np.array(database['meantime'][legcond])
    
    #%% rotate db coords

    rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
    rot_db_lon = rot_db_lon + 360
    
    
    #%%
    t_indices = [27]#[28,29,30]
    cmap = 'cividis'
    for index in t_indices:
        fig = workshop.plot_vertical_xsection(maincube, thetacube, orogcube,
                                              index = index, Case = Case, cmap = cmap, vmin = vmin,vmax = vmax, levels = 20, 
                                              theta_contour = False,alt_limit = 2000)
        plt.scatter(rot_db_lon, db_alt,c=db_q, s = 300, cmap = cmap, 
                    vmin = vmin,vmax = vmax, edgecolor = 'k')
        path = f'../../Figures/PNG/301/u-bu807/Vertical/{varname}/Alt2000/Obs/{res}'
        workshop.save_vert_xsection(fig, maincube,index, path, Case, fileleg = leg, obs = True, obs_detail = 'l2_LICOR')
        #plt.close()
        plt.show()
#%%     
plot_series =True
if plot_series:
    ## plot a spatial-time series
    
   
    varname = 'specific_humidity'
    vmin = 0;vmax = 0.004
    D = '4p4km'; leg = 3
    tindex = 30
    hindex = 20
    plotleg = 6
    Case = {'res' : D, 'flight' : 301, 'varname' : varname, 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
    res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
    fileleg = Case['fileleg']
    
    
    maincube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc',varname)[:,:30]
    
    thetacube = maincube#iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
    orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
    print(maincube, '\n', orogcube)
    print(maincube.units)
    print('maximum',np.nanmax(maincube.data))
    compcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/RA1M_0p5km_umi1_flt{flight}.nc', varname)
    polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
    polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
    print(polelat)
    
    orogmax = orogcube.collapsed('grid_latitude', iris.analysis.MAX)
    orogdata = orogmax.data
    oroglon = orogmax.coord('grid_longitude').points
   ## load quality controlled and interval-meaned observational data from Database
    # using pandas
    database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_60s_301.txt', delimiter = ' ')
    
    db_legno = database['legno']
    #legcond = np.array(db_legno == 3) + np.array(db_legno == 4) + np.array(db_legno == 5) + np.array(db_legno == 6)
    legcond = db_legno == plotleg
    db_q_buck = np.array(database['q_BUCK'][legcond])/1000 # load and convert from g/kg to kg/kg
    db_q_licor = np.array(database['q_LICOR'][legcond])/1000
   
    db_lon = np.array(database['lon'][legcond])
    db_lat = np.array(database['lat'][legcond])
    db_alt = np.array(database['altgps'][legcond])
    #%% rotate db coords

    rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
    rot_db_lon = rot_db_lon + 360
    
    #%%
    ## do some plotting
    # print(maincube)
    
    # for i in range(30):
    #     print(i,np.mean(maincube.coord('altitude').points, axis = 1)[i],'\n' )
    # plt.plot(np.arange(1,31),np.mean(maincube.coord('altitude').points, axis = 1))
    # xrow = rot_db_lon
    # yrow = rot_db_lat
    # zrow = db_alt
    # sample_points = [('altitude', zrow),('grid_longitude',xrow),('grid_latitude',yrow)]
    #linecube = trajectory.interpolate(maincube,sample_points, method = 'nearest')
    #print(linecube)
    
    
    
    #%%

    mlon = maincube[tindex,hindex].coord('grid_longitude').points
    mdata = maincube[tindex,hindex].data
    timecoord = maincube.coord('time')
    datetimes = timecoord.units.num2date(timecoord.points)
    timepoint = datetimes[tindex].strftime('H%HM%M')
    mean_alt = np.mean(maincube.coord('altitude').points, axis = 1)[hindex]//1
    from mpl_toolkits.axes_grid1 import host_subplot
    import mpl_toolkits.axisartist as AA
    import matplotlib.pyplot as plt
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
    
    host.plot(mlon,mdata, label = 'UM - Control')
    host.plot(rot_db_lon,db_q_buck, label = 'BUCK')
    host.plot(rot_db_lon,db_q_licor, label = 'LICOR')
    host.set_ylim(bottom = vmin, top =vmax+0.001)
    host.legend()
    par2.plot(oroglon,orogdata, color = 'k', linestyle = '--',alpha = 0.5)
    plt.xticks([])
    host.set_ylabel(maincube.units)
    par2.set_ylabel('Orographic altitude [m]')
    plt.title(f'Spatial series for {varname} T{timepoint} model mean alt {mean_alt} res {res} leg{plotleg} flt301')
    host.grid(True)
    plt.draw()
    #plt.show()
    fig.savefig(f'../../Figures/PNG/301/{suite}/Series/{exp}/{varname}/{res}/{exp}_{res}_{varname}_spatial_series_obsl{plotleg}_{timepoint}_flt301.png')