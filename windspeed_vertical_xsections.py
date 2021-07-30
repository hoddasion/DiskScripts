# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 12:23:16 2020

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
#from iris.maths import add
#%% execute function
legs = [8]
resolutions = ['0p5km', '1p5km', '4p4km']
varname1 = 'm01s15i002'
varname2 = 'm01s15i003'
varname = 'windspeed'
vmin = 0; vmax = 20
cmap = 'Blues'
#%%
Test =True
if Test:
    D = '0p5km'; leg = 8
    Case = {'res' : D, 'flight' : 301, 'varname' : varname, 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
    res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
    fileleg = Case['fileleg']
    
    ucube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname1}/{config}_{res}_{varname1}_flt{flight}_leg{fileleg}.nc', varname1)
    vcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname2}/{config}_{res}_{varname2}_flt{flight}_leg{fileleg}.nc', varname2)
    maincube = (ucube**2 + vcube**2)**0.5
    thetacube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname1}/{config}_{res}_{varname1}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
    orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
    print(maincube, '\n', thetacube, '\n', orogcube)
    print(maincube.units)
    #%%
    fig = workshop.plot_vertical_xsection(maincube, thetacube, orogcube,index = 6, 
                                          Case = Case, cmap = cmap, 
                                          vmin = vmin,vmax = vmax, levels = 30, alt_limit = 2000,
                                          unit = r'ms$^{-1}$')
    plt.show()
#%%
run_loop = False
if run_loop:
    for D in resolutions:
        for leg in legs:
            try:
                #%%
            
                Case = {'res' : D, 'flight' : 301, 'varname' : varname, 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
                res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
                fileleg = Case['fileleg']
                print(Case.values())
                #%%
                
                ucube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname1}/{config}_{res}_{varname1}_flt{flight}_leg{fileleg}.nc', varname1)
                vcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname2}/{config}_{res}_{varname2}_flt{flight}_leg{fileleg}.nc', varname2)
                maincube = (ucube**2 + vcube**2)**0.5
                thetacube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname1}/{config}_{res}_{varname1}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
                orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
                print(maincube, '\n', thetacube, '\n', orogcube)
                print(maincube.units)
            except:
                print('File error')
                continue
            path = f'../../Figures/PNG/301/u-bu807/Vertical/{varname}/Alt2000/leg{fileleg}/{res}'
            for index in range(47):
                
        
                
                fig = workshop.plot_vertical_xsection(maincube, thetacube, orogcube,index, Case, cmap = cmap, vmin = vmin,vmax = vmax, levels = 30, alt_limit = 2000, unit = r'ms$^{-1}$')
                workshop.save_vert_xsection(fig, maincube,index, path, Case, fileleg = leg)
                #plt.close()
                plt.show()
                
obs_overlay = False
if obs_overlay:
    #%%
    for D in resolutions:
    
        leg = 8
        Case = {'res' : D, 'flight' : 301, 'varname' : varname, 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
        res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
        fileleg = Case['fileleg']
        
        #tcubes = iris.load(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname1}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc')[2]
        #print(tcubes)
        #%%
        ucube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname1}/{config}_{res}_{varname1}_flt{flight}_leg{fileleg}.nc', varname1)
        vcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname2}/{config}_{res}_{varname2}_flt{flight}_leg{fileleg}.nc', varname2)
        maincube = (ucube**2 + vcube**2)**0.5
        
        thetacube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname1}/{config}_{res}_{varname1}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
        orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
        print(maincube, '\n', orogcube)
        print(maincube.units)
        print('maximum',np.nanmax(maincube.data))
        compcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/RA1M_0p5km_umh1_flt{flight}.nc', varname1)
        polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
        polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
        print(polelat)
        #%%
        ## load quality controlled and interval-meaned observational data from Database
        # using pandas
        database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_60s_301.txt', delimiter = ' ')
        
        db_legno = database['legno']
        #legcond = np.array(db_legno == 3) + np.array(db_legno == 4) + np.array(db_legno == 5) + np.array(db_legno == 6)
        #legcond = db_legno == leg
        legcond = np.array(db_legno == 8) + np.array(db_legno == 9)
        db_u = np.array(database['wsp'][legcond]) 
        
       
        db_lon = np.array(database['lon'][legcond])
        db_lat = np.array(database['lat'][legcond])
        db_alt = np.array(database['altgps'][legcond])
        
        db_mtime = np.array(database['meantime'][legcond])
        
        #%% rotate db coords
    
        rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
        rot_db_lon = rot_db_lon + 360
        
        
        #%%
        t_indices = [31,32]#[28,29,30]
        
        for index in t_indices:
            fig = workshop.plot_vertical_xsection(maincube, thetacube, orogcube,
                                                  index = index, Case = Case, cmap = cmap, vmin = vmin,vmax = vmax, levels = 30, 
                                                  theta_contour = True,alt_limit = 2000, unit = r'ms$^{-1}$')
            plt.scatter(rot_db_lon, db_alt,c=db_u, s = 300, cmap = cmap, 
                        vmin = vmin,vmax = vmax, edgecolor = 'k')
            path = f'../../Figures/PNG/301/u-bu807/Vertical/{varname}/Alt2000/Obs/{res}'
            workshop.save_vert_xsection(fig, maincube,index, path, Case, fileleg = leg, obs = True, obs_detail = 'l89')
            #plt.close()
            plt.show()

#%%     
plot_series = False
if plot_series:
    ## plot a spatial-time series
    
    vmin = 0;vmax = 12
    
    
    D = '0p5km'; leg = 8
    tindex = 31
    hindex = 15
    plotleg = 8
    Case = {'res' : D, 'flight' : 301, 'varname' : varname, 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
    res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
    fileleg = Case['fileleg']
    
    
    ucube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname1}/{config}_{res}_{varname1}_flt{flight}_leg{fileleg}.nc', varname1)
    vcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname2}/{config}_{res}_{varname2}_flt{flight}_leg{fileleg}.nc', varname2)
    maincube = (ucube**2 + vcube**2)**0.5
    
    #thetacube = maincube#iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
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
    database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_60s_301.txt', delimiter = ' ')
    
    db_legno = database['legno']
    #legcond = np.array(db_legno == 3) + np.array(db_legno == 4) + np.array(db_legno == 5) + np.array(db_legno == 6)
    legcond = db_legno == plotleg
    db_w = np.array(database['wsp'][legcond]) # load and convert from g/kg to kg/kg
   
   
    db_lon = np.array(database['lon'][legcond])
    db_lat = np.array(database['lat'][legcond])
    db_alt = np.array(database['altgps'][legcond])
    #%% rotate db coords

    rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
    rot_db_lon = rot_db_lon + 360
    #%% load qc 1hz data
    import xarray as xr
    masin_file = xr.open_dataset('../../Obvs_Data/UEA_qc_MASIN_data/MASIN_flightdata_301.nc')
    print(masin_file)
    mas_u = np.array(masin_file['u_50hz'])[::50]
    mas_v = np.array(masin_file['v_50hz'])[::50]
    mas_time = np.array(masin_file['time_1hz'])
    mas_lon = np.array(masin_file['lat_50hz'])[::50]
    mas_lat = np.array(masin_file['lon_50hz'])[::50]  
    mas_wind = (mas_u**2 + mas_v**2)**0.5
    
    rot_mas_lon, rot_mas_lat = rotate_pole(mas_lon, mas_lat, polelon, polelat)
    rot_mas_lon = rot_mas_lon + 360
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
    host.plot(rot_mas_lon,db_w, label = 'Obs')
    host.set_ylim(bottom = vmin, top = vmax)
    #host.plot(rot_db_lon,db_q_licor, label = 'LICOR')
    plt.legend()
    par2.plot(oroglon,orogdata, color = 'k', linestyle = '--',alpha = 0.5)
    plt.xticks([])
    host.set_ylabel(r'ms$^{-1}$')
    par2.set_ylabel('Orographic altitude [m]')
    plt.title(f'Spatial series for {varname} T{timepoint} model mean alt {mean_alt} res {res} leg{plotleg} flt301')
    plt.grid(True)
    plt.draw()
    plt.show()
    #plt.savefig(f'../../Figures/PNG/301/{suite}/Series/{exp}/{varname}/{res}/{exp}_{res}_{varname}_spatial_series_obsl{plotleg}_{timepoint}_flt301.png')