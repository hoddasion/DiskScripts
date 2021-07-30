# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 15:05:01 2020

@author: kse18nru
"""

#%% Import modules and scripts

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import iris
import iris.coord_categorisation
import xsection_workshop as workshop
from iris.analysis.cartography import rotate_pole

#%% execute function
legs = [8]
resolutions = ['0p5km', '1p5km', '4p4km']
varname = 'air_potential_temperature'
vmin = 265; vmax = 285
#%%
test = False
if test:
    D = '0p5km'; leg = 8
    Case = {'res' : D, 'flight' : 301, 'varname' : varname, 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
    res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
    fileleg = Case['fileleg']
    
    cubes = iris.load(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc')
    maincube = cubes[0][:,:30]
    thetacube = maincube#iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
    orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
    print(maincube, '\n', orogcube)
    print(maincube.units)
    print('maximum',np.nanmax(maincube.data))
#%%
    fig = workshop.plot_vertical_xsection(maincube, thetacube, orogcube,index = 6, Case = Case, cmap = 'plasma', vmin = vmin,vmax = vmax, levels = 20, theta_contour = False,alt_limit = 2000)
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
                
                cubes = iris.load(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc'); maincube = cubes[0][:,:30]
                thetacube = maincube#iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
                orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
                #print(maincube, '\n', thetacube, '\n', orogcube)
                
            except:
                print('File error')
                continue
            path = f'../../Figures/PNG/301/u-bu807/Vertical/{varname}/Alt2000/leg{fileleg}/{res}'
            for index in range(47):
                
        
                
                fig = workshop.plot_vertical_xsection(maincube, thetacube, orogcube,index, Case, cmap = 'plasma', vmin = vmin,vmax = vmax, levels = 20, theta_contour = False,alt_limit = 2000)
                workshop.save_vert_xsection(fig, maincube,index, path, Case, fileleg = leg)
                #plt.close()
                plt.show()
                
#%%
obs_overlay = False
if obs_overlay:
    #%%
    D = '4p4km'; leg = 8
    Case = {'res' : D, 'flight' : 301, 'varname' : varname, 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
    res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
    fileleg = Case['fileleg']
    
    
    cubes = iris.load(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc')
    maincube = cubes[0][:,:30]
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
    #legcond = db_legno == leg
    legcond = np.array(db_legno == 8) + np.array(db_legno == 9)
    db_theta = np.array(database['theta'][legcond])
    
   
    db_lon = np.array(database['lon'][legcond])
    db_lat = np.array(database['lat'][legcond])
    db_alt = np.array(database['altgps'][legcond])
    
    db_mtime = np.array(database['meantime'][legcond])
    
    #%% rotate db coords

    rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
    rot_db_lon = rot_db_lon + 360
    
    
    #%%
    t_indices = [31,32]#[28,29,30]
    cmap = 'plasma'
    for index in t_indices:
        fig = workshop.plot_vertical_xsection(maincube, thetacube, orogcube,
                                              index = index, Case = Case, cmap = cmap, vmin = vmin,vmax = vmax, levels = 30, 
                                              theta_contour = False,alt_limit = 2000)
        plt.scatter(rot_db_lon, db_alt,c=db_theta, s = 300, cmap = cmap, 
                    vmin = vmin,vmax = vmax, edgecolor = 'k')
        path = f'../../Figures/PNG/301/u-bu807/Vertical/{varname}/Alt2000/Obs/{res}'
        workshop.save_vert_xsection(fig, maincube,index, path, Case, fileleg = leg, obs = True, obs_detail = 'l89')
        #plt.close()
        plt.show()
        
#%%     
plot_series = False
if plot_series:
    ## plot a spatial-time series
    
   
    varname = 'air_potential_temperature'
    
    D = '0p5km'; leg = 8
    tindex = 31
    hindex = 16
    plotleg = 8
    Case = {'res' : D, 'flight' : 301, 'varname' : varname, 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
    res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
    fileleg = Case['fileleg']
    
    
    maincube = iris.load(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc')[0][:,:30]
    
    #thetacube = maincube#iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
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
    db_theta = np.array(database['theta'][legcond]) # load and convert from g/kg to kg/kg
    
   
    db_lon = np.array(database['lon'][legcond])
    db_lat = np.array(database['lat'][legcond])
    db_alt = np.array(database['altgps'][legcond])
    #%% rotate db coords

    rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
    rot_db_lon = rot_db_lon + 360
    
  
    
    
    
    #%%
    
    yb = vmin; yt = vmax
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
    host.plot(rot_db_lon,db_theta, label = 'Obs')
    host.set_ylim(bottom = yb, top = yt)
    plt.legend()
    par2.plot(oroglon,orogdata, color = 'k', linestyle = '--',alpha = 0.5)
    plt.xticks([])
    host.set_ylabel(maincube.units)
    par2.set_ylabel('Orographic altitude [m]')
    plt.title(f'Spatial series for {varname} T{timepoint} model mean alt {mean_alt} res {res} leg{plotleg} flt301')
    plt.grid(True)
    plt.draw()
    plt.savefig(f'../../Figures/PNG/301/{suite}/Series/{exp}/{varname}/{res}/{exp}_{res}_{varname}_spatial_series_obsl{plotleg}_{timepoint}_flt301.png')
