"""
To compute cross-sections and resolved maps of wind fields from model output data.

Created on 07/10/2019 by Wilhelm Hodder
"""

### module imports
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import iris
import iris.coord_categorisation
#import CubeCrossSectioner_UK as ccs
import iris.quickplot as qplt
import cartopy
import cartopy.feature as cfeat
import cartopy.crs as ccrs

#%% global presets
font = {'size' : 22}
matplotlib.rc('font', **font)
matplotlib.rcParams.update({'figure.max_open_warning':False})

#%% File loading

### pre-define path to data files
# variable naming convention: datapath_[res]_[file-suffix with redundant (2) zeros removed]
#datapath_1p5_pj0 = '../cylc-run/u-bk574/share/cycle/20180312T0000Z/IcelandGreenlandSeas/1p5km/RA1M/um/umnsaa_pj000' 
#datapath_1p5_cb1 = '../cylc-run/u-bk574/share/cycle/20180312T0000Z/IcelandGreenlandSeas/1p5km/RA1M/um/umnsaa_cb001'
#datapath_1p5_cb2 = '../cylc-run/u-bk574/share/cycle/20180312T0000Z/IcelandGreenlandSeas/1p5km/RA1M/um/umnsaa_cb002'

#datapath_4p4_ph06 = '../Model_Data/u-bk574/20180312T0000Z_IcelandGreenlandSeas_4p4km_RA1M_ph006.pp'
datapath_4p4_pg06 = '../Model_Data/u-bk574/nc/20180312T0000Z_IcelandGreenlandSeas_4p4km_RA1M_pg006.nc'
#datapath_4p4_pi06 = '../Model_Data/u-bk574/20180312T0000Z_IcelandGreenlandSeas_4p4km_RA1M_pi006.pp'
datapath_4p4_pj06 = '../Model_Data/u-bk574/nc/20180312T0000Z_IcelandGreenlandSeas_4p4km_RA1M_pj006.nc'
datapath_4p4_pf00 = '../Model_data/u-bk574/nc/20180312T0000Z_IcelandGreenlandSeas_4p4km_RA1M_pf000.nc'
#%%
datapath_4p4_pj = ['../Model_Data/u-bk574/nc/20180312T0000Z_IcelandGreenlandSeas_4p4km_RA1M_pj000.nc',
                   '../Model_Data/u-bk574/nc/20180312T0000Z_IcelandGreenlandSeas_4p4km_RA1M_pj006.nc',
                   '../Model_Data/u-bk574/nc/20180312T0000Z_IcelandGreenlandSeas_4p4km_RA1M_pj012.nc',
                   '../Model_Data/u-bk574/nc/20180312T0000Z_IcelandGreenlandSeas_4p4km_RA1M_pj018.nc']


#%%
### load cubes from file
cont_4p4_pf_height = iris.load(datapath_4p4_pf00, 'Height')
cont_4p4_pj06_xwind = iris.load(datapath_4p4_pj06, 'x wind component (with respect to grid)')
cont_4p4_pj06_ywind = iris.load(datapath_4p4_pj06, 'y wind component (with respect to grid)')
print(cont_4p4_pj06_ywind)
cont_4p4_pj_xwind = iris.cube.CubeList([])
cont_4p4_pj_ywind = iris.cube.CubeList([])
for i in range(4):
    # make list of cubes
    cont_4p4_pj_xwind.append(iris.load(datapath_4p4_pj[i], 'x wind component (with respect to grid)'))
    cont_4p4_pj_ywind.append(iris.load(datapath_4p4_pj[i], 'y wind component (with respect to grid)'))
cont_4p4_pj_xwind = cont_4p4_pj_xwind.concatenate()

#%%
## create list of pressure level values
p_levels = [950,850,800,500,300,200] # units in mbar
# extract cubes
cube_xwind = cont_4p4_pj06_xwind[0]
cube_ywind = cont_4p4_pj06_ywind[0]
cube_topo = cont_4p4_pf_height[0]
print(cube_ywind[0][0][0][::3])
#print(cube_ywind[2][0][99][0])

## rotate pole back
ywind_coords = cube_ywind[0][0]
print(ywind_coords)
#iris.analysis.cartography.rotate_pole(cube_ywind('longitude'), cube_ywind('latitude'), 157.2, 24.2)


plt.ion() # for interactive plotting


## prepare normalised wind data for quiver plot
topo_xcoor = cube_topo.coord('longitude').points
topo_ycoor = cube_topo.coord('latitude').points
x_coor = cube_xwind.coord('longitude').points # load coordinates
y_coor = cube_xwind.coord('latitude').points
windmagn = (cube_xwind ** 2 + cube_ywind ** 2) ** 0.5 # create cube of wind magnitude
windmagn.rename('windspeed')
print(windmagn)
u_wind = cube_xwind.data
v_wind = cube_ywind.data

u_norm = u_wind / np.sqrt(u_wind ** 2 + v_wind ** 2)
v_norm = v_wind / np.sqrt(u_wind ** 2 + v_wind ** 2)

## translate coordinates into rough distance scale
dist_lon_dom = (x_coor-np.nanmin(x_coor))/0.04 * 4.444 
dist_lat_dom = (y_coor-np.nanmin(y_coor))/0.04 * 4.444
dist_lat_topo = (topo_ycoor - np.nanmin(topo_ycoor))/0.04 * 4.444
for t in range(len(cube_xwind.coord('t').points)):
    for p in range(len(cube_xwind.coord('p').points)):
        ## detailed plot
        fig0, ax0 = plt.subplots(1,1, figsize = (17,15))
        n = 11 # <11 and the plotting messes up, I don't know why
        u_n = u_wind[t][p][::n]#; u_n = u_n[0:len(y_coor) - 1:n]
        v_n = v_wind[t][p][::n]#; v_n = v_n[0:len(y_coor) - 1:n]
        
        u_n2 = []
        v_n2 = []
        for j in range(len(u_n)):
            u_n2.append(u_n[j][0::n])
            v_n2.append(v_n[j][0::n])
        
        q1 = ax0.contourf(dist_lon_dom, dist_lat_dom, windmagn[t][p].data)#, cmap = plt.cm.YlGn)
        q2 = ax0.quiver(dist_lon_dom[::n], dist_lat_dom[:-1:n], u_n2, v_n2, pivot = 'middle') #u_wind[0][0][::n][0:len(y_coor)-1:n], v_wind[0][0][::n][0:len(y_coor)-1:n], pivot = 'middle')
        q3 = ax0.contour(dist_lon_dom, dist_lat_topo, cube_topo[0][0].data, cmap = plt.cm.gray)
        #pos0 = ax0.imshow(windmagn[t][p].data)
        cbar0 = fig0.colorbar(q1, ax = ax0)
        cbar0.ax.set(xlabel = r'U [ms$^{-1}$]')
        ax0.set(title = f'Wind speed and direction at {p_levels[p]}hPa, {t + 6}hrs', 
                                                       ylabel = 'km', xlabel = 'km')
        fig0.savefig(f'../Figures/PNG/P_levels/{p_levels[p]}hPa/Wind/RA1M_4p4km_301_windmagn_fulldomain_time{6 + t}_plevel{p_levels[p]}.png')
        fig0.savefig(f'../Figures/PDF/P_levels/{p_levels[p]}hPa/Wind/RA1M_4p4km_301_windmagn_fulldomain_time{6 + t}_plevel{p_levels[p]}.pdf')

### subset 4p4km data to focus on Snaesfellness Peninsula and Westfjords
N_limit = 200; S_limit = 120
W_limit = 120; E_limit = 200        

#sub_u = sub_cube_u.data
#sub_v = sub_cube_v.data


# prepare coordinate arrays of sub domain
sub_x_coor = x_coor[W_limit:E_limit] # load coordinates
sub_y_coor = y_coor[S_limit:N_limit]
sub_topo_xcoor = topo_xcoor[W_limit:E_limit]
sub_topo_ycoor = topo_ycoor[S_limit:N_limit]
domain_topo = cube_topo[0][0].data 
## rotate longitude and latotude back to original
print('x',len(sub_x_coor),'y',len(sub_y_coor))
Pole_lon = 223; Pole_lat = -47
rot_lon, rot_lat = iris.analysis.cartography.rotate_pole(sub_x_coor, sub_y_coor, Pole_lon, Pole_lat) # taking negatives of pole coordinates to rotate back to original
## convert coordinates into rough distance (km)
dist_lon = (sub_x_coor-np.nanmin(sub_x_coor))/0.04 * 4.444 
dist_lat = (sub_y_coor-np.nanmin(sub_y_coor))/0.04 * 4.444


for t in range(len(cube_xwind.coord('t').points)):
    for p in range(len(cube_xwind.coord('p').points)):
        
        ## prepare constrained data
        
        domain_windmagn = windmagn[t][p].data # extract data arrays
        
       
        temp_m = domain_windmagn[S_limit:N_limit] # temporary array
        u_wind_loc = u_wind[t][p]# t-p local wind data
        v_wind_loc = v_wind[t][p]
        temp_u = u_wind_loc[S_limit:N_limit]
        temp_v = v_wind_loc[S_limit:N_limit]
        temp_topo = domain_topo[S_limit:N_limit]
        sub_wind = []
        sub_u = []
        sub_v = []
        sub_topo = []
        for i in range(len(temp_m)):
            temp_ms = temp_m[i] # single out rows
            temp_us = temp_u[i]
            temp_vs = temp_v[i]
            temp_topos = temp_topo[i]
            sub_wind.append(temp_ms[W_limit:E_limit]) # append delimited rows
            sub_u.append(temp_us[W_limit:E_limit])
            sub_v.append(temp_vs[W_limit:E_limit])
            sub_topo.append(temp_topos[W_limit:E_limit])
        
        
        ## detailed plot
        fig1, ax1 = plt.subplots(1,1, figsize = (17,15))
        n = 3 # <11 and the plotting messes up, I don't know why
        #print(cube_xwind[t][p])
        u_n = sub_u[::n]#; u_n = u_n[0:len(y_coor) - 1:n]
        v_n = sub_v[::n]#; v_n = v_n[0:len(y_coor) - 1:n]
        
        u_n2 = []
        v_n2 = []
        for j in range(len(u_n)):
            u_n2.append(u_n[j][0::n])
            v_n2.append(v_n[j][0::n])
        
        q1 = ax1.contourf(dist_lon, dist_lat, sub_wind)#, cmap = plt.cm.YlGn)
        q2 = ax1.quiver(dist_lon[::n], dist_lat[::n], u_n2, v_n2, pivot = 'middle') 
        q3 = ax1.contour(dist_lon, dist_lat, sub_topo, cmap = plt.cm.gray)
        
        cbar1 = fig1.colorbar(q1, ax = ax1)
        cbar1.ax.set(xlabel = r'U [ms$^{-1}$]')
        ax1.set(title = f'Wind speed and direction at {p_levels[p]}hPa, {t + 6}hrs', 
                                                       ylabel = 'km', xlabel = 'km')
        fig1.savefig(f'../Figures/PNG/P_levels/{p_levels[p]}hPa/Wind/RA1M_4p4km_301_windmagn_subdomain_time{6 + t}_plevel{p_levels[p]}.png')
        fig1.savefig(f'../Figures/PDF/P_levels/{p_levels[p]}hPa/Wind/RA1M_4p4km_301_windmagn_subdomain_time{6 + t}_plevel{p_levels[p]}.pdf')

