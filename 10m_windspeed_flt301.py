# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 09:51:09 2020

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
from scipy.interpolate import interp2d

plt.ioff()
#%% global definitions
Case = '4p4km'


#%% first load and extract data
# main data cube 
xcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/windvector_10m/{Case}_10m_x_wind_24hrs_vera_301.nc', 'x_wind')
xdata = xcube.data
ycube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/windvector_10m/{Case}_10m_y_wind_24hrs_vera_301.nc', 'y_wind')
ydata = ycube.data

print(xcube)
#%%
magdata = np.sqrt(xdata**2 + ydata**2)
magdata[np.where(magdata == 0)] = np.nan


# land mask data cube
lmcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/land_binary_mask/{Case}_land_binary_mask_flt301.nc')
lmdata = lmcube.data
lmlon = lmcube.coord('grid_longitude').points
lmlat = lmcube.coord('grid_latitude').points
# mean-sealevel-pressure cube
mspcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/air_pressure_at_sea_level/{Case}_air_pressure_at_sea_level_24hrs_verb_301.nc')
mspdata = mspcube.data

msplon = mspcube.coord('grid_longitude').points
msplat = mspcube.coord('grid_latitude').points

gphcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/geopotential_height/{Case}_geopotential_height_24hrs_verd_301.nc')
gphdata = gphcube.data
gphlon = gphcube.coord('grid_longitude').points
gphlat = gphcube.coord('grid_latitude').points
gphcube = None
model_lat = xcube.coord('grid_latitude').points
model_lon = xcube.coord('grid_longitude').points
polelat = xcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
polelon = xcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
## free up memory/ unpoint cubes
#thcube = None, mspcube = None
#print(xcube)
pressure_labels = ['100hPa','150hPa', '200hPa','25hPa','300hPa','400hPa', '500hPa','600hPa', '650hPa','700hPa','750hPa','800hPa','850hPa','925hPa','950hPa','1000hPa']

#%%
## load quality controlled and interval-meaned observational data from Database
# using pandas
database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_sr_301.txt', delimiter = ' ')
db_u10m = database['u10m']
db_lon = database['lon']
db_lat = database['lat']
#%%
lons, lats = foundry.modf.unrotate_coords(model_lon, model_lat, polelon, polelat)
#%%
#foundry.toolkit.interp_obs_coords(lons,lats, model_lon, model_lat, db_lon,db_lat)
## meshgrid model x-y arrays
model_xx, model_yy = np.meshgrid(model_lon,model_lat)
## build 2d array as if lat-lon grid were flat, not curved
flat_lons = np.tile(lons[0], (len(lons[:,0]),1)) # taking edge of curved lons grid, and copy over all points in y
## lats need to be reshaped to create an array of distinct 1-element arrays before expanding
lats_border = lats[:,0]  #take edge of curved lats grid
lats_border = np.reshape(lats_border, (len(lats_border),1))
flat_lats = np.tile(lats_border, (len(lats[0])))
# make lwarg dict
keywords = {}
## before interpolation, which reorders the arrays, 
## we need to save the original order of elements for later using numpy argsort()
x_idx = np.argsort(np.argsort(db_lon)) # return list of indices from before ordering
y_idx = np.argsort(np.argsort(db_lat))

## interpolate using lats and lons edges (i.e. interpolating onto a flat lat-lon grid)
fx_flat = interp2d(lons[0], lats[:,0], model_xx)

fy_flat = interp2d(lons[0], lats[:,0], model_yy)
print(db_lon,db_lat)
# after interpolation, also resort arrays back into original order
dbx_flat = fx_flat(db_lon,db_lat, **keywords) # input database observation coordinates
dby_flat = fy_flat(db_lon,db_lat, **keywords) # this produces two 2d arrays from which the first index is what I want
dbx_flat = dbx_flat[:,x_idx]
dby_flat = dby_flat[y_idx,:]
#plt.plot(model_lon)
plt.plot(dbx_flat[0])
plt.show()
#plt.plot(model_lat)
plt.plot(dby_flat[:,0])
plt.show()
plt.plot(db_lon)
plt.show()
plt.plot(db_lat)
plt.show()

## interpolate onto deviation grid using newly found model x-y coordinates
lats_dev = lats - flat_lats # take grid deviation in latitude between flat and curved grid
lons_dev = lons - flat_lons # take grid deviation in longitude between flat and curved grid
fx_dev = interp2d(model_lon,model_lat,lons_dev) # now we're flipping the mapping 
fy_dev = interp2d(model_lon,model_lat,lats_dev)
# after interpolation, also resort arrays back into original order
dblon_dev = fx_dev(dbx_flat[0],dby_flat[:,0], **keywords) # and input the interpolated coordinates
dblat_dev = fy_dev(dbx_flat[0],dby_flat[:,0], **keywords)
dblon_dev = dblon_dev[:,x_idx]
dblat_dev = dblat_dev[y_idx,:]
#print(dblon_dev[:,0])
## now we use this to modify the measured observations coordinates
## and repeat the above process for one round of interpolations
db_lon_crd = db_lon - dblon_dev[0] # calculate 'corrected' coordinate for database obs
db_lat_crd = db_lat - dblat_dev[:,0]
# after interpolation, also resort arrays back into original order
dbx = fx_flat(db_lon_crd,db_lat_crd, **keywords)
dby = fy_flat(db_lon_crd,db_lat_crd, **keywords)
dbx = dbx[:,x_idx]
dby = dby[y_idx,:]
print('mean deviation in x =',np.mean(-dbx_flat[0] + dbx[0]))
print('mean deviation in y =',np.mean(-dby_flat[0] + dby[0]))
plt.scatter(db_lon,db_lat)
plt.show()
q = plt.contour(model_lon,model_lat,lons)
p = plt.contour(model_lon,model_lat,lats)
plt.clabel(q)
plt.clabel(p)
plt.scatter(dbx_flat[0],dby_flat[:,0])
plt.show()
plt.contour(model_lon,model_lat,lons,colors = 'k')
plt.contour(model_lon,model_lat,lats,colors = 'k')
#plt.scatter(dbx[0],dby[:,0])
plt.show()
fig,axes = plt.subplots(2,2, figsize = (18,16))
axes[0,0].contour(model_lon,model_lat,flat_lons,colors = 'k')
axes[0,0].contour(model_lon,model_lat,flat_lats,colors = 'k')
axes[0,1].contour(model_lon,model_lat,lons,colors = 'k')
axes[0,1].contour(model_lon,model_lat,lats,colors = 'k')
q2a = axes[1,0].contour(model_lon,model_lat,lons_dev)
q2b = axes[1,1].contour(model_lon,model_lat,lats_dev)
plt.clabel(q2a);plt.clabel(q2b)
#plt.savefig('grid_comparisons_with_deviations.png')
plt.show()



#%%

#fy = interpolate.interp2d(lons,lats, model_yy)
#fy = interpolate.interp2d(lons,lats, model_xx)
#db_y = fy(db_lon,db_lat)

#%%

domain = 'fulldomain'
if domain == 'fulldomain':
    lons, lats = foundry.modf.unrotate_coords(model_lon, model_lat, polelon, polelat)
    msp_levels = np.arange(600, 1200, 1)
    gph_levels = np.arange(0,12000,10)
    
    counter = 0
    for level in range(1):
        for time in range(24):
           
            if (time) % 1 == 0:
                mspdata_at_t = mspdata[time]
                #mspdata_at_t[np.where(lmdata == 1)] = np.nan
                counter += 1; print('counter:',counter, 'time:', time, 'level:', level)
                fig0,ax0,norm0 = foundry.plot_hoz_xsec_mass(magdata, lmdata, model_lon, model_lat, lons, lats, lmlon, lmlat, time, level,
                                                                 data_label = r'Windspeed, $ms^{-1}$', domain = domain, 
                                                                 flight = 301, savefig = True, res = Case, unrotated = True,
                                                                 variable_name = 'Surface horizontal windspeed', time_norm = True,
                                                                 variable_in_file_name = '10m_windvector', 
                                                                 figure_path_pdf = f'../../Figures/PDF/301/u-bu807/P_levels/windvector_pl/{Case}/{domain}/verstash/surface_10m',
                                                                 figure_path_png = f'../../Figures/PNG/301/u-bu807/P_levels/windvector_pl/{Case}/{domain}/verstash/surface_10m',
                                                                 contour_type = 'pcolor', colourmap = 'cool',sampling_rate = 'hourly',
                                                                 msp_contour = True, msp_data = mspdata_at_t/100, msp_lon = msplon, msp_lat = msplat,
                                                                 surface_variable = True, msp_levels = msp_levels, pressure_levels = True,
                                                                 gph_contour = False, gph_data = gphdata, gph_lon = gphlon, gph_lat = gphlat, gph_levels = gph_levels,
                                                                 quiver = True, quiver_u = xdata, quiver_v = ydata,
                                                                 quiver_n = 15, gph_style = 'solid',figsize = (18,15),
                                                                 cscale_override = False, cscale_min = 0, cscale_max = 16,
                                                                 obvs_on = False, obvs_z = db_u10m, obvs_x = dbx[0], obvs_y = dby[:,0], PDF = False, verstash = True)
                
                                                                 

    
                plt.show()