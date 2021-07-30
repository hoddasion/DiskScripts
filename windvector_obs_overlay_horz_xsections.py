# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 19:49:27 2020

@author: kse18nru
"""


#%% Import modules and scripts
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import iris
from iris.analysis.cartography import rotate_pole
import iris.coord_categorisation
import iris.plot as iplt
import model_foundry as foundry
import pandas as pd

#%% 


#%% globally set variables

Case = '1p5km'
domain = 'snaesfellness'
pressure_labels = ['100hPa','150hPa', '200hPa','250hPa','300hPa','400hPa', '500hPa','600hPa', '650hPa','700hPa','750hPa','800hPa','850hPa','925hPa','950hPa','1000hPa']

#%% load data cubes

## x wind (u)
xcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/windvector_pl/{Case}_x_wind_24hrs_verc_301.nc', 'x_wind')
xdata = xcube.data
ulat = xcube.coord('grid_latitude').points
ulon = xcube.coord('grid_longitude').points
polelat = xcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
polelon = xcube.coord('grid_longitude').coord_system.grid_north_pole_longitude

## y wind (v)
ycube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/windvector_pl/{Case}_y_wind_24hrs_verc_301.nc', 'y_wind')
ydata = ycube.data

# ## magnitude wind
# magcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/windvector_10m/{Case}_mag_wind_24hrs_vera_301.nc', '10m_mag_wind')
# magdata = magcube.data
# maglon = magcube.coord('grid_longitude').points
# maglat = magcube.coord('grid_latitude').points

## land binary mask (coastline)
lmcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/land_binary_mask/{Case}_land_binary_mask_flt301.nc')
lmdata = lmcube.data
lmlon = lmcube.coord('grid_longitude').points
lmlat = lmcube.coord('grid_latitude').points

## geopotential height
gphcube = iris.load_cube(f'../../Model_Data/u-bu807/nc/Control/geopotential_height/{Case}_geopotential_height_24hrs_verd_301.nc')
gphdata = gphcube.data
gphlon = gphcube.coord('grid_longitude').points
gphlat = gphcube.coord('grid_latitude').points

#%%

print(xcube.coord('time'))

#%% subset le data

if domain == 'snaesfellness':
    ## boundaries (in points):
    if Case == '4p4km':
        South = 143; North = 155; West = 150; East = 178
        quiver_n = 1
    if Case == '1p5km':
        South = 150; North = 220; West = 90; East = 160
        quiver_n = 2
    if Case == '0p5km':
        South = 280; North = 370; West = 50; East = 250
        quiver_n = 6

if domain == 'peninsulas':
    ## boundaries (in points):
    if Case == '4p4km':
        South = 115; North = 200; West = 150; East = 185
        quiver_n = 2
    if Case == '1p5km':
        South = 110; North = 350; West = 90; East = 185
        quiver_n = 5
    if Case == '0p5km':
        South = 80; North = -1; West = 50; East = 350
        quiver_n = 15

print('xdata', np.shape(xdata))
xwind = foundry.modf.subset_4d(xdata, South, North, West, East)
ywind = foundry.modf.subset_4d(ydata, South, North, West, East)

magwind = (xwind**2 + ywind**2)**0.5
gph = foundry.modf.subset_4d(gphdata, South, North, West, East) 
# two-dimensional subsetting of land-mask
lm = foundry.modf.subset_2d(lmdata, South, North, West, East )

ulon = foundry.modf.subset_1d(ulon, West, East)
ulat = foundry.modf.subset_1d(ulat, South, North)
lmlon = foundry.modf.subset_1d(lmlon, West, East)
lmlat = foundry.modf.subset_1d(lmlat, South, North)
gphlon = foundry.modf.subset_1d(gphlon, West, East)
gphlat = foundry.modf.subset_1d(gphlat, South, North)

print('xwind', np.shape(xwind))
print('gph', np.shape(gph))

#%%
## load quality controlled and interval-meaned observational data from Database
# using pandas
database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_60s_301.txt', delimiter = ' ')
db_legno = database['legno']; legno = 8
legcond = db_legno == legno
db_u = np.array(database['u'][legcond])
db_v = np.array(database['v'][legcond])
db_wsp = np.array(database['wsp'][legcond])
db_lon = np.array(database['lon'][legcond])
db_lat = np.array(database['lat'][legcond])
db_alt = np.array(database['altgps'][legcond])

db_mtime = np.array(database['meantime'][legcond])

#%% rotate db coords

rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
rot_db_lon = rot_db_lon + 360

#%% set normalisation  for db and model data
level = 13
time = 5
vmax = 20
vmin = 0
hour = time * 3
#model _max = np.nanmax(magwind[:])

#%% ploet
quiver_kwargs = {}

cmap = 'cool'
fig, ax = plt.subplots(1,1, figsize = (15,15))
norm = plt.Normalize()
q0 = ax.pcolormesh(ulon, ulat, magwind[time, level], shading = 'auto', vmax = vmax, vmin = vmin, cmap = cmap)
n = 3; q01 = ax.quiver(ulon[::n],ulat[::n],xwind[time,level,::n,::n], ywind[time,level,::n,::n])
q1 = ax.contour(lmlon,lmlat, lm, levels = [0.5], colors = 'k', linewidths = 3)
q2 = ax.scatter(rot_db_lon, rot_db_lat, c = db_wsp, vmax = vmax, vmin = vmin, edgecolors = 'k', marker = 'o', s = 200, cmap = cmap )
n = 1; q02 = ax.quiver(rot_db_lon[::n], rot_db_lat[::n], db_u[::n],db_v[::n])
q3 = ax.contour(gphlon,gphlat,gph[time,level], colors = 'k', levels = np.arange(0,2000,5))  
plt.clabel(q3, inline=False, fontsize = 14)
cbar = plt.colorbar(q0); cbar.ax.set(ylabel = r'Horizontal windspeed [ms$^{-1}$]')
q4 = ax.set(title = f'{Case} Windspeed on {pressure_labels[level]} level at TH{time*3}')
q8 = ax.set_yticklabels([]); q9 = ax.set_xticklabels([])
plt.savefig(f'../../Figures/PNG/301/u-bu807/P_levels/windvector_pl/{Case}/{domain}/verstash/obs/{Case}_windvector_dir_{pressure_labels[level]}_TH{hour}_flt301leg{legno}_60s.png')
