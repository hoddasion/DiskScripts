"""
Created on Sun Nov 10 14:39:23 2019

@author: Wilhelm Hodder

Script to create an orographic map of the 1.5km domain
"""
#%% module imports
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import iris
import iris.coord_categorisation
import iris.analysis.cartography as cairis
#import CubeCrossSectioner_UK as ccs

import model_functions as modf
import itertools
import datetime
import xarray
import pandas as pd
#%% global configurations 

import warnings
warnings.filterwarnings("ignore") ##BAD! If something doesn't work, re-enable warnings to check
font = {'size' : 18}
matplotlib.rc('font', **font)



#%% load land data
land_file_name = '.nc'
land_path = f'../Model_Data/u-bk574/nc/{land_file_name}'

modf.file_info(land_path)

stash = 'surface_altitude'
orography_set, land_latitude, land_longitude = modf.return_cube_components(land_path, stash, 
                                                        lat_name = 'grid_latitude',
                                                        lon_name = 'grid_longitude')
land_mask, land_latitude, land_longitude = modf.return_cube_components(land_path, 'land_binary_mask', 
                                                        lat_name = 'grid_latitude',
                                                        lon_name = 'grid_longitude')
orography = orography_set[0][0]
coastline = land_mask[0][0]

orocube = iris.load_cube(land_path, stash)

#%% 

lons, lats = modf.unrotate_coords(land_longitude, land_latitude,
                                  NPoleLon = orocube.coord('grid_longitude').coord_system.grid_north_pole_longitude,
                                  NPoleLat = orocube.coord('grid_latitude').coord_system.grid_north_pole_latitude)

#%% create dataset for region boxes

South = 60; North = 250; West = 130; East = 230
box_reg1 = np.array([]) 
BOT1 = np.linspace(West, East)
#BOT1_data = np.
RHS1 = np.linspace(South, North)

#%% replace ocean points with nan

orography_set[orography_set <= 0] = np.nan

#%% define dictionary of keyword arguments for contour plots

keywords = {'contour_number' : 10, 'contour_type' : 'contour'}


#%% plot full domain map
plt.close()

fig0 = modf.make_horizontal_model_contourf_linear(orography_set, coastline, land_longitude, land_latitude, 
                                                  lons, lats, land_longitude, land_latitude, time = 0, level = 0, 
                                                  data_label = 'Elevation above sealevel, [m]', domain = 'fulldomain', flight =301, 
                                                  variable_in_file_name = 'surface_altitude', savefig = False, colourmap = 'summer',
                                                  title = 'Orography over 1p5km full domain',
                                                  figure_path_pdf = '../Figures/PDF/Orography',
                                                  figure_path_png = '../Figures/PNG/Orography',
                                                  contour_number = 10,
                                                  contour_type = 'contour', vmin = -1000)
plt.close()

#%% subset 1

## boundaries (in points):
South = 60; North = 250; West = 130; East = 230

# subset 4d data
sub_data_region1 = modf.subset_4d(orography_set, South, North, West, East)

       
# two-dimensional subsetting of land-mask
sub_coast_region1 = modf.subset_2d(coastline, South, North, West, East )

# subset coordinates and unrotate
sub_lon1 = modf.subset_1d(land_longitude, West, East)
sub_lat1 = modf.subset_1d(land_latitude, South, North)

sub_lons1, sub_lats1 = modf.unrotate_coords(sub_lon1,
                                          sub_lat1,
                                          orocube.coord('grid_longitude').coord_system.grid_north_pole_longitude,
                                          NPoleLat = orocube.coord('grid_latitude').coord_system.grid_north_pole_latitude)

#%%

fig0 = modf.make_horizontal_model_contourf_linear(sub_data_region1, sub_coast_region1, sub_lon1, sub_lat1, 
                                                  sub_lons1, sub_lats1, sub_lon1, sub_lat1, time = 0, level = 0, 
                                                  data_label = 'Elevation above sealevel, [m]', domain = 'subset1', flight =301, 
                                                  variable_in_file_name = 'surface_altitude', savefig = False, colourmap = 'summer',
                                                  title = 'Orography over 1p5km region 1',
                                                  figure_path_pdf = '../Figures/PDF/Orography',
                                                  figure_path_png = '../Figures/PNG/Orography',
                                                  contour_number = 10,
                                                  contour_type = 'contour')
plt.close()

#%% subset 2

## boundaries (in points):
South = 60; North = 150; West = 130; East = 210

# subset 4d data
sub_data_region1 = modf.subset_4d(orography_set, South, North, West, East)

       
# two-dimensional subsetting of land-mask
sub_coast_region1 = modf.subset_2d(coastline, South, North, West, East )

# subset coordinates and unrotate
sub_lon1 = modf.subset_1d(land_longitude, West, East)
sub_lat1 = modf.subset_1d(land_latitude, South, North)

sub_lons1, sub_lats1 = modf.unrotate_coords(sub_lon1,
                                          sub_lat1,
                                          orocube.coord('grid_longitude').coord_system.grid_north_pole_longitude,
                                          NPoleLat = orocube.coord('grid_latitude').coord_system.grid_north_pole_latitude)

#%%

fig, ax = modf.make_horizontal_model_contourf_linear(sub_data_region1, sub_coast_region1, sub_lon1, sub_lat1, 
                                                  sub_lons1, sub_lats1, sub_lon1, sub_lat1, time = 0, level = 0, 
                                                  data_label = 'Elevation above sealevel, [m]', domain = 'subset2', flight =301, 
                                                  variable_in_file_name = 'surface_altitude', savefig = False, colourmap = 'summer',
                                                  title = 'Orography over 1p5km region 2',
                                                  figure_path_pdf = '../Figures/PDF/Orography',
                                                  figure_path_png = '../Figures/PNG/Orography',
                                                  contour_number = 10,
                                                  contour_type = 'contour')


#%% load observational data
obvs_vars = ['gps_alt_50hz','lon_50hz','lat_50hz']
obvs_file = xarray.open_dataset('../Obvs_Data/UEA_qc_MASIN_data/MASIN_flightdata_301.nc')

#%% extract observational data file contents

obvs_alt = np.array(obvs_file['gps_alt_50hz']) # altitude for orographic map plot for comparison purposes
obvs_lon = np.array(obvs_file['lon_50hz'])
obvs_lat = np.array(obvs_file['lat_50hz'])

## rotate observational data coordinates

obvs_rot_lon, obvs_rot_lat = cairis.rotate_pole(obvs_lon, obvs_lat,
                                  orocube.coord('grid_longitude').coord_system.grid_north_pole_longitude,
                                  orocube.coord('grid_latitude').coord_system.grid_north_pole_latitude)

#%% region 3

## boundaries (in points):
South = 60; North = 200; West = 130; East = 210

# subset 4d data
sub_data_region1 = modf.subset_4d(orography_set, South, North, West, East)

       
# two-dimensional subsetting of land-mask
sub_coast_region1 = modf.subset_2d(coastline, South, North, West, East )

# subset coordinates and unrotate
sub_lon1 = modf.subset_1d(land_longitude, West, East)
sub_lat1 = modf.subset_1d(land_latitude, South, North)

sub_lons1, sub_lats1 = modf.unrotate_coords(sub_lon1,
                                          sub_lat1,
                                          orocube.coord('grid_longitude').coord_system.grid_north_pole_longitude,
                                          NPoleLat = orocube.coord('grid_latitude').coord_system.grid_north_pole_latitude)




#%% load 60s datset from database

# using pandas
database = pd.read_csv('../Obvs_Data/Databases/IGP_flights_database_obs_60s_301.txt', delimiter = ' ')
print(database[['altgps', 'lon', 'lat']])

db_alt = np.array(database['altgps'])
db_lon = np.array(database['lon']); db_lat = np.array(database['lat'])

db_lon_rot, db_lat_rot = cairis.rotate_pole(db_lon, db_lat,
                                  orocube.coord('grid_longitude').coord_system.grid_north_pole_longitude,
                                  orocube.coord('grid_latitude').coord_system.grid_north_pole_latitude)
#%% plot database gps altitude over orography map
fig, ax = modf.make_horizontal_model_contourf_linear(sub_data_region1, sub_coast_region1, sub_lon1, sub_lat1, 
                                                  sub_lons1, sub_lats1, sub_lon1, sub_lat1, time = 0, level = 0, 
                                                  data_label = 'Elevation above sealevel, [m]', domain = 'subset2', flight =301, 
                                                  variable_in_file_name = 'surface_with_flight_altitude', savefig = True, colourmap = 'gist_earth',
                                                  title = 'Orography over 1p5km with flight 301 gps altitude',
                                                  figure_path_pdf = '../Figures/PDF/Orography',
                                                  figure_path_png = '../Figures/PNG/Orography',
                                                  contour_number = 30,
                                                  contour_type = 'pcolor',
                                                  obvs_on = True, obvs_z = db_alt, obvs_size = 100.0,
                                                  obvs_x = db_lon_rot + 360, obvs_y = db_lat_rot,
                                                  figsize = (15,20), coast = False, vmin = -1000)
plt.show()