# -*- coding: utf-8 -*-
"""
Created on Sun Nov 10 14:39:23 2019

@author: Wilhelm Hodder

Script to create various horizontal cross-sections from model data on model levels
"""
#%% module imports
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
import model_functions as modf
import itertools
import datetime

#%% global configurations 

import warnings
warnings.filterwarnings("ignore") ##BAD! If something doesn't work, re-enable warnings to check
font = {'size' : 18}
matplotlib.rc('font', **font)

figure_path_pdf = '../Figures/PDF/301/Model_levels/Vertical_velocity/'
figure_path_png = '../Figures/PNG/301/Model_levels/Vertical_velocity/'

# Turn interactive plotting off
plt.ioff()
#%% load land data
land_file_name = '20180312T0000Z_IcelandGreenlandSeas_1p5km_RA1M_pf000.nc'
land_path = f'../Model_Data/u-bk574/nc/{land_file_name}'

modf.file_info(land_path)

stash = 'land_binary_mask'
land_mask, land_latitude, land_longitude = modf.return_cube_components(land_path, stash, 
                                                        lat_name = 'grid_latitude',
                                                        lon_name = 'grid_longitude')

#%% load first variable to analyse: Theta on theta levels
file_name = '20180312T0000Z_IcelandGreenlandSeas_1p5km_RA1M_pi000.nc'
path1 = f'../Model_Data/u-bk574/nc/'
stash1 = 'air_potential_temperature'
modf.file_info(f'{path1}{file_name}')#, stash1)
modf.cube_info(f'{path1}{file_name}', stash1)
modf.level_info(f'{path1}{file_name}', stash1, 3)
model_height = modf.get_model_levels(f'{path1}{file_name}', stash1)
print(model_height)
arguments = {'Hybrid height' : True}
theta, theta_lat, theta_lon = modf.return_concatenated_cube_components(path1, 'i', '1p5km',
                                                          stash1,
                                                          file_number =4,
                                                          lat_name = 'grid_latitude',
                                                          lon_name = 'grid_longitude')

#%% load second variable to analyse: 
file_name = '20180312T0000Z_IcelandGreenlandSeas_1p5km_RA1M_ph000.nc'

stash2 = 'upward_air_velocity'
modf.file_info(f'{path1}{file_name}')#, stash1)
modf.cube_info(f'{path1}{file_name}', stash2)
modf.level_info(f'{path1}{file_name}', stash2, 3)
#model_height = modf.get_model_levels(f'{path1}{file_name}', stash2)
#print(model_height)
arguments = {'Hybrid height' : True}
w, w_lat, w_lon = modf.return_concatenated_cube_components(path1, 'h', '1p5km',
                                                          stash2,
                                                          file_number = 4,
                                                          lat_name = 'grid_latitude',
                                                          lon_name = 'grid_longitude')

#%%
# load the cube itself for good measure 
paths = [f'{path1}{file_name}', land_path]
w_cube = iris.load_cube(paths, stash2)

#print('\n\nw_cube with orography =\n',w_cube, '\n\n')


#%% unrotate coordinate fields

rotated_2D_lon, rotated_2D_lat = np.meshgrid(land_longitude, land_latitude)
NPoleLon = w_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude
NPoleLat = w_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude
print('Pole lon, lat =', NPoleLon, NPoleLat)
import iris.analysis.cartography as icarto
lons, lats = icarto.unrotate_pole(rotated_2D_lon, rotated_2D_lat,
                                NPoleLon, NPoleLat)
print(rotated_2D_lat)
print(NPoleLon, NPoleLat)

#lon = lons[0]
#lat = lats[:,0]

#%% create figure for vertical velocity fields




for i in range(len(w)): # temporal range
    for j in range(0):  # vertical range


        ## get min and max of vertical velocity data at eacht time and level
        day_wmin, day_wmax = modf.data_extrema(w[:][j]) # over day
        
        
        # make list of contour levels
        
        w_min = np.nanmin(w[i][j])
        levels = np.linspace(day_wmin - 0.1, day_wmax + 0.1)# 10)
        # formatting time string/datetime object
        hour = i//2
        minute = (i/2 - hour)*60
        print(hour, minute); minute = int(minute); hour = int(hour)
        TOD = datetime.time(hour, minute)  # time of day
        # 
        fig, ax = plt.subplots(1,1, figsize = (15,15))
        q1 = ax.contourf(w_lon, w_lat, w[i][j], levels = levels,
                         vmin = day_wmin, vmax = day_wmax, cmap = plt.cm.seismic)
        # add coastlines from model data
        q2 = ax.contour(w_lon, w_lat, land_mask[0][0])
        # add real latitude and longitude lines
        q4 = ax.contour(w_lon, w_lat, lons, colors = 'k')
        q5 = ax.contour(w_lon, w_lat, lats, colors = 'k')
        # label latitude and longitude lines
        label_keys = {'fontsize' : 9, 'inline' : 1}
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid' # make contours for negative values solid, not dashed
        q6 = ax.clabel(q4, **label_keys)
        q7 = ax.clabel(q5, **label_keys)
        # hide axis labels
        q8 = ax.set_yticklabels([]); q9 = ax.set_xticklabels([])
        # titles and text labels
        q3 = ax.set(Title = f'Vertical velocity at model level {j + 1} at {TOD}UTC',
                    ylabel = 'Latitude',
                    xlabel = 'Longitude')
        # colourbar
        cbar0 = fig.colorbar(q1, ax = ax)
        cbar0.ax.set(ylabel = r'w [ms$^{-1}$]')
        # save figures
        figure_name = f'1p5km_fulldomain_{stash2}_20180312-301_TH{hour}M{minute}_ML{j+1}'
        #fig.savefig(f'{figure_path_pdf}fulldomain/{j + 1}/{figure_name}.pdf')
        #fig.savefig(f'{figure_path_png}fulldomain/{j + 1}/{figure_name}.png')

#%% region 1 subsetting
## boundaries (in points):
South = 60; North = 250; West = 130; East = 230
sub_rot_lon = land_longitude[West:East] # subset rotated longitude
sub_rot_lat = land_latitude[South:North] # subset rotated latitude
sub_rot_lon_2d, sub_rot_lat_2d = np.meshgrid(sub_rot_lon, sub_rot_lat)
sub_lons, sub_lats =  icarto.unrotate_pole(sub_rot_lon_2d, sub_rot_lat_2d,
                                NPoleLon, NPoleLat)
i = 0; j = 5
for i in range(len(w)): # temporal range
    for j in range(0):  # vertical range
        ## in-loop subsetting
        full_domain_at_level = w[i][j]
        print(len(full_domain_at_level))
        sub_data = []
        # two-dimensional subsetting of data
        sub_data_rows = full_domain_at_level[South:North]
        for row in range(len(sub_data_rows)):
            sub_data.append(sub_data_rows[row][West:East])
        # two-dimensional subsetting of land-mask
        sub_land = []
        sub_land_rows = land_mask[0][0][South:North]
        for row in range(len(sub_land_rows)):
            sub_land.append(sub_land_rows[row][West:East]) 
            
            
        ## get min and max of vertical velocity data at eacht time and level
        day_wmin, day_wmax = modf.data_extrema(sub_data) # over day
        print(len(sub_rot_lon), len(sub_rot_lat))
        
        # make list of contour levels
        
        w_min = np.nanmin(sub_data)
        levels = np.linspace(day_wmin - 0.1, day_wmax + 0.1)#, 50)
        # formatting time string/datetime object
        hour = i//2
        minute = (i/2 - hour)*60
        print(hour, minute)
        TOD = datetime.time(int(hour), int(minute))  # time of day
        # 
        fig, ax = plt.subplots(1,1, figsize = (10,15))
        q1 = ax.contourf(sub_rot_lon, sub_rot_lat, sub_data, levels = levels,
                         vmin = day_wmin, vmax = day_wmax, cmap = plt.cm.seismic)
        # add coastlines from model data
        q2 = ax.contour(sub_rot_lon, sub_rot_lat, sub_land)
        # add real latitude and longitude lines
        print(len(sub_rot_lon), len(sub_lons))
        q4 = ax.contour(sub_rot_lon, sub_rot_lat, sub_lons, colors = 'k')
        q5 = ax.contour(sub_rot_lon, sub_rot_lat, sub_lats, colors = 'k')
        # label latitude and longitude lines
        label_keys = {'fontsize' : 9, 'inline' : 1}
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid' # make contours for negative values solid, not dashed
        q6 = ax.clabel(q4, **label_keys)
        q7 = ax.clabel(q5, **label_keys)
        # hide axis labels
        q8 = ax.set_yticklabels([]); q9 = ax.set_xticklabels([])
        # titles and text labels
        q3 = ax.set(Title = f'Vertical velocity at model level {j + 1} at {TOD}UTC',
                    ylabel = 'Latitude',
                    xlabel = 'Longitude')
        # colourbar
        cbar0 = fig.colorbar(q1, ax = ax)
        cbar0.ax.set(ylabel = r'w [ms$^{-1}$]')
        # save figures
        figure_name = f'1p5km_subset1_{stash2}_20180312-301_TH{hour}M{minute}_ML{j+1}'
        #fig.savefig(f'{figure_path_pdf}subset1/{j + 1}/{figure_name}.pdf')
        #fig.savefig(f'{figure_path_png}subset1/{j + 1}/{figure_name}.png')
#%% region 2 subsetting
## boundaries (in points):
South = 60; North = 150; West = 130; East = 210
sub_rot_lon = land_longitude[West:East] # subset rotated longitude
sub_rot_lat = land_latitude[South:North] # subset rotated latitude
sub_rot_lon_2d, sub_rot_lat_2d = np.meshgrid(sub_rot_lon, sub_rot_lat)
sub_lons, sub_lats =  icarto.unrotate_pole(sub_rot_lon_2d, sub_rot_lat_2d,
                                NPoleLon, NPoleLat)

## animation prep
images = []
for i in range(len(w)): # temporal range
    for j in range(0):  # vertical range
        ## in-loop subsetting
        full_domain_at_level = w[i][j]
        print(len(full_domain_at_level))
        sub_data = []
        # two-dimensional subsetting of data
        sub_data_rows = full_domain_at_level[South:North]
        for row in range(len(sub_data_rows)):
            sub_data.append(sub_data_rows[row][West:East])
        # two-dimensional subsetting of land-mask
        sub_land = []
        sub_land_rows = land_mask[0][0][South:North]
        for row in range(len(sub_land_rows)):
            sub_land.append(sub_land_rows[row][West:East]) 
            
            
        ## get min and max of vertical velocity data at eacht time and level
        day_wmin, day_wmax = modf.data_extrema(sub_data) # over day
        print(len(sub_rot_lon), len(sub_rot_lat))
        
        # make list of contour levels
        
        w_min = np.nanmin(sub_data)
        levels = np.linspace(day_wmin - 0.1, day_wmax + 0.1)#, 50)
        # formatting time string/datetime object
        hour = i//2
        minute = (i/2 - hour)*60
        print(hour, minute)
        TOD = datetime.time(int(hour), int(minute))  # time of day
        # 
        fig, ax = plt.subplots(1,1, figsize = (15,15))
        q1 = ax.contourf(sub_rot_lon, sub_rot_lat, sub_data, levels = levels,
                         vmin = day_wmin, vmax = day_wmax, cmap = plt.cm.seismic)
        # add coastlines from model data
        q2 = ax.contour(sub_rot_lon, sub_rot_lat, sub_land)
        # add real latitude and longitude lines
        print(len(sub_rot_lon), len(sub_lons))
        q4 = ax.contour(sub_rot_lon, sub_rot_lat, sub_lons, colors = 'k')
        q5 = ax.contour(sub_rot_lon, sub_rot_lat, sub_lats, colors = 'k')
        # label latitude and longitude lines
        label_keys = {'fontsize' : 9, 'inline' : 1}
        matplotlib.rcParams['contour.negative_linestyle'] = 'solid' # make contours for negative values solid, not dashed
        q6 = ax.clabel(q4, **label_keys)
        q7 = ax.clabel(q5, **label_keys)
        # hide axis labels
        q8 = ax.set_yticklabels([]); q9 = ax.set_xticklabels([])
        # titles and text labels
        q3 = ax.set(Title = f'Vertical velocity at model level {j + 1} at {TOD}UTC',
                    ylabel = 'Latitude',
                    xlabel = 'Longitude')
        # colourbar
        cbar0 = fig.colorbar(q1, ax = ax)
        cbar0.ax.set(ylabel = r'w [ms$^{-1}$]')
        # append figure to images list for animation purposes
        images.append(fig)
        # save figures
        figure_name = f'1p5km_subset2_{stash2}_20180312-301_TH{hour}M{minute}_ML{j+1}'
        #fig.savefig(f'{figure_path_pdf}subset2/{j + 1}/{figure_name}.pdf')
        #fig.savefig(f'{figure_path_png}subset2/{j + 1}/{figure_name}.png')

#%% test compact functions
South = 60; North = 150; West = 130; East = 210
sub_data_region2 = modf.subset_4d(w, South, North, West, East)
       
# two-dimensional subsetting of land-mask
sub_land_region2 = modf.subset_2d(land_mask[0][0], South, North, West, East )

# subset coordinates and unrotate
sub_lon = modf.subset_1d(w_lon, West, East)
sub_lat = modf.subset_1d(w_lat, South, North)
sub_lons, sub_lats = modf.unrotate_coords(sub_lon,
                                          sub_lat,
                                          w_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude,
                                          NPoleLat = w_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude)
#%% test compact functions
import matplotlib.animation as animation
anifig = plt.gcf() # initialise animation figure
figure_name = f'1p5km_subset2_{stash2}_20180312-301_TH{hour}M{minute}_ML{j+1}'   
print(len(sub_land_region2[0])); print(len(land_mask[0][0][0]))
from model_functions import make_horizontal_model_contourf 
counter = 0


for level in range(10):
    anifig = plt.figure() # animated figure initialisation
    images = []
    for time in range(48):
        counter += 1
        print(counter)
        fig = make_horizontal_model_contourf(sub_data_region2, sub_land_region2, sub_lon, sub_lat, sub_lons, sub_lats,
                                       time = time, level = level,
                                       data_label =  r'w [ms$^{-1}$]',
                                       figure_path_pdf = f'{figure_path_pdf}subset2',
                                       figure_path_png = f'{figure_path_png}subset2',
                                       file_name = figure_name,
                                       variable_name = 'Vertical velocity')
        images.append([fig])
    animated_day = animation.ArtistAnimation(anifig, images, interval = 48, blit = True, repeat_delay = 200)
    print(len(images[0]))
    figure_path_gif = '../Figures/GIF/301/Model_levels/Vertical_velocity/'
    animated_day.save(f'{figure_path_gif}{level + 1}.gif')
