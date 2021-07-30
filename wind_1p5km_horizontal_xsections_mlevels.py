# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 13:57:10 2019

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
import matplotlib.animation as animation

#%% global configurations 

import warnings
warnings.filterwarnings("ignore") ##BAD! If something doesn't work, re-enable warnings to check
font = {'size' : 18}
matplotlib.rc('font', **font)

figure_path_pdf = '../Figures/PDF/301/Model_levels/Horizontal_velocity/Horizontal/'
figure_path_png = '../Figures/PNG/301/Model_levels/Horizontal_velocity/Horizontal/'
figure_path_mp4 = '../Figures/MP4/301/Model_levels/Horizontal_velocity/Horizontal/'
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


#%% load both wind components
data_path = '../Model_data/u-bk574/nc/wind_components_ml/'
stash1 = 'm01s15i002'
u, u_lat, u_lon = modf.return_concatenated_cube_components(data_path, 'i', '1p5km',
                                                          stash1,
                                                          file_number =4,
                                                          lat_name = 'grid_latitude',
                                                          lon_name = 'grid_longitude')
stash2 = 'm01s15i003'
v, v_lat, v_lon = modf.return_concatenated_cube_components(data_path, 'i', '1p5km',
                                                          stash2,
                                                          file_number =4,
                                                          lat_name = 'grid_latitude',
                                                          lon_name = 'grid_longitude')

#%% combine windspeed components to make wind magnitude

magnitude = np.sqrt(u**2 + v**2)

#%% plot

fig0, ax0 = plt.subplots(1,1, figsize = (12,12))
q0 = ax0.contourf(u_lon, u_lat, magnitude[0][0])
# colourbar
cbar0 = fig0.colorbar(q0, ax = ax0)

#%% make lons and lats arrays
land_cube = iris.load_cube(f'{land_path}', stash)
lons, lats = modf.unrotate_coords(u_lon, u_lat, 
                                  land_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude,
                                  land_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude)

#%%
print(len(land_latitude), len(land_longitude), '\n',
      len(u_lat), len(u_lon))

#%%
counter = 0
full_domain_on = True
if full_domain_on == True:
    for level in range(40):
        for time in range(48):
            counter += 1; 
            print(counter, 'model-level:', level + 1, 'time-step:', time + 1)
            fig = modf.make_horizontal_model_contourf_linear(magnitude, land_mask[0][0], u_lon, u_lat,  lons, lats, land_longitude, land_latitude,
                                                      level = level, time = time, colourmap = 'viridis',
                                                      quiver = True, quiver_u = u, quiver_v = v, quiver_n = 11,
                                                      variable_in_file_name = 'horizontal_velocity',
                                                      figure_path_pdf = '../Figures/PDF/301/Model_levels/Horizontal_velocity',
                                                      figure_path_png = '../Figures/PNG/301/Model_levels/Horizontal_velocity',
                                                      data_label = r'ms$^{-1}$', domain = 'fulldomain', 
                                                      flight = 301, savefig = True, res = '1p5km',
                                                      variable_name = 'Horizontal velocity',
                                                      contour_type = 'contourf')
            plt.close()
            
    
    plt.close()        
#%% subset for region 1
## boundaries (in points):
South = 60; North = 250; West = 130; East = 230

# subset 4d data
sub_data_region1 = modf.subset_4d(magnitude, South, North, West, East)
sub_u_region1 = modf.subset_4d(u, South, North, West, East)
sub_v_region1 = modf.subset_4d(v, South, North, West, East)
       
# two-dimensional subsetting of land-mask
sub_land_region1 = modf.subset_2d(land_mask[0][0], South, North, West, East )

# subset coordinates and unrotate
sub_lon1 = modf.subset_1d(u_lon, West, East)
sub_lat1 = modf.subset_1d(u_lat, South, North)
sub_land_lon1 = modf.subset_1d(land_longitude, West, East)
sub_land_lat1 = modf.subset_1d(land_latitude, South, North)
sub_lons1, sub_lats1 = modf.unrotate_coords(sub_lon1,
                                          sub_lat1,
                                          land_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude,
                                          NPoleLat = land_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude)

#%%
counter = 0
region1_on = True
if region1_on == True:
    for level in range(40):
        for time in range(48):
            counter += 1; 
            print(counter, 'model-level:', level + 1, 'time-step:', time + 1)
            fig = modf.make_horizontal_model_contourf_linear(sub_data_region1, sub_land_region1, sub_lon1, sub_lat1,  sub_lons1, sub_lats1, 
                                                              sub_land_lon1, sub_land_lat1,
                                                      level = level, time = time, colourmap = 'viridis',
                                                      quiver = True, quiver_u = sub_u_region1, quiver_v = sub_v_region1, quiver_n = 5,
                                                      variable_in_file_name = 'horizontal_velocity',
                                                      figure_path_pdf = '../Figures/PDF/301/Model_levels/Horizontal_velocity',
                                                      figure_path_png = '../Figures/PNG/301/Model_levels/Horizontal_velocity',
                                                      data_label = r'ms$^{-1}$', domain = 'subset1', 
                                                      flight = 301, savefig = True, res = '1p5km',
                                                      variable_name = 'Horizontal velocity',
                                                      contour_type = 'pcolor')
    
            plt.close()
    plt.close()   
#%% subset for region 2
## boundaries (in points):
South = 60; North = 150; West = 130; East = 210
# subset 4d data
sub_data_region2 = modf.subset_4d(magnitude, South, North, West, East)
sub_u_region2 = modf.subset_4d(u, South, North, West, East)
sub_v_region2 = modf.subset_4d(v, South, North, West, East)
       
# two-dimensional subsetting of land-mask
sub_land_region2 = modf.subset_2d(land_mask[0][0], South, North, West, East )

# subset coordinates and unrotate
sub_lon2 = modf.subset_1d(u_lon, West, East)
sub_lat2 = modf.subset_1d(u_lat, South, North)
sub_land_lon2 = modf.subset_1d(land_longitude, West, East)
sub_land_lat2 = modf.subset_1d(land_latitude, South, North)
sub_lons2, sub_lats2 = modf.unrotate_coords(sub_lon2,
                                          sub_lat2,
                                          land_cube.coord('grid_longitude').coord_system.grid_north_pole_longitude,
                                          NPoleLat = land_cube.coord('grid_latitude').coord_system.grid_north_pole_latitude)

#%%
counter = 0
region2_on = False
if region2_on == True:
    for level in range(40):
        for time in range(48):
            counter += 1; 
            print(counter, 'model-level:', level + 1, 'time-step:', time + 1)
            fig, ax = modf.make_horizontal_model_contourf_linear(sub_data_region2, sub_land_region2, sub_lon2, sub_lat2,  sub_lons2, sub_lats2, 
                                                              sub_land_lon2, sub_land_lat2,
                                                      level = level, time = time, colourmap = 'viridis',
                                                      quiver = True, quiver_u = sub_u_region2, quiver_v = sub_v_region2, quiver_n = 5,
                                                      variable_in_file_name = 'horizontal_velocity',
                                                      figure_path_pdf = '../Figures/PDF/301/Model_levels/Horizontal_velocity',
                                                      figure_path_png = '../Figures/PNG/301/Model_levels/Horizontal_velocity',
                                                      data_label = r'ms$^{-1}$', domain = 'subset2', 
                                                      flight = 301, savefig = True, res = '1p5km',
                                                      variable_name = 'Horizontal velocity',
                                                      contour_type = 'pcolor')
    
            plt.close()
    plt.close()




#%% try animation again, but better
# Set up formatting for the movie files -  may need to adjust this - e.g. more or less frames per second.
Writer = animation.writers['ffmpeg']
writer = Writer(fps=1, metadata=dict(artist='Wilhelm Hodder'), bitrate=-1)
figsize = (18,15)


norm_offset = 0; level = 9
## set colourmap and normalise colourmap to data
my_cmap = plt.cm.viridis
my_cmap.set_under('w')

## normalise the colourmap to be between a certain min and max
data_min = np.nanmin(sub_data_region2[:][level]) - norm_offset;
data_max = np.nanmax(sub_data_region2[:][level]) + norm_offset
my_norm= matplotlib.colors.Normalize(vmin = data_min,vmax= data_max, clip=True)

    
#%%
print(len(u[:, level]), len(v))

#%% animate quiver
import animation_functions as anif
animation_on = False
if animation_on == True:
    for j in range(40):
        level = j
        print(f'Animation {level} initiated')
        fig = plt.figure(figsize = figsize)
        ax = fig.add_subplot(111)
        ax.set_aspect('equal',adjustable='box') #ensures the correct plot aspect ratio
        a = ax.pcolor(sub_lon2, sub_lat2, sub_data_region2[0][level], cmap = my_cmap, norm = my_norm, rasterized = True) # pcolor base
        #Add colorbar
        cbar = plt.colorbar( a, ax=ax, orientation='vertical'); cbar.ax.set_xlabel('Windspeed [ms$^{-1}$]')
        fargs = (sub_u_region2, sub_v_region2, 5, sub_lon2, sub_lat2,sub_land_region2, level, 0.5, my_cmap, ax)
        ani = animation.FuncAnimation(fig,anif.animate_pcolor_quiver,48,interval=5,blit=False, fargs = fargs)
        ani.save(f'{figure_path_mp4}/subset2/1p5km_subset2_horizontal_velocity_20180312-301_24hrs_ML{level+1}.mp4', writer=writer, dpi=300)