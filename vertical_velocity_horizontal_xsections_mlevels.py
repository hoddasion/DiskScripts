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
import foundry
#%% global configurations 

import warnings
warnings.filterwarnings("ignore") ##BAD! If something doesn't work, re-enable warnings to check
font = {'size' : 18}
matplotlib.rc('font', **font)
variable = 'Vertical_velocity'
figure_path_pdf = f'../../Figures/PDF/301/Model_levels/{variable}/Horizontal/'
figure_path_png = f'../../Figures/PNG/301/Model_levels/{variable}/Horizontal/'
figure_path_mp4 = f'../../Figures/MP4/301/Model_levels/{variable}/Horizontal/'
# Turn interactive plotting off
plt.ioff()
#%% load land data
land_file_name = '20180312T0000Z_IcelandGreenlandSeas_1p5km_RA1M_pf000.nc'
land_path = f'../../Model_Data/u-bk574/nc/Control/land_binary_mask/{land_file_name}'

modf.file_info(land_path)

stash = 'land_binary_mask'
land_mask, land_latitude, land_longitude = modf.return_cube_components(land_path, stash, 
                                                        lat_name = 'grid_latitude',
                                                        lon_name = 'grid_longitude')
 
#%% load first variable to analyse: Theta on theta levels
file_name = '1p5km_upward_air_velocity_1.nc'
path1 = f'../../Model_Data/u-bk574/nc/Control/upward_air_velocity/'
stash1 = 'upward_air_velocity'
modf.file_info(f'{path1}{file_name}')#, stash1)
modf.cube_info(f'{path1}{file_name}', stash1)
#modf.level_info(f'{path1}{file_name}', stash1, 3)
model_height = modf.get_model_levels(f'{path1}{file_name}', stash1)
#print(model_height)
arguments = {'Hybrid height' : True}
data, lat, lon = modf.return_concatenated_cube_components(path1, 'i', '1p5km',
                                                          stash1,
                                                          file_number =4,
                                                          lat_name = 'grid_latitude',
                                                          lon_name = 'grid_longitude')

#%% unortate coordinate grid

cube = iris.load_cube(f'{path1}{file_name}', stash1)

lons, lats = foundry.modf.unrotate_coords(lon, lat, 
                                  cube.coord('grid_longitude').coord_system.grid_north_pole_longitude, 
                                  cube.coord('grid_latitude').coord_system.grid_north_pole_latitude)

#%%



#%% full domain plot
time = 28; level = 5
counter = 0
domain_on = False
if domain_on == True:
    for level in range(70):
        for time in range(48):
            counter += 1; print(counter)
            fig = modf.foundry.make_horizontal_model_contourf_linear(data, land_mask[0][0], lon, lat, lons, lats, land_longitude, land_latitude,
                                                             data_label = r'w, [ms$^{-1}$]', domain = 'fulldomain', 
                                                             flight = 301, savefig = True, res = '1p5km',
                                                             variable_name = 'Vertical velocity', level = level, time = time,
                                                             variable_in_file_name = 'upward_air_velocity',
                                                             figure_path_pdf = f'../../Figures/PDF/301/Model_levels/{variable}',
                                                             figure_path_png = f'../../Figures/PNG/301/Model_levels/{variable}',
                                                             contour_type = 'pcolor', colourmap = 'seismic')
            plt.close()


#%% subset for region 1
## boundaries (in points):
South = 60; North = 250; West = 130; East = 230

# subset 4d data
sub_data_region1 = modf.subset_4d(data, South, North, West, East)
       
# two-dimensional subsetting of land-mask
sub_land_region1 = modf.subset_2d(land_mask[0][0], South, North, West, East )

# subset coordinates and unrotate
sub_lon1 = modf.subset_1d(lon, West, East)
sub_lat1 = modf.subset_1d(lat, South, North)
sub_lons1, sub_lats1 = modf.unrotate_coords(sub_lon1,
                                          sub_lat1,
                                          cube.coord('grid_longitude').coord_system.grid_north_pole_longitude,
                                          NPoleLat = cube.coord('grid_latitude').coord_system.grid_north_pole_latitude)


#%% plot for region 1
counter = 0
region1_on = False
if region1_on == True:
    for level in range(40):
        for time in range(48):
            counter += 1
            print(counter)
            fig = modf.foundry.make_horizontal_model_contourf_linear(sub_data_region1, sub_land_region1, sub_lon1, sub_lat1, sub_lons1, sub_lats1, 
                                                                     sub_lon1, sub_lat1,
                                                                     data_label =r'w, [ms$^{-1}$]', domain = 'subset1', 
                                                                     flight = 301, savefig = True, res = '1p5km',
                                                                     variable_name = 'Vertical velocity', level = level, time = time,
                                                                     variable_in_file_name = 'upward_air_velocity',
                                                                     figure_path_pdf = f'../Figures/PDF/301/Model_levels/{variable}',
                                                                     figure_path_png = f'../Figures/PNG/301/Model_levels/{variable}',
                                                                     contour_type = 'pcolor', colourmap = 'seismic')
            plt.close()       

#%% subset for region 2
## boundaries (in points):
South = 60; North = 150; West = 130; East = 210
box_reg1 = np.array([])


# subset 4d data
sub_data_region2 = modf.subset_4d(data, South, North, West, East)
       
# two-dimensional subsetting of land-mask
sub_land_region2 = modf.subset_2d(land_mask[0][0], South, North, West, East )

# subset coordinates and unrotate
sub_lon2 = modf.subset_1d(lon, West, East)
sub_lat2 = modf.subset_1d(lat, South, North)
sub_lons2, sub_lats2 = modf.unrotate_coords(sub_lon2,
                                          sub_lat2,
                                          cube.coord('grid_longitude').coord_system.grid_north_pole_longitude,
                                          NPoleLat = cube.coord('grid_latitude').coord_system.grid_north_pole_latitude)       

#%% plot for region 2
counter = 0
region2_on = False
if region2_on == True:
    for level in range(40):
        for time in range(48):
            counter += 1
            print(counter)
            fig = modf.foundry.make_horizontal_model_contourf_linear(sub_data_region2, sub_land_region2, sub_lon2, sub_lat2, sub_lons2, sub_lats2, 
                                                                     sub_lon2, sub_lat2,
                                                                     data_label = r'w, [ms$^{-1}$]', domain = 'subset2', 
                                                                     flight = 301, savefig = True, res = '1p5km',
                                                                     variable_name = 'Vertical velocity', level = level, time = time,
                                                                     variable_in_file_name = 'upward_air_velocity',
                                                                     figure_path_pdf = f'../Figures/PDF/301/Model_levels/{variable}',
                                                                     figure_path_png = f'../Figures/PNG/301/Model_levels/{variable}',
                                                                     contour_type = 'pcolor', colourmap = 'seismic')
            plt.close()   
            
            
#%% try animation again, but better
# Set up formatting for the movie files -  may need to adjust this - e.g. more or less frames per second.
from matplotlib import animation
Writer = animation.writers['ffmpeg']
writer = Writer(fps=1, metadata=dict(artist='Wilhelm Hodder'), bitrate=-1)
figsize = (18,15)


norm_offset = 0; level = 9
## set colourmap and normalise colourmap to data
my_cmap = plt.cm.seismic
my_cmap.set_under('w')

## normalise the colourmap to be between a certain min and max
data_min = np.nanmin(sub_data_region2[:][level]) - norm_offset;
data_max = np.nanmax(sub_data_region2[:][level]) + norm_offset
my_norm= matplotlib.colors.Normalize(vmin = data_min,vmax= data_max, clip=True)

    


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
        cbar = plt.colorbar( a, ax=ax, orientation='vertical'); cbar.ax.set_xlabel(r'w, [ms$^{-1}$]')
        fargs = (sub_data_region2, sub_lon2, sub_lat2,sub_land_region2, level, 0.5, my_cmap, ax, 'Vertical velocity')
        ani = animation.FuncAnimation(fig,anif.animate_pcolor,48,interval=5,blit=False, fargs = fargs)
        ani.save(f'{figure_path_mp4}/subset2/1p5km_subset2_vertical_velocity_20180312-301_24hrs_ML{level+1}.mp4', writer=writer, dpi=300)