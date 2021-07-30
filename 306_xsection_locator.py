# -*- coding: utf-8 -*-
"""
Created on Thu May 13 10:23:03 2021

@author: kse18nru
"""

import matplotlib.pyplot as plt
from iris.analysis.cartography import rotate_pole
import numpy as np
import iris
import iris.coord_categorisation
import iris.plot as iplt
import pandas as pd
import xsection_workshop as workshop
import model_foundry as foundry
import matplotlib.colors as mcolors
#%% load obs data

database = pd.read_csv('../../Obvs_Data/Databases/IGP_flights_database_obs_60s_306.txt', delimiter = ' ')
                
db_legno = np.array(database['legno'])
db_lon = np.array(database['lon'])
db_lat = np.array(database['lat'])

#%% loadvert xsection cubes
cubes = []
legs = [1,2,3,4,6,7,8,9,11,13,14,15]
for i in range(12):
    cube = iris.load_cube(f'D:/Project/Model_Data/u-cc134/Vertical/air_potential_temperature/RA1M_vertslice_0p5km_air_potential_temperature_flt306_leg{legs[i]}.nc', 'air_potential_temperature')
    cubes.append(cube)
    

#%% rotate db coords
compcube = iris.load_cube('D:/Project/Model_Data/u-cc134/RA1M_0p5km_um_upward_air_velocity_24hrs_ph_306.nc', 'upward_air_velocity')[0,0]
polelat = compcube.coord('grid_latitude').coord_system.grid_north_pole_latitude
polelon = compcube.coord('grid_longitude').coord_system.grid_north_pole_longitude
rot_db_lon, rot_db_lat = rotate_pole(db_lon, db_lat, polelon, polelat)
rot_db_lon = rot_db_lon + 360
db_coor = rot_db_lat
orogcube = iris.load_cube('D:/Project/Model_Data/u-cc134/RA1M_0p5km_um_land_binary_mask_24hrs_pi_306.nc', 'land_binary_mask')
orog_data = orogcube.data
orog_lon = orogcube.coords('grid_longitude')[0].points
orog_lat = orogcube.coords('grid_latitude')[0].points


#%% commence plotting

fig, ax = plt.subplots(1,1, figsize = (10,16))

ax.contour(orog_lon, orog_lat, orog_data, colors = 'k')
colors = np.array(['firebrick', 'springgreen', 'royalblue', 'gold','crimson','yellow','blueviolet','indigo','peachpuff','fuchsia','cadetblue','forestgreen','darkorange','salmon','aqua'])

for i in range(15):
    ax.scatter(rot_db_lon[db_legno == i + 1], rot_db_lat[db_legno == i + 1], color = colors[i], alpha = 0.8, label = f'leg {i+1}')
#ax.scatter(rot_db_lon, rot_db_lat, c= colors[db_legno-1])           
plt.legend()   
for i in range(len(legs)):
    print(i + 1)
    xsec_lon = cubes[i].coords('grid_longitude')[0].points
    xsec_lat = cubes[i].coords('grid_latitude')[0].points
    ax.plot(xsec_lon,xsec_lat, color = colors[legs[i]-1])
        
fig.suptitle(f'Vert. X-section and Obs point locations')
fig.tight_layout()

ax.set_xticks([])
ax.set_yticks([])
plt.savefig('D:/Project/Figures/PNG/306/u-cc134/xsection_obs_locations_basic.png')
plt.show()

leglocs = True
if leglocs:
    revert_counter = 0
    for i in range(15):
        fig, ax = plt.subplots(1,1, figsize = (10,16))
    
        ax.contour(orog_lon, orog_lat, orog_data, colors = 'k')
        colors = np.array(['firebrick', 'springgreen', 'royalblue', 'gold','crimson','yellow','blueviolet','indigo','peachpuff','fuchsia','cadetblue','forestgreen','darkorange','salmon','aqua'])
        
        
        ax.scatter(rot_db_lon[db_legno == i + 1], rot_db_lat[db_legno == i + 1], color = colors[i], alpha = 0.8, label = f'leg {i+1}')
        #ax.scatter(rot_db_lon, rot_db_lat, c= colors[db_legno-1])           
        #plt.legend()   
        if i + 1 in legs:
            xsec_lon = cubes[revert_counter].coords('grid_longitude')[0].points
            xsec_lat = cubes[revert_counter].coords('grid_latitude')[0].points
            ax.plot(xsec_lon,xsec_lat, color = colors[i])
            
        fig.suptitle(f'Vert. X-section and Obs point location - Leg {i+1}')
        fig.tight_layout()
        
        ax.set_xticks([])
        ax.set_yticks([])
        plt.savefig(f'D:/Project/Figures/PNG/306/u-cc134/leglocs/xsection_obs_locations_leg{i+1}.png')
        if i + 1 in legs:
            revert_counter += 1
        plt.close()