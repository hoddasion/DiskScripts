# -*- coding: utf-8 -*-
"""
Created on Fri Aug 14 15:33:21 2020

@author: kse18nru
"""

#%% Import modules and scripts

import numpy as np
import matplotlib.pyplot as plt

import iris
import iris.coord_categorisation

import xsection_workshop as workshop

#%% execute function
legs = [1,2,3]
resolutions = ['0p5km', '1p5km', '4p4km']
varname = 'air_pressure'
vmax = 105000; vmin = 80000
#%%
D = '0p5km'; leg = 1
Case = {'res' : D, 'flight' : 301, 'varname' : varname, 'experiment' : 'Control', 'config' : 'RA1M', 'suite' : 'u-bu807', 'filestream' : 'i', 'fileleg' : leg}
res = Case['res']; suite = Case['suite']; exp = Case['experiment']; varname = Case['varname']; config = Case['config']; stream = Case['filestream']; flight = Case['flight']
fileleg = Case['fileleg']

maincube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', varname)[:,:23]
thetacube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
print(maincube, '\n', thetacube, '\n', orogcube)
print(maincube.units)
print('maximum',np.nanmax(maincube.data))
#%%
fig = workshop.plot_vertical_xsection(maincube, thetacube, orogcube,index = 6, 
                                      Case = Case, cmap = 'viridis', vmin = vmin,vmax = vmax, levels = 20,
                                      theta_color = 'r')
plt.show()
#%%
run_loop = True
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
                
                maincube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', varname)[:,:23]
                thetacube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/Vertical_slices/{varname}/{config}_{res}_{varname}_flt{flight}_leg{fileleg}.nc', 'air_potential_temperature')
                orogcube = iris.load_cube(f'../../Model_Data/{suite}/nc/{exp}/surface_altitude/{config}_{res}_surface_altitude_snaesfellness_flt301.nc')
                #print(maincube, '\n', thetacube, '\n', orogcube)
                
            except:
                print('File error')
                continue
            path = f'../../Figures/PNG/301/u-bu807/Vertical/{varname}/leg{fileleg}/{res}'
            for index in range(47):
                
        
                
                fig = workshop.plot_vertical_xsection(maincube, thetacube, orogcube,index = index, 
                                      Case = Case, cmap = 'viridis', vmin = vmin,vmax = vmax, levels = 20,
                                      theta_color = 'r')
                workshop.save_vert_xsection(fig, maincube,index, path, Case, fileleg = leg)
                #plt.close()
                plt.show()