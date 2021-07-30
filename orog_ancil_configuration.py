# -*- coding: utf-8 -*-
"""
Created on Fri May 28 15:34:30 2021

@author: kse18nru
"""
#%%
import ants
import iris
#%%
mean_orog_4p4  = ants.fileformats.namelist.load_grid('D:/Project/Model_Data/u-cf117/Ancils/4p4km/qrparm.orog.mn')
print(mean_orog_4p4)
#%%
mean_orog_1p5  = iris.load('D:/Project/Model_Data/u-cf117/Ancils/1p5km/qrparm.orog.mn')
mean_orog_0p5  = iris.load('D:/Project/Model_Data/u-cf117/Ancils/0p5km/qrparm.orog.mn')

#%%
mean_orog_4p4.coord(axis='x').guess_bounds()
mean_orog_4p4.coord(axis='y').guess_bounds()
mean_orog_0p5.coord(axis='x').guess_bounds()
mean_orog_0p5.coord(axis='y').guess_bounds()
mean_orog_1p5.coord(axis='x').guess_bounds()
mean_orog_1p5.coord(axis='y').guess_bounds()
#%%
new_orog_0p5 = mean_orog_4p4.regrid(mean_orog_0p5,iris.analysis.AreaWeighted())
new_orog_1p5 = mean_orog_4p4.regrid(mean_orog_1p5,iris.analysis.AreaWeighted())
#%%
ants.save(new_orog_0p5, 'D:/Project/Model_Data/u-cf117/Ancils/0p5km/qrparm.orog.mn')
ants.save(new_orog_1p5, 'D:/Project/Model_Data/u-cf117/Ancils/1p5km/qrparm.orog.mn')


