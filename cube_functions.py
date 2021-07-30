# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 13:36:51 2019

@author: kse18nru
"""

import iris

def add_hybrid_height(orography_cube, cubes):
    ''' Add hybrid height to a cube where it's missing - thanks Becky Stretton '''
    regrid_cache = {}
    for cube in cubes:
        print(type(orography_cube), type(cube))
        orog = iris.fileformats.rules._ensure_aligned(regrid_cache, orography_cube, cube)
        orog = cube.regrid(orography_cube, iris.analysis.Linear())
        print(orog)
        new_coord = iris.coords.AuxCoord(orog.data,
                                         orog.standard_name,
                                         orog.long_name,
                                         orog.var_name,
                                         orog.units,
                                         attributes=orog.attributes)
        dims = [cube.coord_dims(src_coord)[0] for src_coord in orog.dim_coords]
        cube.add_aux_coord(new_coord, dims)
        cube.add_aux_factory(iris.aux_factory.HybridHeightFactory(cube.coord('level_height'),
                                                                  cube.coord('sigma'),
                                                                  cube.coord('surface_altitude')))
    return cubes
