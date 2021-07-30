
"""
Created on Wed Jul  1 10:59:11 2020

@author: kse18nru
"""
import iris
## defining resolution domain, directory, and filenames

var_dir = 'windvector_pl'
res = '4p4km'
suite = 'u-bu807'
nc_path = f'../../Model_Data/{suite}/nc/Control/'
xname = 'x_wind'; yname = 'y_wind'
filenamex = f'{var_dir}/{res}_{xname}_'
filenamey = f'{var_dir}/{res}_{yname}_'
name = 'mag_wind'
filename = f'{var_dir}/{res}_{name}_'
flight = 301
stream = 'verc'
#dump_cube = iris.load_cube(f'{nc_path}{var_dir}/{res}_x_wind_verc1_flt{flight}.nc')
cubex = iris.load_cube(f'{nc_path}{filenamex}24hrs_{stream}_301.nc', 'x_wind')
cubey = iris.load_cube(f'{nc_path}{filenamey}24hrs_{stream}_301.nc', 'y_wind')

cubemag = (cubex**2 + cubey**2)**0.5
cubemag.rename(name)
try:
    delta = cubex.coord('atmosphere_hybrid_height_coordinate') # load hybrid height
    sigma = cubex.coord('sigma') # load sigma
    orography = cubex.coord('surface_altitude') # load orography
    factory = iris.aux_factory.HybridHeightFactory(delta, sigma, orography)
    factory.rename('altitude') # rename derived coordinate from 'unknown' to 'altitude'
    cubemag.add_aux_factory(factory) # integrate into cube
    ## save new cube to file
    iris.save(cubemag, f'{nc_path}{filename}24hrs_{stream}_{flight}.nc')
except:
    print('Exception occured')
    iris.save(cubemag, f'{nc_path}{filename}24hrs_{stream}_{flight}.nc')
    


## test everything was successfull by inspecting cube in new file
testcube = iris.load_cube(f'{nc_path}{filename}24hrs_{stream}_{flight}.nc', name) 
print(testcube)


