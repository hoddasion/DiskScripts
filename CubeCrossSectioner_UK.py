import numpy as np
import iris
from iris.analysis import trajectory
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
from matplotlib import mlab,colors,cm
from matplotlib.font_manager import FontProperties
import matplotlib.patheffects as PathEffects
from scipy import interpolate,stats
#import scipy
import os
import sys
from copy import deepcopy
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from math import sin, cos, sqrt, atan2, log10
import datetime
from falklands_data_reader import get_aws_header
#import argparse
#import trui
'''
Plot vertical or horizontal slices through data cubes

Written by Peter Sheridan, Met Office, UK 
'''

#p = argparse.ArgumentParser()
sopdatapath='/data/local/fran/OLD_PROJECTS/FALKLANDS/SURFACE_DATA/10min_sop_data/'


def find_closest_point(loclatlondict,xvals,orogslice):
    ''' 
    Crude minimum distance of a gridded line from a point, returns position along the line in 
    distance and point number  
    '''
    iout,xout=[],[] 
    for location in loclatlondict.keys():
        lldist2=(orogslice.coord('grid_longitude').points-loclatlondict[location][0])**2 + \
                (orogslice.coord('grid_latitude').points-loclatlondict[location][1])**2
        iclosest=np.argmin(lldist2)
        xout+=[xvals[iclosest]]
        iout+=[iclosest]
    
    return iout,xout
        

def get_station_positions(locations):
    ''' get dictionary of long lat from list of station numbers using AWS file headers''' 
    stnnums,filenums=np.loadtxt('station_key.txt',skiprows=1,dtype='int').T
    stndict={}
    for i,stnnum in enumerate(stnnums): stndict[stnnum]='%02i' % filenums[i]
    latlondict={}
    locations=[int(location) for location in locations]
    for location in locations:
        sopdatafile=sopdatapath+'10min_sop_data.b'+stndict[location]
        headerdict=get_aws_header(sopdatafile)
        latlondict[location]=[headerdict['station_longitude'], headerdict['station_latitude']]

    return locations,latlondict


def get_locs_from_file(locationsfile):
    ''' read from a file of name,lon,lat, locations to plot as points on cross-sections '''
    fopen=open(locationsfile,'r')
    filelines=fopen.readlines()
    markname,marklon,marklat=[],[],[]
    for fline in filelines:
        markname+=[fline.split()[0]]
        marklon+=[float(fline.split()[1])]
        marklat+=[float(fline.split()[2])]

    latlondict={}
    for name in markname: latlondict[name]=[marklon,marklat]

    return markname,latlondict


def getslicecompt(uslice,vslice,lonstart,latstart,lonend,latend,resolvedir=None):
    '''
    Convert horizontal velocities into the component in-plane using the initial bearing of 
    a vertical slice, or defined by resolvedir 
    ''' 
    dlong=(lonend-lonstart)*np.pi/180.
    if resolvedir: slicedir=np.pi*resolvedir/180.
    else: slicedir=np.arctan2( np.sin(dlong)*np.cos(latend*np.pi/180.) , \
                               np.cos(latstart*np.pi/180.)*np.sin(latend*np.pi/180.) -
                               np.sin(latstart*np.pi/180.)*cos(latend*np.pi/180.)*cos(dlong) )

    udir=90.*np.pi/180.
    vdir=0.
    comptslice=deepcopy(uslice)
    comptslice.data=uslice.data*np.cos(udir-slicedir) + vslice.data*np.cos(vdir-slicedir)

    return comptslice 


def roundit(floatNumber,updown):
    """ 
    rounds numbers up or down as specified
    """
    if updown == 'up':
        if round(floatNumber) < floatNumber:
            return int(round(floatNumber)+1)
        else:
            return int(round(floatNumber))
    elif updown == 'down':
        if round(floatNumber) > floatNumber:
            return int(round(floatNumber - 1))
        else:
            return int(round(floatNumber))


def get_cont_lims(axis1,axis2,xsecdata,xsec,axlims):
    """ 
    Get contour limits given coords and values
    """
    mindat=100000.
    maxdat=-100000.
    if len(axis1.shape) == 2:
        iax=[ i for i in range(len(axis1[0])) if axis1[0,i] >= axlims[0][0] and axis1[0,i] <= axlims[0][1] ]
        jax=[ j for j in range(len(axis2)) if axis2[j,len(axis2[0])/2] >= axlims[1][0] and axis2[j,len(axis2[0])/2] <= axlims[1][1] ]
    else:
        iax=[ i for i in range(len(axis1)) if axis1[i] >= axlims[0][0] and axis1[i] <= axlims[0][1] ]
        jax=[ j for j in range(len(axis2)) if axis2[j] >= axlims[1][0] and axis2[j] <= axlims[1][1] ]
    for i in iax:
        for j in jax:
            maxdat=max([maxdat,xsecdata[j][i]])
            mindat=min([mindat,xsecdata[j][i]])
    datrange=maxdat-mindat
    nmin=5
    nmax=2*nmin
    powr=1+roundit(log10(nmin/datrange),'down')
    facs=(1.0,2.5,5.)
    fac=facs[[i for i in range(len(facs)) if nmin < datrange*10**powr/facs[i] <= nmax][0]]
    div=fac/10**powr
    mindat=div*roundit(mindat/div,'down')
    maxdat=div*roundit(maxdat/div,'up')
    contlevs=np.arange(mindat,maxdat+div,div)

    return(contlevs)


def tickspacecalc(extent):
    '''calculate tick separation'''
    extent=float(extent)
    nmin=4
    nmax=8#5*nmin/2
    powr=min([1+roundit(log10(nmin/extent),'down'),1+roundit(log10(nmax/extent),'down')])
    facs=(1.,2.,3.,5.)
    fac=facs[[i for i in range(len(facs)) if nmin <= extent*10**powr/facs[i] <= nmax][-1]]
    if fac == 3.: fac=facs[[i for i in range(len(facs)) if nmin <= extent*10**powr/facs[i] <= nmax][0]]
    div=fac/10**powr
    mindat=0.
    maxdat=extent
    mindat=div*roundit(mindat/div,'up')
    maxdat=div*roundit(maxdat/div,'down')
    ticks=np.arange(mindat,maxdat+div,div)

    return ticks


def lininterp(coord,data,coordval):
    '''simple lightweight 1-D linear interpolation''' 
    ibottom=np.max(np.where(coord <= coordval))
    coordstep=coord[ibottom+1]-coord[ibottom]
    dcoord=coordval-coord[ibottom]
    ddata=data[ibottom+1]-data[ibottom]
    outdatum=dcoord*ddata/coordstep

    return outdatum


def trajectory_to_lonlat(trajfile):
    ''' Converts a cartesian distance trajectory to a trajectory in lon/lat '''

    R = 6371008.8 # mean radius of earth in m
    
    opentraj=open(trajfile,'r')
    releaselon,releaselat=[float(val) for val in opentraj.readline().split(' ')[2:4]]
    dateline=opentraj.readline()
    opentraj.close()
    ryy,rmm,rdd=[int(val) for val in dateline.split(' ')[2].split('/')]
    rhh,rmn=[int(val) for val in dateline.split(' ')[3].split(':')]
    timesecs,x,y,z=np.loadtxt(trajfile,skiprows=3).T
    
    trajtime=[datetime.datetime(ryy,rmm,rdd,rhh,rmn)]    
    trajlon=[releaselon]    
    trajlat=[releaselat]    
    for i,timenow in enumerate(timesecs[1:],1):
        dlat=(180./np.pi)*(y[i]-y[i-1]) / R
        trajlat+=[trajlat[i-1]+dlat]
        meanlat=(trajlat[i]+trajlat[i-1])/2.
        dlon=(180./np.pi)*(x[i]-x[i-1]) / (R*cos(np.pi*meanlat/180.))
        trajlon+=[trajlon[i-1]+dlon]
        trajtime+=[ trajtime[i-1] + datetime.timedelta(seconds=timenow-timesecs[i-1]) ]
                                                      
    return timesecs,trajtime,trajlon,trajlat,z


def haversine(lat1,lon1,lat2,lon2):
    '''distance between points on a great circle'''

    R = 6371.0088 # mean radius of earth (km) 

    lat1,lon1,lat2,lon2=[element*np.pi/180. for element in lat1,lon1,lat2,lon2] 
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = (sin(dlat/2))**2 + cos(lat1) * cos(lat2) * (sin(dlon/2))**2
    c = 2 * atan2(sqrt(a), sqrt(1-a))
    distance = R * c

    return distance


def seghsvcmap(hsvcols,segvals):
    '''
    defines a segmented colour map
    hsvcols: numpy array of HSV colours bounding segments of the colour map, e.g.
             np.array([[[140.,0.,1.],[80.,0.8,1.],[0.,1.,1.],[0.,1.,0.4]]]) for 3 segments with 4 bounds
    segvals: location of these bounds between 0 and 1 on the colour scale, starting and ending with 0 and 1 themselves
             and have the same number of elements as there are colours defined in hsvcols
    '''
    #hsvcols[:,:,0]=hsvcols[:,:,0]/360.
    #rgbcols=colors.hsv_to_rgb(hsvcols)
    #cdict = {'red':   [(segval,  rgbcols[0,i,0], rgbcols[0,i,0]) for i,segval in enumerate(segvals)],
    #         'green': [(segval,  rgbcols[0,i,1], rgbcols[0,i,1]) for i,segval in enumerate(segvals)],
    #         'blue':  [(segval,  rgbcols[0,i,2], rgbcols[0,i,2]) for i,segval in enumerate(segvals)]}
    #colmap=colors.LinearSegmentedColormap('velocity',cdict)

    finesegs=np.linspace(0,1.,101)
    hsv=[interpolate.griddata(segvals,hsvcols[:,0],finesegs),\
         interpolate.griddata(segvals,hsvcols[:,1],finesegs),\
         interpolate.griddata(segvals,hsvcols[:,2],finesegs)]
    hsvfine=np.array([np.array(hsv).T])
    hsvfine[:,:,0]=hsvfine[:,:,0]/360.
    rgbcols=colors.hsv_to_rgb(hsvfine)
    cdict = {'red':   [(segval,  rgbcols[0,i,0], rgbcols[0,i,0]) for i,segval in enumerate(finesegs)],
             'green': [(segval,  rgbcols[0,i,1], rgbcols[0,i,1]) for i,segval in enumerate(finesegs)],
             'blue':  [(segval,  rgbcols[0,i,2], rgbcols[0,i,2]) for i,segval in enumerate(finesegs)]}
    colmap=colors.LinearSegmentedColormap('velocity',cdict)

    return colmap


def velcmap():
    """
    defines a colour map for velocities, 2/3 of the map is for positive velocities, 1/3 for negative
    """
    getmap=cm.ScalarMappable(cmap='jet').get_cmap()
    poscols=[getmap(3*i/2) for i in range(2*255/3)]
    getmap=cm.ScalarMappable(cmap='bone').get_cmap()
    negcols=[getmap(255-2*i) for i in range(255-2*255/3)]

    colmap=colors.LinearSegmentedColormap.from_list('velocity',negcols+poscols)

    return colmap


def wcmap():
    """
    defines a colour map for vertical velocities
    """
    getmap=cm.ScalarMappable(cmap='RdBu_r').get_cmap()
    getwhite=cm.ScalarMappable(cmap='Greys').get_cmap()
    cols=[getmap(i) for i in range(255)]
    cols[256*7/16:256*9/16]=[getwhite(0)]*len(cols[256*7/16:256*9/16])

    colmap=colors.LinearSegmentedColormap.from_list('vertvel',cols)

    return colmap


def vom_cmap():
        ''' gets w colour scheme as 3DVOM '''
        colsegs=[0,15./32.,15./32.+0.001,17./32.-0.001,17./32.,1.]
        hsvcols=np.array([[240.,1.,0.5],[120.,1.,1.],[120.,0.,1.],[60.,0.,1.],[60.,1.,1.],[0.,1.,0.5]])
        colmap=seghsvcmap(hsvcols,colsegs)
        
        return colmap


def add_hybrid_height(orography_cube, cubes):
    ''' Add hybrid height to a cube where it's missing - thanks Becky Stretton '''
    regrid_cache = {}
    for cube in cubes:
        orog = iris.fileformats.rules._ensure_aligned(regrid_cache, orography_cube, cube)
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


def get_1D_profile(thppname,velppname,oppname,ptlon0,ptlat0,ppyy,ppmm,ppdd,pphh,outfstem):
    ''' Output vertical model profile above a point at a given time '''
     
    ppyy,ppmm,ppdd,pphh,ptlon0,ptlat0=int(ppyy),int(ppmm),int(ppdd),int(pphh),float(ptlon0),float(ptlat0)

    thcube=iris.load_cube(thppname,iris.AttributeConstraint(STASH='m01s00i004'))
    wcube=iris.load_cube(velppname,iris.AttributeConstraint(STASH='m01s00i150')) 
    ucube=iris.load_cube(velppname,iris.AttributeConstraint(STASH='m01s00i002')) 
    vcube=iris.load_cube(velppname,iris.AttributeConstraint(STASH='m01s00i003')) 
    ocube=iris.load_cube(oppname,iris.AttributeConstraint(STASH='m01s00i033'))
    ocube=ocube.regrid(thcube,iris.analysis.Linear())
    ucube=ucube.regrid(thcube,iris.analysis.Linear())
    vcube=vcube.regrid(thcube,iris.analysis.Linear())
    thcube,wcube,ucube,vcube=add_hybrid_height(ocube,[thcube,wcube,ucube,vcube])

    polelon=ocube.coord_system('RotatedGeogCS').grid_north_pole_longitude
    polelat=ocube.coord_system('RotatedGeogCS').grid_north_pole_latitude

    ppdatetime=datetime.datetime(ppyy,ppmm,ppdd,pphh)

    with iris.FUTURE.context(cell_datetime_objects = True):
        timeconstraint=iris.Constraint(time=lambda cell: cell.point == ppdatetime)
        thcube=thcube.extract(timeconstraint)
        ucube=ucube.extract(timeconstraint)
        vcube=vcube.extract(timeconstraint)
        wcube=wcube.extract(timeconstraint)

    if ocube.coord_system('RotatedGeogCS') != None and polelon != 0 and polelat != 90: rotated=True

    if rotated: 
        ptlon,ptlat=iris.analysis.cartography.rotate_pole(np.array(ptlon0),np.array(ptlat0),polelon,polelat)
        thpt=thcube.interpolate([('grid_latitude',ptlat),('grid_longitude',ptlon)], iris.analysis.Linear())
        upt=ucube.interpolate([('grid_latitude',ptlat),('grid_longitude',ptlon)], iris.analysis.Linear())
        vpt=vcube.interpolate([('grid_latitude',ptlat),('grid_longitude',ptlon)], iris.analysis.Linear())
        wpt=wcube.interpolate([('grid_latitude',ptlat),('grid_longitude',ptlon)], iris.analysis.Linear())

    else:
        ptlon,ptlat=ptlon0,ptlat0
        thpt=thcube.interpolate([('latitude',ptlat),('longitude',ptlon)], iris.analysis.Linear())
        upt=ucube.interpolate([('latitude',ptlat),('longitude',ptlon)], iris.analysis.Linear())
        vpt=vcube.interpolate([('latitude',ptlat),('longitude',ptlon)], iris.analysis.Linear())
        wpt=wcube.interpolate([('latitude',ptlat),('longitude',ptlon)], iris.analysis.Linear())
        
    thpt,upt,vpt,wpt=iris.util.squeeze(thpt),iris.util.squeeze(upt),iris.util.squeeze(vpt),iris.util.squeeze(wpt)
    outdat=np.ndarray((len(thpt.data),6))
    datalist=thpt.coord('altitude').points,wpt.coord('altitude').points,thpt.data,upt.data,vpt.data,wpt.data
    for i in range(6): outdat[:,i]=np.array(datalist[i])
    
    outfname=outfstem+'_%7.3f_%7.3f.dat' % (ptlon0,ptlat0)
    np.savetxt(outfname,outdat)#,fmt='%5i %8i %8.3 %7.2 %7.2 %6.2')



def trace_1D_trajectory(thppname,velppname,oppname,trajfile,ppyy,ppmm,ppdd,pphh,outfname):
    ''' Trace a 1D trajectory (e.g. radiosonde) through model and interpolate data fields ''' 
    
    timesecs,trajtime,trajlon,trajlat,z=trajectory_to_lonlat(trajfile)
    trajtime,trajlon,trajlat=np.array(trajtime),np.array(trajlon),np.array(trajlat)
    ppyy,ppmm,ppdd,pphh=int(ppyy),int(ppmm),int(ppdd),int(pphh)
    
    thcube=iris.load_cube(thppname,iris.AttributeConstraint(STASH='m01s00i004'))
    wcube=iris.load_cube(velppname,iris.AttributeConstraint(STASH='m01s00i150')) 
    ucube=iris.load_cube(velppname,iris.AttributeConstraint(STASH='m01s00i002')) 
    vcube=iris.load_cube(velppname,iris.AttributeConstraint(STASH='m01s00i003')) 
    ocube=iris.load_cube(oppname,iris.AttributeConstraint(STASH='m01s00i033'))
    ocube=ocube.regrid(thcube,iris.analysis.Linear())
    ucube=ucube.regrid(thcube,iris.analysis.Linear())
    vcube=vcube.regrid(thcube,iris.analysis.Linear())
    thcube,wcube,ucube,vcube=add_hybrid_height(ocube,[thcube,wcube,ucube,vcube])
    
    polelon=ocube.coord_system('RotatedGeogCS').grid_north_pole_longitude
    polelat=ocube.coord_system('RotatedGeogCS').grid_north_pole_latitude

    ppdatetime=datetime.datetime(ppyy,ppmm,ppdd,pphh)

    with iris.FUTURE.context(cell_datetime_objects = True):
        timeconstraint=iris.Constraint(time=lambda cell: cell.point == ppdatetime)
        thcube=thcube.extract(timeconstraint)
        ucube=ucube.extract(timeconstraint)
        vcube=vcube.extract(timeconstraint)
        wcube=wcube.extract(timeconstraint)
    
    if ocube.coord_system('RotatedGeogCS') != None and polelon != 0 and polelat != 90: rotated=True

    if rotated: 
        trajlon,trajlat=iris.analysis.cartography.rotate_pole(np.array(trajlon),np.array(trajlat),polelon,polelat)
        iltzero=np.where(trajlon < 0)
        trajlon[iltzero]+=360.

        idomain=np.where( (trajlon <= thcube.coord('grid_longitude').points[-1]) & \
                          (trajlon >= thcube.coord('grid_longitude').points[0]) & \
                          (trajlat <= thcube.coord('grid_latitude').points[-1]) & \
                          (trajlat >= thcube.coord('grid_latitude').points[0]) & \
                          (z <= np.max(thcube.coord('altitude').points) ) )

        timesecs,trajtime,trajlon,trajlat,z = \
            timesecs[idomain],trajtime[idomain],trajlon[idomain],trajlat[idomain],z[idomain]
        
        thslice=trajectory.interpolate(thcube, \
                                      (('grid_longitude', trajlon), ('grid_latitude', trajlat)), \
                                      method='linear')        
        uslice=trajectory.interpolate(ucube, \
                                      (('grid_longitude', trajlon), ('grid_latitude', trajlat)), \
                                      method='linear')        
        vslice=trajectory.interpolate(vcube, \
                                      (('grid_longitude', trajlon), ('grid_latitude', trajlat)), \
                                      method='linear')        
        wslice=trajectory.interpolate(wcube, \
                                      (('grid_longitude', trajlon), ('grid_latitude', trajlat)), \
                                      method='linear')        
    
    else: 
        idomain=np.where( (trajlon <= thcube.coord('longitude').points[-1]) & \
                          (trajlon >= thcube.coord('longitude').points[0]) & \
                          (trajlat <= thcube.coord('latitude').points[-1]) & \
                          (trajlat >= thcube.coord('latitude').points[0]) & \
                          (z <= thcube.coord('altitude').points[-1]) )

        timesecs,trajtime,trajlon,trajlat,z = \
            timesecs[idomain],trajtime[idomain],trajlon[idomain],trajlat[idomain],z[idomain]

        thslice=trajectory.interpolate(thcube, \
                                      (('longitude', trajlon), ('latitude', trajlat)), \
                                      method='linear')        
        uslice=trajectory.interpolate(ucube, \
                                      (('longitude', trajlon), ('latitude', trajlat)), \
                                      method='linear')        
        vslice=trajectory.interpolate(vcube, \
                                      (('longitude', trajlon), ('latitude', trajlat)), \
                                      method='linear')        
        wslice=trajectory.interpolate(wcube, \
                                      (('longitude', trajlon), ('latitude', trajlat)), \
                                      method='linear')        

    thtraj,utraj,vtraj,wtraj=np.ndarray((len(z),)),np.ndarray((len(z),)),np.ndarray((len(z),)),np.ndarray((len(z),))
    for i,zval in enumerate(z):
        thtraj[i]=np.interp(zval,thslice.coord('altitude')[:,i].points , thslice.data[:,i]) 
        utraj[i]=np.interp(zval,uslice.coord('altitude')[:,i].points , uslice.data[:,i]) 
        vtraj[i]=np.interp(zval,vslice.coord('altitude')[:,i].points , vslice.data[:,i]) 
        wtraj[i]=np.interp(zval,wslice.coord('altitude')[:,i].points , wslice.data[:,i]) 
    

    outdat=np.ndarray((len(timesecs),6))
    datalist=timesecs,z,thtraj,utraj,vtraj,wtraj
    for i in range(6): outdat[:,i]=np.array(datalist[i])
    
    np.savetxt(outfname,outdat)#,fmt='%5i %8i %8.3 %7.2 %7.2 %6.2')
    


def plot_hoz_xsec(velppnames,oppname,sliceheight,swlonin,swlatin,nelonin,nelatin,velcontmax,prefix,affix,outdir,\
                  slicedttm=None,locate=None,locations=None,locfile=None,drawpaths=None,pathcols=None,linestyles=None,\
                  orogint=100.,ocolour='k',coastcolour='k',vectcol='k',xlims=None,ylims=None):
    '''
    Plot a horizontal cross section zoomed on part of a domain for a known file type 
    and specific plot set-up:

    "plot_hoz_xsec(velppnames,oppname,sliceheight,swlonin,swlatin,nelonin,nelatin,velcontmax,prefix,affix,outdir,\
                  slicedttm=None,locate=None,locations=None,locfile=None,drawpaths=None,pathcols=None,linestyles=None,\
                  orogint=100.,ocolour='k',coastcolour='k',vectcol='k',xlims=None,ylims=None):"   
    velppnames: list of model velocity  data filenames - will plot all times unless a specific time is given by slicedttm 
    oppname: model orography data filename
    slicedttm: list containing [year,month,day,hour,minutes]
    sliceheight: height of desired horizontal slice - used to choose a single level to plot (no interpolation)
    swlonin,swlatin,nelonin,nelatin: location of slice corners
    velcontmax: colour contour maximum
    prefix: prefix to add to image filenames
    affix: affix to add to image filenames
    outdir: directory to output image files
    locations: AWS numbers of which to plot location  
    locate: 'list' indicates locations are instead listed in a file given by 'locfile'
    drawpaths: lines to draw, list of numpy arrays of lon/lat (e.g. drawpaths[0][0,:] is longitude of first path)
    pathcols: colours of drawn paths (default black)
    linestyles: line styles of drawn paths (default dashed)
    orogint: orography contour interval 
    ocolour: orography contour colour
    coastcolour: coastline colour
    vectcol: wind vector colour
    xlims: x axis limits in unrotated longitude
    ylims: y axis limits in unrotated latitude
    '''
    
    affix='%06.2f_%06.2f_to_%06.2f_%06.2f_%05im' % ( swlonin,swlatin,nelonin,nelatin,sliceheight ) + affix

    if locate == 'list': locations,loclatlondict=get_locs_from_file(locfile)
    elif locations: locations,loclatlondict=get_station_positions(locations)
    else: loclatlondict=None

    ocube=iris.load_cube(oppname,iris.AttributeConstraint(STASH='m01s00i033'))#/data/local/fran/COLPEX/MASS-R_RETRIEVALS/orog_UKV.pp')

    polelon=ocube.coord_system('RotatedGeogCS').grid_north_pole_longitude
    polelat=ocube.coord_system('RotatedGeogCS').grid_north_pole_latitude
    if ocube.coord_system('RotatedGeogCS') != None and polelon != 0 and polelat != 90: rotated=True

    print 'corners',(swlonin,nelonin),(swlatin,nelatin)
    
    if rotated: 
        (swlon,nelon),(swlat,nelat)=iris.analysis.cartography.rotate_pole(np.array([swlonin,nelonin]),np.array([swlatin,nelatin]),polelon,polelat)
#        if swlon < 0: swlon+=360. 
#        if nelon < 0: nelon+=360. 
        lon_extent=iris.coords.CoordExtent(ocube.coord('grid_longitude'),swlon,nelon)
        lat_extent=iris.coords.CoordExtent(ocube.coord('grid_latitude'),swlat,nelat)
        if locations:
            for location in locations:
                loclatlondict[location][0],loclatlondict[location][1]=iris.analysis.cartography.rotate_pole(\
                 np.array([loclatlondict[location][0]]),np.array([loclatlondict[location][1]]),polelon,polelat)
                if loclatlondict[location][0] < 0: loclatlondict[location][0]+=360.
        if drawpaths:
            for i,drawpath in enumerate(drawpaths):
                drawpaths[i][0,:],drawpaths[i][1,:]=\
                 iris.analysis.cartography.rotate_pole(drawpath[0,:],drawpath[1,:],polelon,polelat)
                ilonover=np.where(drawpaths[i][0,:] < 0)
                drawpaths[i][0,ilonover]+=360.

    print 'rotated corners',(swlon,nelon),(swlat,nelat)

    orog_zoom=ocube.intersection(lon_extent,lat_extent)
    ocontlevs=np.linspace(orogint,10000,10000/orogint)

    for velppname in velppnames:

        print 'velppname',velppname
        
        try: ucube=iris.load_cube(velppname,iris.AttributeConstraint(STASH='m01s00i002'))
        except: ucube=iris.load(velppname,iris.AttributeConstraint(STASH='m01s00i002'))[0]
        try: vcube=iris.load_cube(velppname,iris.AttributeConstraint(STASH='m01s00i003'))
        except: vcube=iris.load(velppname,iris.AttributeConstraint(STASH='m01s00i003'))[0]
        try: wcube=iris.load_cube(velppname,iris.AttributeConstraint(STASH='m01s00i150'))
        except: wcube=iris.load(velppname,iris.AttributeConstraint(STASH='m01s00i150'))[0]
    
        wzoom=wcube.regrid(orog_zoom,iris.analysis.Linear())
        uzoom=ucube.regrid(orog_zoom,iris.analysis.Linear())
        vzoom=vcube.regrid(orog_zoom,iris.analysis.Linear())
        wzoom,uzoom,vzoom=add_hybrid_height(orog_zoom,[wzoom,uzoom,vzoom])
    
        print 'cubes regridded to orography'

        ilev=np.argmin((wzoom.coord('level_height').points-sliceheight)**2)
        levelconstraint=iris.Constraint(model_level_number=ilev+1)
        if slicedttm:
            slicedttm=datetime.datetime(slicedttm[0],slicedttm[1],slicedttm[2],slicedttm[3],slicedttm[4])
            print 'slicedttm',slicedttm

            with iris.FUTURE.context(cell_datetime_objects = True):
                # nearest level to desired height
                timeconstraint=iris.Constraint(time=lambda cell: cell.point == slicedttm)
                wslice=wzoom.extract(levelconstraint & timeconstraint)
                uslice=uzoom.extract(levelconstraint & timeconstraint)
                vslice=vzoom.extract(levelconstraint & timeconstraint)
                validdttm=wslice.coord('time').units.num2date(uslice.coord('time').points[0]).strftime("%Y%m%d%H%M")
            create_hoz_xsec_plot(wslice,uslice,vslice,orog_zoom,ocontlevs,velcontmax,linestyles,drawpaths,pathcols,ocolour,coastcolour,vectcol,\
                                 polelon,polelat,swlonin,nelonin,swlatin,nelatin,xlims,ylims,validdttm,outdir,prefix,affix,locations,loclatlondict)

        else: 
            # plot all slices in the file
            wlvl=wzoom.extract(levelconstraint)
            ulvl=uzoom.extract(levelconstraint)
            vlvl=vzoom.extract(levelconstraint)
            for i,wslice in enumerate(wlvl.slices(['grid_latitude', 'grid_longitude'])):
                if len(ulvl.shape) == 3: uslice,vslice=ulvl[i],vlvl[i]
                else: uslice,vslice=ulvl,vlvl
                validdttm=wslice.coord('time').units.num2date(uslice.coord('time').points[0]).strftime("%Y%m%d%H%M")
                create_hoz_xsec_plot(wslice,uslice,vslice,orog_zoom,ocontlevs,velcontmax,linestyles,drawpaths,pathcols,ocolour,coastcolour,vectcol,\
                                     polelon,polelat,swlonin,nelonin,swlatin,nelatin,xlims,ylims,validdttm,outdir,prefix,affix,locations,loclatlondict)



def create_hoz_xsec_plot(wslice,uslice,vslice,orog_zoom,ocontlevs,velcontmax,linestyles,drawpaths,pathcols,ocolour,coastcolour,vectcol,\
                         polelon,polelat,swlonin,nelonin,swlatin,nelatin,xlims,ylims,validdttm,outdir,prefix,affix,locations,loclatlondict):
        
    print 'wslice,uslice,vslice',wslice,uslice,vslice
    if velcontmax == None: velcontmax=( int(np.max(np.abs(wslice.data))) / 2 ) * 2
    wcontlevs=np.linspace(-velcontmax,velcontmax,8*velcontmax+1)

    meanlon=(swlonin+nelonin)/2.
    meanlat=(swlatin+nelatin)/2.
    ewlength=haversine(meanlat,swlonin,meanlat,nelonin)
    nslength=haversine(swlatin,meanlon,nelatin,meanlon)
    figsize=( 7.5*ewlength/stats.mstats.gmean([ewlength,nslength]) , 7.*nslength/stats.mstats.gmean([ewlength,nslength]) )

    fig=plt.figure(figsize=figsize)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines('10m', color = coastcolour, linewidth=1.)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linestyle = ':', draw_labels = True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    theplot=iplt.contourf(wslice,levels=wcontlevs,cmap=vom_cmap(),extend="both")
    iplt.contour(orog_zoom,levels=ocontlevs,colors=ocolour,linewidths=0.5)
    iplt.contour(orog_zoom,levels=[0,10000],colors=ocolour,linewidths=1)

    xvals,yvals=orog_zoom.coord('grid_longitude').points,orog_zoom.coord('grid_latitude').points
    
    ax.set_aspect( (nslength/(ylims[1]-ylims[0])) / (ewlength/(xlims[1]-xlims[0])) )
    
    axextentinches=ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    axleftedge=axextentinches.x0/fig.get_figwidth()
    axrightedge=axextentinches.x1/fig.get_figwidth()
    axbottedge=axextentinches.y0/fig.get_figheight()
    axtopedge=axextentinches.y1/fig.get_figheight()

    nvect=(15,20)
    nskipx,nskipy=max([uslice.data.shape[1]/nvect[0],1]),np.max([uslice.data.shape[0]/nvect[1],1])
    lon,lat = orog_zoom.coord('grid_longitude').points[::nskipx],orog_zoom.coord('grid_latitude').points[::nskipy]
    u,v = uslice.data[::nskipy,::nskipx],vslice.data[::nskipy,::nskipx]
    rotated_pole = ccrs.RotatedPole(pole_longitude=polelon, pole_latitude=polelat)
    ax.set_extent((swlonin,nelonin,swlatin,nelatin), crs=ccrs.PlateCarree())
    vectscale=nvect[0]*20.
    arrowobj=ax.quiver(lon, lat, u, v, transform = rotated_pole, scale = vectscale, scale_units = 'width')

    plt.quiverkey(arrowobj,0.05,-0.13,20.,'20 m/s',color='k')

    if drawpaths:
        print 'adding drawn paths...'
        for i,drawpath in enumerate(drawpaths):
            if pathcols: pathcol=pathcols[i]
            else: pathcol='k'
            if linestyles: linestyle=linestyles[i]
            else: linestyle='-'
            if linestyle == '-': alpha=0.6
            elif linestyle == ':': alpha=1.
            else: alpha=0.8
            if len(drawpath[0,:]) > 2: 
                # implies a trajectory, e.g. sonde
                ax.plot(drawpath[0,:400:5],drawpath[1,:400:5],pathcol,transform = rotated_pole) # has trouble with more than 100 points in some cases
                ax.plot([drawpath[0,0]],[drawpath[1,0]],'k*',markersize=16,markeredgecolor='w',\
                                               markeredgewidth=1,transform = rotated_pole)
            else: ax.plot(drawpath[0,:],drawpath[1,:],pathcol,linestyle=linestyle,linewidth=2.5,transform = rotated_pole,alpha=alpha)


    if locations:
        print 'adding station locations....'
        for location in loclatlondict.keys():
            ax.plot(loclatlondict[location][0],loclatlondict[location][1],'ko',markeredgecolor='w',\
                    markeredgewidth=1,markersize=7,transform = rotated_pole)
            txt=ax.text(loclatlondict[location][0]+(xvals[-1]-xvals[0])/100.,\
                        loclatlondict[location][1]+(xvals[-1]-xvals[0])/100.,\
                        str(location),color='k',fontsize=13,transform = rotated_pole)
            txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='w')])

    if xlims: ax.set_xlim(xlims)
    if ylims: ax.set_ylim(ylims)
    
    rightspace=1.-axrightedge
    colorbar_axes = plt.gcf().add_axes([axrightedge-rightspace/8., axbottedge, rightspace/4., axtopedge-axbottedge])
    bar=plt.colorbar(theplot,colorbar_axes,orientation='vertical')
    barlabel='vertical velocity (ms$^{-1})$'
    bar.set_label(barlabel)

    plt.savefig(outdir+'/'+prefix+'_'+validdttm+'_'+affix+'.png')
    plt.close()



def plot_single_level_vel(velppnames,oppname,swlonin,swlatin,nelonin,nelatin,velcontmax,prefix,affix,outdir,\
                          slicedttm=None,locate=None,locations=None,locfile=None,drawpaths=None,pathcols=None,linestyles=None,\
                          orogint=100.,ocolour='k',coastcolour='k',vectcol='k',xlims=None,ylims=None):
    '''
    Plot single level (e.g. 10 m wind) data zoomed on part of a domain for a known file type 
    and specific plot set-up:
    Plotted with respect to lon/lat coords with terrain contours and hoz vectors using quiver
    Specified area, specified time    

    "plot_single_level_vel(velppnames,oppname,swlonin,swlatin,nelonin,nelatin,velcontmax,prefix,affix,outdir,\
                          slicedttm=None,locate=None,locations=None,locfile=None,drawpaths=None,pathcols=None,linestyles=None,\
                          orogint=100.,ocolour='k',coastcolour='k',vectcol='k',xlims=None,ylims=None)"   

    velppnames: list of model velocity  data filenames - will plot all times unless a specific time is given by slicedttm 
    oppname: model orography data filename
    slicedttm: list containing [year,month,day,hour,minutes]
    swlonin,swlatin,nelonin,nelatin: location of slice corners
    velcontmax: colour contour maximum
    prefix: prefix to add to image filenames
    affix: affix to add to image filenames
    outdir: directory to output image files
    locations: AWS numbers of which to plot location  
    locate: 'list' indicates locations are instead listed in a file given by 'locfile'
    drawpaths: lines to draw, list of numpy arrays of lon/lat (e.g. drawpaths[0][0,:] is longitude of first path)
    pathcols: colours of drawn paths (default black)
    linestyles: line styles of drawn paths (default dashed)
    orogint: orography contour interval 
    ocolour: orography contour colour
    coastcolour: coastline colour
    vectcol: wind vector colour
    xlims: x axis limits in unrotated longitude
    ylims: y axis limits in unrotated latitude
    '''

    affix='%06.2f_%06.2f_to_%06.2f_%06.2f' % ( swlonin,swlatin,nelonin,nelatin ) + affix

    if locate == 'list': locations,loclatlondict=get_locs_from_file(locfile)
    elif locations: locations,loclatlondict=get_station_positions(locations)
    else: loclatlondict=None

    ocube=iris.load_cube(oppname,iris.AttributeConstraint(STASH='m01s00i033'))
    
    polelon=ocube.coord_system('RotatedGeogCS').grid_north_pole_longitude
    polelat=ocube.coord_system('RotatedGeogCS').grid_north_pole_latitude
    if ocube.coord_system('RotatedGeogCS') != None and polelon != 0 and polelat != 90: rotated=True

    print 'corners',(swlonin,nelonin),(swlatin,nelatin)

    if rotated: 
        (swlon,nelon),(swlat,nelat)=iris.analysis.cartography.rotate_pole(np.array([swlonin,nelonin]),np.array([swlatin,nelatin]),polelon,polelat)
        #if swlon < 0: swlon+=360. 
        #if nelon < 0: nelon+=360. 
        lon_extent=iris.coords.CoordExtent(ocube.coord('grid_longitude'),swlon,nelon)
        lat_extent=iris.coords.CoordExtent(ocube.coord('grid_latitude'),swlat,nelat)
        if locations:
            for location in locations:
                loclatlondict[location][0],loclatlondict[location][1]=iris.analysis.cartography.rotate_pole(\
                 np.array([loclatlondict[location][0]]),np.array([loclatlondict[location][1]]),polelon,polelat)
                if loclatlondict[location][0] < 0: loclatlondict[location][0]+=360.
        if drawpaths:
            for i,drawpath in enumerate(drawpaths):
                drawpaths[i][0,:],drawpaths[i][1,:]=\
                 iris.analysis.cartography.rotate_pole(drawpath[0,:],drawpath[1,:],polelon,polelat)
                ilonover=np.where(drawpaths[i][0,:] < 0)
                drawpaths[i][0,ilonover]+=360.

    print 'rotated corners',(swlon,nelon),(swlat,nelat)

    orog_zoom=ocube.intersection(lon_extent,lat_extent)

    ocontlevs=np.linspace(orogint,10000,10000/orogint)


    for velppname in velppnames:

        #analdttm=[filepart for filepart in velppname.split('/')[1:] if filepart[-1] == 'Z'][0]# and filepart[-6] == 'T'][0]
        
        print 'velppname',velppname
    
        try: ucube=iris.load_cube(velppname,iris.AttributeConstraint(STASH='m01s03i225'))
        except: ucube=iris.load(velppname,iris.AttributeConstraint(STASH='m01s03i225'))[0]
        try: vcube=iris.load_cube(velppname,iris.AttributeConstraint(STASH='m01s03i226'))
        except: vcube=iris.load(velppname,iris.AttributeConstraint(STASH='m01s03i226'))[0]
        
        uzoom=ucube.regrid(orog_zoom,iris.analysis.Linear()) # NB this is an interpolation, B grid is not the same as orog grid 
        vzoom=vcube.regrid(orog_zoom,iris.analysis.Linear())

        print 'cubes regridded to orography'
    
        if slicedttm:
            #validdttm=str(slicedttm[0])+''.join(['%02i' % dt for dt in slicedttm[1:]])#''.join(slicedttm[1:])
            slicedttm=datetime.datetime(slicedttm[0],slicedttm[1],slicedttm[2],slicedttm[3],slicedttm[4])
            print 'slicedttm',slicedttm
            with iris.FUTURE.context(cell_datetime_objects = True):
                timeconstraint=iris.Constraint(time=lambda cell: abs(cell.point-slicedttm) < datetime.timedelta(0,180)) # within 2 minutes
                uslice=uzoom.extract(timeconstraint)
                vslice=vzoom.extract(timeconstraint)
                validdttm=uslice.coord('time').units.num2date(uslice.coord('time').points[0]).strftime("%Y%m%d%H%M")
            create_single_level_plot(uslice,vslice,orog_zoom,ocontlevs,velcontmax,linestyles,drawpaths,pathcols,ocolour,coastcolour,vectcol,\
                                     polelon,polelat,swlonin,nelonin,swlatin,nelatin,xlims,ylims,validdttm,outdir,prefix,affix,locations,loclatlondict)
        else:
            # plot all slices in the file
            for i,uslice in enumerate(uzoom.slices(['grid_latitude', 'grid_longitude'])):
                vslice=vzoom[i]
                validdttm=uslice.coord('time').units.num2date(uslice.coord('time').points[0]).strftime("%Y%m%d%H%M")
                create_single_level_plot(uslice,vslice,orog_zoom,ocontlevs,velcontmax,linestyles,drawpaths,pathcols,ocolour,coastcolour,vectcol,\
                                         polelon,polelat,swlonin,nelonin,swlatin,nelatin,xlims,ylims,validdttm,outdir,prefix,affix,locations,loclatlondict)
            #else:
            #    validdttm=uzoom.coord('time').units.num2date(uzoom.coord('time').points[0]).strftime("%Y%m%d%H%M")
            #    create_single_level_plot(uzoom,vzoom,orog_zoom,ocontlevs,velcontmax,linestyles,drawpaths,pathcols,ocolour,coastcolour,vectcol,\
            #                             polelon,polelat,swlonin,nelonin,swlatin,nelatin,xlims,ylims,validdttm,outdir,prefix,affix,locations,loclatlondict)
    


def create_single_level_plot(uslice,vslice,orog_zoom,ocontlevs,velcontmax,linestyles,drawpaths,pathcols,ocolour,coastcolour,vectcol,\
                             polelon,polelat,swlonin,nelonin,swlatin,nelatin,xlims,ylims,validdttm,outdir,prefix,affix,locations,loclatlondict):
        
    speedslice=(uslice**2 + vslice**2)**0.5
    if velcontmax == None: velcontmax=( int(np.max(speedslice.data)) / 2 ) * 2
    spdcontlevs=np.linspace(0,velcontmax,1+velcontmax)        

    meanlon=(swlonin+nelonin)/2.
    meanlat=(swlatin+nelatin)/2.
    ewlength=haversine(meanlat,swlonin,meanlat,nelonin)
    nslength=haversine(swlatin,meanlon,nelatin,meanlon)
    figsize=( 7.5*ewlength/stats.mstats.gmean([ewlength,nslength]) , 7.*nslength/stats.mstats.gmean([ewlength,nslength]) )

    fig=plt.figure(figsize=figsize)
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines('10m', color = coastcolour, linewidth=1.)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linestyle = ':', draw_labels = True)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER


    theplot=iplt.contourf(speedslice,levels=spdcontlevs,cmap="jet",extend="both")
    iplt.contour(orog_zoom,levels=ocontlevs,colors=ocolour,linewidths=1.)

    xvals,yvals=orog_zoom.coord('grid_longitude').points,orog_zoom.coord('grid_latitude').points
    
#    if "km0" in prefix: ax.set_aspect(1.7*(xvals[-1]-xvals[0])/(yvals[-1]-yvals[0]))
#    else: ax.set_aspect(1.1*(xvals[-1]-xvals[0])/(yvals[-1]-yvals[0]))
#
#    if "km0" in prefix: newposn=[ax.get_position().x0,ax.get_position().y0,ax.get_position().width*0.95,ax.get_position().height]
#    else: newposn=[ax.get_position().x0*0.7,ax.get_position().y0,ax.get_position().width,ax.get_position().height]
#    ax.set_position(newposn)
    
    #ax.set_aspect(1.1*(xvals[-1]-xvals[0])/(yvals[-1]-yvals[0]))
    #xlims=ax.get_xlim()
    #ylims=ax.get_ylim()
    #print 'xlims,ylims',xlims,ylims
    ax.set_aspect( (nslength/(ylims[1]-ylims[0])) / (ewlength/(xlims[1]-xlims[0])) )
    
    axextentinches=ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    axleftedge=axextentinches.x0/fig.get_figwidth()
    axrightedge=axextentinches.x1/fig.get_figwidth()
    axbottedge=axextentinches.y0/fig.get_figheight()
    axtopedge=axextentinches.y1/fig.get_figheight()

    nvect=(15,20)
    nskipx,nskipy=max([uslice.data.shape[1]/nvect[0],1]),np.max([uslice.data.shape[0]/nvect[1],1])
    lon,lat = orog_zoom.coord('grid_longitude').points[::nskipx],orog_zoom.coord('grid_latitude').points[::nskipy]
    u,v = uslice.data[::nskipy,::nskipx],vslice.data[::nskipy,::nskipx]
    rotated_pole = ccrs.RotatedPole(pole_longitude=polelon, pole_latitude=polelat)
    #transform = orog_zoom.coord('grid_longitude').coord_system.as_cartopy_projection()
    #vectcol,spdscale='k',15
    ax.set_extent((swlonin,nelonin,swlatin,nelatin), crs=ccrs.PlateCarree())
    vectscale=nvect[0]*20.
    #arrowobj=ax.quiver(lon, lat, u, v, transform = rotated_pole, width = 0.001, headwidth = 8, headlength = 10, headaxislength = 10, scale =100, scale_units = 'xy') 
    arrowobj=ax.quiver(lon, lat, u, v, transform = rotated_pole, scale = vectscale, scale_units = 'width')#, width = 0.005, headwidth = 4, headlength = 5, \
                       #headaxislength = 10)

    plt.quiverkey(arrowobj,0.05,-0.13,10.,'10 m/s',color='k')

    if drawpaths:
        print 'adding drawn paths...'
        for i,drawpath in enumerate(drawpaths):
            if pathcols: pathcol=pathcols[i]
            else: pathcol='k'
            if linestyles: linestyle=linestyles[i]
            else: linestyle='-'
            if linestyle == '-': alpha=0.6
            elif linestyle == ':': alpha=1.
            else: alpha=0.8
            if len(drawpath[0,:]) > 2: 
                # implies a trajectory, e.g. sonde
                ax.plot(drawpath[0,:400:5],drawpath[1,:400:5],pathcol,transform = rotated_pole) # has trouble with more than 100 points in some cases
                ax.plot([drawpath[0,0]],[drawpath[1,0]],'k*',markersize=16,markeredgecolor='w',\
                                               markeredgewidth=1,transform = rotated_pole)
            else: ax.plot(drawpath[0,:],drawpath[1,:],pathcol,linestyle=linestyle,linewidth=2.5,transform = rotated_pole,alpha=alpha)


    if locations:
        print 'adding station locations....'
        for location in loclatlondict.keys():
            ax.plot(loclatlondict[location][0],loclatlondict[location][1],'ko',markeredgecolor='w',\
                    markeredgewidth=1,markersize=7,transform = rotated_pole)
            txt=ax.text(loclatlondict[location][0]+(xvals[-1]-xvals[0])/100.,\
                        loclatlondict[location][1]+(xvals[-1]-xvals[0])/100.,\
                        str(location),color='k',fontsize=13,transform = rotated_pole)
            txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='w')])

    #ax.set_xlim((xlims[0]+(xlims[1]-xlims[0])/20. , xlims[1]-(xlims[1]-xlims[0])/20.))
    #ax.set_ylim((ylims[0]+(ylims[1]-ylims[0])/20. , ylims[1]))#-(ylims[1]-ylims[0])/20.))
    if xlims: ax.set_xlim(xlims)
    if ylims: ax.set_ylim(ylims)
    
    rightspace=1.-axrightedge
    colorbar_axes = plt.gcf().add_axes([axrightedge-rightspace/8., axbottedge, rightspace/4., axtopedge-axbottedge])
    bar=plt.colorbar(theplot,colorbar_axes,orientation='vertical')
    barlabel='wind speed (ms$^{-1})$'
    bar.set_label(barlabel)

    plt.savefig(outdir+'/'+prefix+'_'+validdttm+'_'+affix+'.png')
    plt.close()



def plot_vert_xsec(velppname,thppname,oppname,slicetype,slicedttm,lonstart,latstart,lonend,latend,ylims,\
                   velcontmax,prefix,outdir,loccoords=None,resolvedir=None):
    '''
    Plot specified vertical slice through 3D model data
    Vertical or hoz velocities, or potentially other field, 
    plotted with relative distance coords and terrain outline  
    (hoz vectors using quiver not yet implemented)
    
    "plot_vert_xsec(velppname,thppname,oppname,slicetype,slicedttm,lonstart,latstart,lonend,latend,ylims,\
                    velcontmax,prefix,outdir,loccoords=None,resolvedir=None)"
    velppname,thppname,oppname: name of data files for velocities, theta and orography
    slicetype: vertvel (vertical velocity), hozcompt (horizontal wind component in-plane or 
                resolved along direction given by resolvedir) or hozspeed (horizontal wind speed)
    resolvedir: direction into which winds are resolved for plotting (default in-plane), 
                given as a vector orientation (as opposed to a normal angular wind direction definition) 
    slicedttm: list containing [year,month,day,hourminutes]
    lonstart,latstart,lonend,latend: start and end coordinates of vertical slice in unrotated lon/lat
    ylims: limits of y axis
    velcontmax: colour contour maximum
    prefix: prefix to add to image filenames
    outdir: directory to output image files
    loccoords: coordinates of a location to plot (projected) in the depicted plane
    '''
    
    npointsinspan=100 #no. points along the section 
    analdttm=[filepart for filepart in velppname.split('/')[1:] if filepart[-1] == 'Z'][0]# and filepart[-6] == 'T'][0]
    validdttm=str(slicedttm[0])+''.join(['%02i' % dt for dt in slicedttm[1:]])#''.join(slicedttm[1:])
    print validdttm
    slicedttm=datetime.datetime(slicedttm[0],slicedttm[1],slicedttm[2],slicedttm[3],slicedttm[4])
    print 'slicedttm',slicedttm

    affix='%06.2f_%06.2f_to_%06.2f_%06.2f_%05i-%05im_contmax%02i' % ( lonstart,latstart,lonend,latend,ylims[0],ylims[1],velcontmax )

    if loccoords: 
        loclatlondict={}
        for loc in loccoords:
            loclatlondict[loc[0]]=[loc[1],loc[2]]

    thcube=iris.load_cube(thppname,iris.AttributeConstraint(STASH='m01s00i004'))
    wcube=iris.load_cube(velppname,iris.AttributeConstraint(STASH='m01s00i150')) 
    ucube=iris.load_cube(velppname,iris.AttributeConstraint(STASH='m01s00i002')) 
    vcube=iris.load_cube(velppname,iris.AttributeConstraint(STASH='m01s00i003')) 
    ocube=iris.load_cube(oppname,iris.AttributeConstraint(STASH='m01s00i033'))
    ocube=ocube.regrid(thcube,iris.analysis.Linear())
    thcube,wcube,ucube,vcube=add_hybrid_height(ocube,[thcube,wcube,ucube,vcube])
    
    polelon=ocube.coord_system('RotatedGeogCS').grid_north_pole_longitude
    polelat=ocube.coord_system('RotatedGeogCS').grid_north_pole_latitude
    if ocube.coord_system('RotatedGeogCS') != None and polelon != 0 and polelat != 90: rotated=True

    print 'slice location',(lonstart,latstart),(lonend,latend)
    
    if rotated: 
        (lonstart,lonend),(latstart,latend)=iris.analysis.cartography.rotate_pole(np.array([lonstart,lonend]),np.array([latstart,latend]),polelon,polelat)
        if loccoords:
            for location in loclatlondict.keys():
                loclatlondict[location][0],loclatlondict[location][1]=iris.analysis.cartography.rotate_pole(\
                 np.array([loclatlondict[location][0]]),np.array([loclatlondict[location][1]]),polelon,polelat)
                if loclatlondict[location][0] < 0: loclatlondict[location][0]+=360.

    if lonstart < 0: lonstart+=360.
    if lonend < 0: lonend+=360.
    slicelons = np.linspace(lonstart,lonend, num=npointsinspan)
    slicelats = np.linspace(latstart,latend, num=npointsinspan)
    
    with iris.FUTURE.context(cell_datetime_objects = True):
        timeconstraint=iris.Constraint(time=lambda cell: cell.point == slicedttm)
        thcube=thcube.extract(timeconstraint)
        wcube=wcube.extract(timeconstraint)
        ucube=ucube.extract(timeconstraint)
        vcube=vcube.extract(timeconstraint)

        #thcube=thcube.extract(timeconstraint).collapsed('time',iris.analysis.MEAN)
        #ucube=ucube.extract(timeconstraint).collapsed('time',iris.analysis.MEAN)
        #vcube=vcube.extract(timeconstraint).collapsed('time',iris.analysis.MEAN)
        #wcube=wcube.extract(timeconstraint).collapsed('time',iris.analysis.MEAN)

    wslice = trajectory.interpolate(wcube, (('grid_longitude', slicelons), ('grid_latitude', slicelats)), method='nearest')
    uslice = trajectory.interpolate(ucube, (('grid_longitude', slicelons), ('grid_latitude', slicelats)), method='nearest')
    vslice = trajectory.interpolate(vcube, (('grid_longitude', slicelons), ('grid_latitude', slicelats)), method='nearest')
    thslice = trajectory.interpolate(thcube, (('grid_longitude', slicelons), ('grid_latitude', slicelats)), method='nearest') 
    orogslice = trajectory.interpolate(ocube, (('grid_longitude', slicelons), ('grid_latitude', slicelats)), method='nearest')

    print "obtained vertical data slice"
    
    # plot figure
    fig=plt.figure(figsize=(9,4))

    affix+='_'+slicetype
    if slicetype == 'hozcompt' and resolvedir: affix+='_%05.1f' % (resolvedir)
    thcontlevs=np.linspace(250,1000,376)

    if slicetype == "vertvel":
        wcontlevs=np.linspace(-velcontmax,velcontmax,4*velcontmax+1)

        plotobj=iplt.contourf(wslice,levels=wcontlevs,cmap=wcmap(),extend="both")

    elif slicetype == "hozcompt": 

        comptslice=getslicecompt(uslice,vslice,lonstart,latstart,lonend,latend,resolvedir=resolvedir)
        if velcontmax == None:
            iviewedlevs=np.where(comptslice.coord('level_height').points <= ylims[1])[0]
            velcontmax=( int(np.max(comptslice[iviewedlevs,:].data)) / 2 ) * 2
        velcontlevs=np.linspace(-velcontmax/2,velcontmax,1+3*velcontmax/2)

        plotobj=iplt.contourf(comptslice,levels=velcontlevs,cmap=velcmap(),extend="both")

    elif slicetype == "hozspeed": 

        spdslice=deepcopy(uslice)
        spdslice.data=(uslice.data**2 + vslice.data**2)**0.5
        if velcontmax == None:
            velcontmax=( int(np.max(spdslice.data)) / 2 ) * 2
        velcontlevs=np.linspace(-velcontmax/2,velcontmax,1+3*velcontmax/2)

        plotobj=iplt.contourf(spdslice,levels=velcontlevs,cmap=velcmap(),extend="both")

    else: raise Exception("unsupported slice type")

    # plot isentropes        
    iplt.contour(thslice,levels=thcontlevs,colors='k',linewidths=.75)

    # plot vectors
    extentkm=haversine(lonstart,latstart,lonend,latend)
    griddx=extentkm/(npointsinspan-1)
    xvals=np.arange(npointsinspan)
    if not(slicetype == "hozcompt"): comptslice=getslicecompt(uslice,vslice,lonstart,latstart,lonend,latend)
    #ilevs=np.where(comptslice.coord('level_height').points <= ylims[1])

    #xgrid=np.tile( xvals , [comptslice.data.shape[0],1] )[ilevs,:]
    #ygridu=comptslice.coord('altitude').points[ilevs,:]
    #ygridw=wslice.coord('altitude').points[ilevs,:]

    #nvect=[5,5]
    #nskipx=max([xgrid.shape[1]/nvect[0],1])
    #print 'nskipx',nskipx
    #print 'xgrid.shape,ygridu.shape',xgrid.shape,ygridu.shape
    #xvect,yvect = np.meshgrid(xvals[::nskipx] , np.linspace(ylims[1]/nvect[1],ylims[1],nvect[1]))
    #uvect=mlab.griddata(xgrid.flatten(),ygridu.flatten(),comptslice.data[ilevs,:].flatten(),xvect,yvect,interp='nn')
    #wvect=mlab.griddata(xgrid.flatten(),ygridw.flatten(),wslice.data[ilevs,:].flatten(),xvect,yvect,interp='nn')

    #ax=plt.gca()
    #spdscale=300.
    #arrowobj=ax.quiver(xvect, yvect, uvect, 1000*wvect,minshaft=0.,\
    #                   angles='xy',width=0.002,scale_units='width',scale=nvect[0]*spdscale,\
    #                   headwidth=10,headlength=14) 
    #plt.quiverkey(arrowobj,0.05,-0.15,20.,'20 m/s',color='k')

    # plot orography
    iplt.plot(orogslice,color='k')
    plt.fill_between(np.arange(npointsinspan),orogslice.data,color='k')#,zorder=zorder)
    plt.ylim(ylims)

    # plot locations
    if loccoords: 
        print 'adding station locations....'
        ilocs,xlocs=find_closest_point(loclatlondict,xvals,orogslice)
        for i,xloc in enumerate(xlocs):
            print 'orogslice.data[ilocs[i]]',orogslice.data[ilocs[i]]
            plt.plot(xloc,orogslice.data[ilocs[i]]+ylims[1]/80.,'k*',markeredgewidth=1,markeredgecolor='w',markersize=14)

    # plot sonde path



    xtickvals=tickspacecalc(extentkm)
    xticklocs=xvals[0]+(xtickvals/griddx)*(xvals[-1]-xvals[0])/(len(xvals)-1)
    plt.xticks(xticklocs,xtickvals)
    plt.xlabel('x (km)')
    plt.ylabel('z (m)')
    
    bar=plt.colorbar(plotobj,orientation='vertical')
    if slicetype == 'vertvel': barlabel='vertical velocity (ms$^{-1})$'
    elif slicetype == 'hozcompt': 
        if not resolvedir: barlabel='in-plane horizontal velocity (ms$^{-1})$'
        else: barlabel='resolved horizontal velocity component (ms$^{-1})$'
    elif slicetype == 'hozspeed': barlabel='horizontal wind speed (ms$^{-1})$'
    bar.set_label(barlabel)
    
    plt.tight_layout()
    
    plt.savefig(outdir+'/'+prefix+'_'+analdttm+'_'+validdttm+'_'+affix+'.png')
    plt.close()
    #plt.show()
    

