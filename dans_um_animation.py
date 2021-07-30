#!/usr/local/sci/bin/python2.7

# simple python program to read in some UM data and plot a map

# import modules which are needed
import iris #all iris routines used for manipulating UM data
#import iris.quickplot as qplt #specific iris routing for plotting data
import matplotlib as mpl #all matplotlib routines used for plotting
import matplotlib.pyplot as plt #general matplotlib plotting routine
import matplotlib.cm as mpl_cm
import numpy as np
import matplotlib.animation as animation
import datetime as dt
from matplotlib import dates

# Set up formatting for the movie files -  may need to adjust this - e.g. more or less frames per second.
Writer = animation.writers['ffmpeg']
writer = Writer(fps=1, metadata=dict(artist='Me'), bitrate=-1)

# select field to plot from stashcode given in stash.txt file
# format is always m01sXXiYYY where XX is the section number
# and YYY is the item number
stash = 'm01s00i254' # qcl

# location of the data files
datapath = '/projects/phdcase/dansm/'
#IOP - case study
iop = 'iop1/'
# identifier of the UM run
filestart = 'lmqva' #London Model-like model

# data stream containing the stashcode - 60=a, 61=b etc
# again, given in the stash.txt file
fileid = 'g' #69

# load the field into python
# if you subsequently do "print field" you will get some info
# about what is in the cube
field = iris.load_cube(datapath+iop+filestart+'a_p'+fileid+'*',iris.AttributeConstraint(STASH=stash))

# select field to plot from stashcode given in stash.txt file
# format is always m01sXXiYYY where XX is the section number
# and YYY is the item number

stash_2 = 'm01s00i033' #Orography
fileid_2 ='g'
orog = iris.load_cube(datapath+iop+filestart+'a_p'+fileid_2+'*',iris.AttributeConstraint(STASH=stash_2))
fileid = 'j'

stash = 'm01s03i225' # U component of 10 m wind
stash_2 = 'm01s03i226' # V component of 10 m wind

u = iris.load_cube(datapath+iop+filestart+'a_p'+fileid+'*',iris.AttributeConstraint(STASH=stash))
v = iris.load_cube(datapath+iop+filestart+'a_p'+fileid+'*',iris.AttributeConstraint(STASH=stash))

x = field.coord('grid_longitude').points
y = field.coord('grid_latitude').points
u_lon = u.coord('grid_longitude').points
u_lat = u.coord('grid_latitude').points
start_date = dt.datetime(2014, 11, 24, 22, 0, 0) #equal to forecast_reference time
z = dates.date2num(start_date)
plot_time_list = z + field.coord('forecast_period').points/24
x_100 = field.coord('grid_longitude').points
y_100 = field.coord('grid_latitude').points

field = field[48:-1,:,:,:].data
orog = orog[48:-1,:,:,].data
u = u[48:-1,:,:,:].data
v = v[48:-1,:,:,:].data

# load a colourmap for the data, others can be found here:
# http://matplotlib.org/examples/color/colormaps_reference.html
my_cmap = plt.cm.gnuplot_r
my_cmap.set_under('w')
# normalise the colourmap to be between a certain min and max
my_norm=mpl.colors.Normalize(vmin=0.01,vmax=0.49, clip=True)
# Define levels for contourf
levels = np.linspace(0.01,0.49,25)
orog_levels = [0,100,200,300,400,500]
# plot that data, rasterising it significantly reduces the plotting time
# and size of plots
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
ax.set_aspect('equal',adjustable='box') #ensures the correct plot aspect ratio
# Plot qcl at model level 1, convert to g/kg. Levels [0] is only to ensure set_under works.
a = ax.contourf(x, y, field[0,0,:,:]*1000,  levels= [0]+levels,  cmap=my_cmap, norm=my_norm, extend="max", rasterized=True)
#Add colorbar
cbar = plt.colorbar( a, ax=ax, orientation='horizontal')
cbar.ax.set_xlabel('5m LWC (gkg$^{-1}$)')
#Define the animation function - this is called by the .FuncAnimation command. 
def animate(i):
       ax.clear()
       a = ax.contourf(x, y, field[i,0,:,:]*1000,  levels,  cmap=my_cmap, norm=my_norm, extend="max", rasterized=True)
       CS = plt.contour(x,y,orog[i,:,:],orog_levels, colors='k', rasterized=True)
       plt.clabel(CS, inline=1,fmt='%.0f', fontsize=10)
       Q = plt.quiver(u_lon[::5], u_lat[::5], u[i, ::5, ::5], v[i, ::5, ::5], pivot='middle', color='black')
       qk = plt.quiverkey(Q, 0.15, 0.15, 2, r'$2 \frac{m}{s}$', labelpos='E', coordinates='figure')
       ax.xaxis.set_major_formatter(plt.NullFormatter())
       ax.yaxis.set_major_formatter(plt.NullFormatter())
       ax.set_title(dates.num2date(plot_time_list[i]).strftime('%H:%M'))
       return ax
#Create animation of figure for 34 frames. May need to change inveral etc until happy. i in animate is given as frames (I think).
ani = animation.FuncAnimation(fig,animate,34,interval=1*1e+3,blit=False)

ani.save('5mqcl_iop1_lm_gnuplot_300dpi.mp4', writer=writer, dpi=300)

