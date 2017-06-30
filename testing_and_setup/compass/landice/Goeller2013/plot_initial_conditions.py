#!/usr/bin/env python
'''
Plot the results for the model spinup.
'''
import numpy as np
import netCDF4
#import datetime
# import math
# from pylab import *
from optparse import OptionParser
import matplotlib.pyplot as plt
from matplotlib import cm
# from matplotlib.contour import QuadContourSet
# import time

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", help="file to visualize", metavar="FILE")
parser.add_option("-t", "--time", dest="time", help="time step to visualize (0 based)", metavar="TIME")
parser.add_option("-s", "--save", action="store_true", dest="saveimages", help="include this flag to save plots as files")
parser.add_option("-n", "--nodisp", action="store_true", dest="hidefigs", help="include this flag to not display plots (usually used with -s)")

options, args = parser.parse_args()

if not options.filename:
	print "No filename provided. Using output.nc."
        options.filename = "landice_grid.nc"

if not options.time:
	print "No time provided. Using time -1."
        time_slice = -1
else:
        time_slice = int(options.time)

#############################################################################################

# Import desired variables from the input file

f = netCDF4.Dataset(options.filename,'r')
xCell = f.variables['xCell'][:]
yCell = f.variables['yCell'][:]
xEdge = f.variables['xEdge'][:]
yEdge = f.variables['yEdge'][:]
H = f.variables['thickness'][0,:]
b = f.variables['bedTopography'][0,:]
Ubsl = f.variables['uReconstructX'][0,:,-1]
SMB = f.variables['sfcMassBal'][0,:]

# Find center row
unique_ys=np.unique(yCell[:])
centerY=unique_ys[np.argmin(abs(unique_ys-(unique_ys.max()-unique_ys.min())/2.))]
print "number of ys, center Y value", len(unique_ys), centerY
ind = np.nonzero(yCell[:] == centerY)[0]

############################################################################################

print "start plotting."

# reduce the number of sampled points and plot all those that are sampled
output_cells = np.zeros(np.shape(xCell[:]))
output_cells[ind] = 1.

### plot the indexed values

markersize = 30.0
gray = np.ones(3)*0.8

fig = plt.figure(1, figsize=(9,6), facecolor='w', dpi=100)
fig.add_subplot(1,1,1)

plt.scatter(xCell[:]/1000.0,yCell[:]/1000.0,markersize,output_cells[:], marker='h', edgecolors='none')
plt.colorbar()
plt.axis('equal')
plt.title('Output Cell Locations')
plt.xlabel('x (km)'); plt.ylabel('y (km)')

#############################################################################################

### plot variables at the indexed locations to be sure of convergence onto a steady state

fig = plt.figure(2, figsize = (9,6), facecolor='w')
nplot=3

# Initial glacier profile
ax1 = fig.add_subplot(nplot,1,1)
plt.plot(xCell[ind]/1000.,H[ind])
plt.plot(xCell[ind]/1000.,b[ind])
plt.xlabel('X-distance (km)')
plt.ylabel('Ice Thickness (m)')
plt.legend(fontsize=6)
plt.grid(True)

# Surface Mass Balance
ax1 = fig.add_subplot(nplot,1,2)
plt.plot(xCell[ind]/1000.,SMB[ind])
plt.xlabel('X-distance (km)')
plt.ylabel('Surface Mass Balance (kg m2 s-1)')
plt.legend(fontsize=6)
plt.grid(True)

# Sliding Speed
ax1 = fig.add_subplot(nplot,1,3)
plt.plot(xCell[ind]/1000.,Ubsl[ind])
plt.xlabel('X-distance (km)')
plt.ylabel('Sliding Speed (m s-1)')
plt.legend(fontsize=6)
plt.grid(True)

#############################################################################################


print "plotting complete"
plt.draw()
if options.saveimages:
        print "Saving figures to files."
        plt.savefig('GL-position.png')




if options.hidefigs:
     print "Plot display disabled with -n argument."
else:
     plt.show()

