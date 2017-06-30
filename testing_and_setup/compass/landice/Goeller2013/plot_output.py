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
        options.filename = "output.nc"

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
H = f.variables['thickness'][:]
Usfc = f.variables['uReconstructX'][:,:,0]
Ubsl = f.variables['uReconstructX'][:,:,-1]
h = f.variables['waterThickness'][:]
N = f.variables['effectivePressure'][:]
days = f.variables['daysSinceStart'][:]
years = days/365.
spy = 3600.0 * 24.0 * 365.0  # seconds per year,  Note: this may be slightly wrong for some calendar types!
xtime = f.variables['xtime'][:]

print "Total number of time levels=", len(days)
print "Using time slice", time_slice, " which is year ", days[time_slice]/365.0
print "xtime=", ''.join(xtime[time_slice,:])

# Find center row
unique_ys=np.unique(yCell[:])
centerY=unique_ys[np.argmin(abs(unique_ys-(unique_ys.max()-unique_ys.min())/2.))]
print "number of ys, center Y value", len(unique_ys), centerY
ind = np.nonzero(yCell[:] == centerY)[0]

############################################################################################

print "start plotting."

# reduce the number of sampled points and plot all those that are sampled
ind = ind[::5]
output_cells = np.zeros(np.shape(xCell[:]))
output_cells[ind] = 1.

### plot the indexed values

markersize = 30.0
gray = np.ones(3)*0.8
maskindices = np.nonzero(H[time_slice,:] >= 0.0)

fig = plt.figure(1, figsize=(9,6), facecolor='w', dpi=100)
fig.add_subplot(1,1,1)

plt.scatter(xCell[maskindices]/1000.0,yCell[maskindices]/1000.0,markersize,output_cells[maskindices], marker='h', edgecolors='none')
plt.colorbar()
plt.axis('equal')
plt.title('Output Cell Locations')
plt.xlabel('x (km)'); plt.ylabel('y (km)')

#############################################################################################

### plot variables at the indexed locations to be sure of convergence onto a steady state

fig = plt.figure(2, figsize = (9,6), facecolor='w')

fig.suptitle('Velocity Variables')

nplot = 3
# thickness
ax1 = fig.add_subplot(nplot,1,1)
for i in ind:
    plt.plot(years, H[:,i], '-',label=xCell[i]/1000.)
plt.xlabel('Time (years)')
plt.ylabel('Ice Thickness (m)')
plt.legend(fontsize=6)
plt.grid(True)

# surface speed
ax = fig.add_subplot(nplot,1,2, sharex=ax1)
for i in ind:
    plt.plot(years, Usfc[:,i]*spy, '-',label=xCell[i]/1000.)
plt.xlabel('Time (years)')
plt.ylabel('Surface Speed (m a-1)')
plt.legend(fontsize=6)
plt.grid(True)

# basal speed
ax = fig.add_subplot(nplot,1,3, sharex=ax1)
for i in ind:
    plt.plot(years, Ubsl[:,i]*spy, '-',label=xCell[i]/1000.)
plt.xlabel('Time (years)')
plt.ylabel('Basal Speed (m a-1)')
plt.legend(fontsize=6)
plt.grid(True)




fig = plt.figure(3, figsize = (9,6), facecolor='w')

fig.suptitle('Hydrology Variables')

nplot = 2
# water thickness
ax1 = fig.add_subplot(nplot,1,1)
for i in ind:
    plt.plot(years, h[:,i], '-',label=xCell[i]/1000.)
plt.xlabel('Time (years)')
plt.ylabel('water Thickness (m)')
plt.legend(fontsize=6)
plt.grid(True)

# Effective Pressure
ax = fig.add_subplot(nplot,1,2, sharex=ax1)
for i in ind:
    plt.plot(years, N[:,i], '-',label=xCell[i]/1000.)
plt.xlabel('Time (years)')
plt.ylabel('Effective Pressure (Pa)')
plt.legend(fontsize=6)
plt.grid(True)

#############################################################################################

### plot variables from the last time step over the entire domain

# thickness
fig = plt.figure(4, figsize=(9,6), facecolor='w', dpi=100)

fig.add_subplot(111)
plt.scatter(xCell[maskindices]/1000.0,yCell[maskindices]/1000.0,markersize,H[time_slice,maskindices], marker='h', edgecolors='none')
plt.colorbar()
plt.axis('equal')
plt.title('Modeled thickness (m) \n at time ' + netCDF4.chartostring(xtime)[time_slice].strip() )
plt.xlabel('x (km)'); plt.ylabel('y (km)')

# surface speed
fig = plt.figure(5, figsize=(9,6), facecolor='w', dpi=100)

fig.add_subplot(111)
plt.scatter(xCell[maskindices]/1000.0,yCell[maskindices]/1000.0,markersize,Usfc[time_slice,maskindices]*spy, marker='h', edgecolors='none')
plt.colorbar()
plt.axis('equal')
plt.title('Modeled Surface Speed (m a-1) \n at time ' + netCDF4.chartostring(xtime)[time_slice].strip() )
plt.xlabel('x (km)'); plt.ylabel('y (km)')

# Effective Pressure
fig = plt.figure(6, figsize=(9,6), facecolor='w', dpi=100)

fig.add_subplot(111)
plt.scatter(xCell[maskindices]/1000.0,yCell[maskindices]/1000.0,markersize,N[time_slice,maskindices], marker='h', edgecolors='none')
plt.colorbar()
plt.axis('equal')
plt.title('Effective Pressure (Pa) \n at time ' + netCDF4.chartostring(xtime)[time_slice].strip() )
plt.xlabel('x (km)'); plt.ylabel('y (km)')




print "plotting complete"
plt.draw()
if options.saveimages:
        print "Saving figures to files."
        plt.savefig('GL-position.png')




if options.hidefigs:
     print "Plot display disabled with -n argument."
else:
     plt.show()

