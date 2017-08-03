#!/usr/bin/env python

"""
Plot the results for the Goeller test case.
"""

import numpy as np
import netCDF4
from optparse import OptionParser
import matplotlib.pyplot as plt
from matplotlib import cm

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
beta = f.variables['beta'][:]
BMB = f.variables['groundedBasalMassBal'][:]
days = f.variables['daysSinceStart'][:]
years = days/365.
spy = 3600.0 * 24.0 * 365.0  # seconds per year,  Note: this may be slightly wrong for some calendar types!
xtime = f.variables['xtime'][:]

print "Total number of time levels=", len(days)
print "Using time slice", time_slice, " which is year ", days[time_slice]/365.0
print "xtime=", ''.join(xtime[time_slice,:])

############################################################################################

# Adjust the x values so that 0 km is the ice sheet margin
x0 = xCell[:].max()-200000.     # ice sheet margin
if xCell[:].min() >= 0.:
    xCell[:] -= x0

# Find center row
unique_ys=np.unique(yCell[:])
unique_xs = np.unique(xCell[:])
centerY=unique_ys[np.argmin(abs(unique_ys-(unique_ys.max()-unique_ys.min())/2.))]
print "number of ys, center Y value", len(unique_ys), centerY
center_ind = np.nonzero(yCell[:] == centerY)[0]

# calculate a mean profile
def MeanProf(variable):
    prof = np.empty_like(unique_xs)
    for i in range(len(unique_xs)):
        ind = np.where(xCell[:]==unique_xs[i])
        prof[i] = np.mean(variable[ind])
    return prof

############################################################################################

print "start plotting."

# reduce the number of sampled points and plot all those that are sampled
output_spacing = 10
output_cells = np.zeros(np.shape(xCell[:]))
start = np.argmin(abs(unique_xs))
columns = range(start,len(unique_xs),output_spacing)
for col in columns:
    ind = np.where(xCell[:]==unique_xs[col])
    output_cells[ind] = 1.

#############################################################################################

### plot variables at the indexed locations to be sure of convergence onto a steady state

fig = plt.figure(2, figsize = (9,6), facecolor='w')
fig.suptitle('Velocity Variables')
nplot = 3
colors = [ cm.jet(x) for x in np.linspace(0.0, 1.0, len(unique_xs)) ]

### Velocity Variables

# thickness
ax1 = fig.add_subplot(nplot,1,1)
for col in columns:
    ind = np.where(xCell[:]==unique_xs[col])[0]
    plt.plot(years, np.mean(H[:,ind],axis=1), '-',color=colors[col],label=unique_xs[col]/1000.)
plt.xlabel('Time (years)')
plt.ylabel('Ice Thickness (m)')
plt.legend(fontsize=6)
plt.grid(True)

# surface speed
ax = fig.add_subplot(nplot,1,2, sharex=ax1)
for col in columns:
    ind = np.where(xCell[:]==unique_xs[col])[0]
    plt.plot(years, np.mean(spy*Usfc[:,ind],axis=1), '-',color=colors[col],label=unique_xs[col]/1000.)
plt.xlabel('Time (years)')
plt.ylabel('Surface Speed (m a-1)')
plt.legend(fontsize=6)
plt.grid(True)

# basal speed
ax = fig.add_subplot(nplot,1,3, sharex=ax1)
for col in columns:
    ind = np.where(xCell[:]==unique_xs[col])[0]
    plt.plot(years, np.mean(spy*Ubsl[:,ind],axis=1), '-',color=colors[col],label=unique_xs[col]/1000.)
plt.xlabel('Time (years)')
plt.ylabel('Basal Speed (m a-1)')
plt.legend(fontsize=6)
plt.grid(True)

### Hydrology Variables

fig = plt.figure(3, figsize = (9,6), facecolor='w')
fig.suptitle('Hydrology Variables')

colors = [ cm.jet(x) for x in np.linspace(0.0, 1.0, len(unique_xs)) ]
nplot = 2

# water thickness
ax1 = fig.add_subplot(nplot,1,1)
for col in columns:
    ind = np.where(xCell[:]==unique_xs[col])[0]
    plt.plot(years, np.mean(h[:,ind],axis=1), '-',color=colors[col],label=unique_xs[col]/1000.)
plt.xlabel('Time (years)')
plt.ylabel('Water Thickness (m)')
plt.legend(fontsize=6)
plt.grid(True)

# Effective Pressure
ax = fig.add_subplot(nplot,1,2, sharex=ax1)
for col in columns:
    ind = np.where(xCell[:]==unique_xs[col])[0]
    plt.plot(years, np.mean(N[:,ind],axis=1), '-',color=colors[col],label=unique_xs[col]/1000.)
plt.xlabel('Time (years)')
plt.ylabel('Effective Pressure (Pa)')
plt.legend(fontsize=6)
plt.grid(True)

#############################################################################################

### plot variables in 2-D from the last time step over the entire domain

# thickness
fig = plt.figure(4, figsize=(9,6), facecolor='w', dpi=100)

fig.add_subplot(111)
plt.scatter(xCell[:]/1000.0,yCell[:]/1000.0,markersize,H[time_slice,:], marker='h', edgecolors='none')
plt.colorbar()
plt.axis('equal')
plt.title('Modeled thickness (m) \n at time ' + netCDF4.chartostring(xtime)[time_slice].strip() )
plt.xlabel('x (km)'); plt.ylabel('y (km)')

# surface speed
fig = plt.figure(5, figsize=(9,6), facecolor='w', dpi=100)

fig.add_subplot(111)
plt.scatter(xCell[:]/1000.0,yCell[:]/1000.0,markersize,Usfc[time_slice,:]*spy, marker='h', edgecolors='none')
plt.colorbar()
plt.axis('equal')
plt.title('Modeled Surface Speed (m a-1) \n at time ' + netCDF4.chartostring(xtime)[time_slice].strip() )
plt.xlabel('x (km)'); plt.ylabel('y (km)')

# Effective Pressure
fig = plt.figure(6, figsize=(9,6), facecolor='w', dpi=100)

fig.add_subplot(111)
plt.scatter(xCell[:]/1000.0,yCell[:]/1000.0,markersize,N[time_slice,:], marker='h', edgecolors='none')
plt.colorbar()
plt.axis('equal')
plt.title('Effective Pressure (Pa) \n at time ' + netCDF4.chartostring(xtime)[time_slice].strip() )
plt.xlabel('x (km)'); plt.ylabel('y (km)')

# Water Thickness
fig = plt.figure(7, figsize=(9,6), facecolor='w', dpi=100)

fig.add_subplot(111)
plt.scatter(xCell[:]/1000.0,yCell[:]/1000.0,markersize,h[time_slice,:], marker='h', edgecolors='none')
plt.colorbar()
plt.axis('equal')
plt.title('Water Thickness (m) \n at time ' + netCDF4.chartostring(xtime)[time_slice].strip() )
plt.xlabel('x (km)'); plt.ylabel('y (km)')



########################################################################################################################

### plot 1-D profiles at initial and final time steps

fig = plt.figure(8, figsize = (9,6), facecolor='w')
nplot=3

### Velocity Variables

# Glacier profile
ax1 = fig.add_subplot(nplot,1,1)
plt.plot(unique_xs/1000.,MeanProf(beta[0]))
plt.xlabel('X-distance (km)')
plt.ylabel('Ice Thickness (m)')
plt.legend(fontsize=6)
plt.grid(True)

# Surface Speed
ax1 = fig.add_subplot(nplot,1,2)
plt.plot(unique_xs/1000.,spy*MeanProf(Usfc[0]))
plt.plot(unique_xs/1000.,spy*MeanProf(Usfc[-1]))
plt.xlabel('X-distance (km)')
plt.ylabel('Surface Speed (m a-1)')
plt.legend(fontsize=6)
plt.grid(True)

# Basal Mass Balance
ax1 = fig.add_subplot(nplot,1,3)
plt.plot(unique_xs/1000.,(spy/rhow)*MeanProf(BMB[0]))
plt.plot(unique_xs/1000.,(spy/rhow)*MeanProf(BMB[-1]))
plt.xlabel('X-distance (km)')
plt.ylabel('Basal Mass Balance (m a-1)')
plt.legend(fontsize=6)
plt.grid(True)

### Hydrology Variables

fig = plt.figure(9, figsize = (9,6), facecolor='w')
nplot=2

# Water Thickness
ax1 = fig.add_subplot(nplot,1,1)
plt.plot(unique_xs/1000.,MeanProf(h[0]),label='initial')
plt.plot(unique_xs/1000.,MeanProf(h[-1]),label='final')
plt.xlabel('X-distance (km)')
plt.ylabel('Water Thickness (m)')
plt.grid(True)

# Effective Pressure
ax1 = fig.add_subplot(nplot,1,2)
plt.plot(unique_xs/1000.,MeanProf(N[0]),label='initial')
plt.plot(unique_xs/1000.,MeanProf(N[-1]),label='final')
plt.xlabel('X-distance (km)')
plt.ylabel('Effective Pressure (Pa)')
plt.legend(loc=4)
plt.grid(True)


print "plotting complete"
plt.draw()
if options.saveimages:
        print "Saving figures to files."
        plt.savefig('GL-position.png')




if options.hidefigs:
     print "Plot display disabled with -n argument."
else:
     plt.show()

