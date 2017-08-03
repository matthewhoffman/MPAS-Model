#!/usr/bin/env python
'''
Plot the initial conditions for the Goeller test case.
'''
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
bmelt = f.variables['basalMeltInput'][0,:]

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

### plot variables at the indexed locations to be sure of convergence onto a steady state

fig = plt.figure(2, figsize = (12,9), facecolor='w')
nplot=4

# Initial glacier profile
ax1 = fig.add_subplot(nplot,1,1)
plt.plot(unique_xs/1000.,MeanProf(H[:]))
plt.plot(unique_xs/1000.,MeanProf(b[:]))
plt.xlabel('X-distance (km)')
plt.ylabel('Ice Thickness (m)')
plt.legend(fontsize=6)
plt.grid(True)

# Surface Mass Balance
ax1 = fig.add_subplot(nplot,1,2)
plt.plot(unique_xs/1000.,MeanProf(SMB[:]))
plt.xlabel('X-distance (km)')
plt.ylabel('Surface Mass Balance (kg m2 s-1)')
plt.legend(fontsize=6)
plt.grid(True)

# Basal Melt Input
ax1 = fig.add_subplot(nplot,1,3)
plt.plot(unique_xs/1000.,MeanProf(bmelt[:]))
plt.xlabel('X-distance (km)')
plt.ylabel('Basal Melt Input (kg m2 s-1)')
plt.legend(fontsize=6)
plt.grid(True)

# Sliding Speed
ax1 = fig.add_subplot(nplot,1,4)
plt.plot(unique_xs/1000.,MeanProf(Ubsl[:]))
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

