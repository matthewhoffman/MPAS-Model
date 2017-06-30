#!/usr/bin/env python
'''
This script sets up initial conditions for an experiment similar to:
Goeller, S., M. Thoma, K. Grosfeld, and H. Miller (2013), A balanced water layer concept for subglacial hydrology in large-scale ice sheet models, Cryosph., 7(4), 1095-1106, doi:10.5194/tc-7-1095-2013.

'''


import sys
from netCDF4 import Dataset
from math import sqrt
import numpy as np
from collections import Counter


# Parse options
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-f", "--file", dest="filename", type='string', help="file in which to set up the goeller initial conditions", metavar="FILE")
options, args = parser.parse_args()
if not options.filename:
   options.filename = 'landice_grid.nc'
   print 'No file specified.  Attempting to use landice_grid.nc'


# Open the file, get needed dimensions
try:
    gridfile = Dataset(options.filename,'r+')
    nCells = len(gridfile.dimensions['nCells'])
    nEdges = len(gridfile.dimensions['nEdges'])
    nVertLevels = len(gridfile.dimensions['nVertLevels'])
    nVertInterfaces = nVertLevels + 1
    maxEdges = len(gridfile.dimensions['maxEdges'])
    #WHL- Maybe do production runs with fewer than 10 levels?
    if nVertLevels != 10:
         print 'nVertLevels in the supplied file was ', nVertLevels, '.  10 levels is a preliminary value to be used with this test case.'
except:
    sys.exit('Error: The grid file specified is missing needed dimensions.')


# Put the domain origin at the center of the lower left grid cell.
# Only do this if it appears this has not already been done:

#WHL - Need [:] to actually read the data into a numpy array, and read it just once. Good for performance.
xCell = gridfile.variables['xCell'][:]
yCell = gridfile.variables['yCell'][:]

if yCell.min() > 0.0:
   print 'Shifting domain origin, because it appears that this has not yet been done.'
   unique_xs = np.array(sorted(list(set(xCell[:]))))
   unique_ys = np.array(sorted(list(set(yCell[:]))))

   print 'unique_ys.min:', unique_ys.min()
   print 'unique_ys.max:', unique_ys.max()
   print 'unique_xs.min:', unique_xs.min()
   print 'unique_xs.max:', unique_xs.max()

   xShift = -1.0 * unique_xs.min()
   yShift = -1.0 * unique_ys.min()
   gridfile.variables['xCell'][:] = xCell + xShift
   gridfile.variables['yCell'][:] = yCell + yShift
   xCell = xCell + xShift
   yCell = yCell + yShift
   gridfile.variables['xEdge'][:] = gridfile.variables['xEdge'][:] + xShift
   gridfile.variables['yEdge'][:] = gridfile.variables['yEdge'][:] + yShift
   gridfile.variables['xVertex'][:] = gridfile.variables['xVertex'][:] + xShift
   gridfile.variables['yVertex'][:] = gridfile.variables['yVertex'][:] + yShift

   # Need to adjust geometry along top and bottom boundaries to get flux correct there.
   # Essentially, we only want to model the interior half of those cells.
   # Adding this here because we only want to do this if it hasn't been done before.
   # This method is assuming a periodic_hex mesh!
   print "Adjusting areaCell and dvEdge for cells along north and south boundaries"

   # Adjust area in half for N/S boundary cells
   unique_ys = np.array(sorted(list(set(yCell[:]))))  # recalculate after above adjustment
   areaCell = gridfile.variables['areaCell']  # Note just getting object here
   areaCell[ np.nonzero(yCell == unique_ys[0]) ] *= 0.5   # cut area in half for south row
   areaCell[ np.nonzero(yCell == unique_ys[-1]) ] *= 0.5  # cut area in half for north row

   # Adjust length in half for edges connecting N/S boundary cells
   yEdge = gridfile.variables['yEdge'][:]
   dvEdge = gridfile.variables['dvEdge']  # Note just getting object here
   unique_ys_edge = np.array(sorted(list(set(yEdge[:]))))
   dvEdge[ np.nonzero(yEdge == unique_ys_edge[0]) ] *= 0.0  # zero out the edges on the boundary (not necessary because velocity will also be zero)
   dvEdge[ np.nonzero(yEdge == unique_ys_edge[1]) ] *= 0.5   # cut length in half for edges between boundary cells
   dvEdge[ np.nonzero(yEdge == unique_ys_edge[-1]) ] *= 0.0  # zero out the edges on the boundary (not necessary because velocity will also be zero)
   dvEdge[ np.nonzero(yEdge == unique_ys_edge[-2]) ] *= 0.5  # cut length in half for edges between boundary cells

#############################################################################################

# The following function computes teh bed according to Goeller et al., (2013)

def computeBedGoeller(x,y):
    Bmax = 1000.
    Bmin = 0.
    r = np.random.random_sample(np.shape(x))
    B = x*(Bmax/np.max(x))*r
    return B

# Set bedTopography (this variable should always be present in the input file)
print "Defining bedTopography"

# Compute the bed topography
bedTopography = np.zeros((nCells,))
#bedTopography = computeBedGoeller(xCell,yCell)
gridfile.variables['bedTopography'][0,:] = bedTopography[:]  # dimensions of gridfile variable are Time and nCells

# Debug: Print the topography along the bottom row.
#unique_xs = np.array(sorted(list(set(xCell[:]))))
#for iCell in range(1,nCells):
#   if yCell[iCell] == yCell.min():
#      if xCell[iCell] in unique_xs:
#         print xCell[iCell], bedTopography[iCell]



# Set the initial thickness
print "Defining thickness"

# Use a 'plastic glacier' shape as in Schoof et al. (2012) equation 5.14
# Solve the ode for shape of the surface
# for now just assume that the bed is flat

tau_d = 1.0e5
rhoi = 910.
g = 9.81
seconds_per_year = 3600.0 * 24.0 * 365.0
H0 = 0.01
dx = xCell[1]-xCell[0]
x0 = xCell[:].max()-200000.-dx     # ice sheet margin

usrf = np.zeros((nCells,))

for i in range(len(xCell[:])):
    # Calculate the constant value after integrating eqn. 5.14
    K = tau_d/(rhoi*g)*x0
    # Calculate the coefficients to solve the quadratic equation for s that remains after integrating eqn 5.14
    a = 0.5
    b = -bedTopography[i]
    c = -tau_d/(rhoi*g)*xCell[i] + K
    # Solve for s - upper surface has a sqrt shape
    if xCell[i] >= x0:
        usrf[i] = np.sqrt(b**2. - 4.0*a*c)/(2.0*a)
    else:
        usrf[i] = bedTopography[i]

thickness = usrf[:] - bedTopography[:]
gridfile.variables['thickness'][0,:] = thickness[:]



# Set the surface mass balance.
print "Defining SMB"

# Goeller et al. (2013) assumes an SMB of 0.5 m/yr.
# Convert from m/yr to kg/m2/s using appropriate ice density.

# Assign a large negative SMB along the margin where x0-dx < x < x0,  to prevent ice advancing.
# then zero mass balance for the rest of the domain to make sure that ice isn't moving past.
SMB = np.zeros((nCells,))
SMB[xCell[:] > x0-dx] = -1000. * rhoi/seconds_per_year
SMB[xCell[:]>x0] = 0.5 * rhoi/seconds_per_year
gridfile.variables['sfcMassBal'][0,:] = SMB[:]

#############################################################################################

# Approximate boundary conditions with a Dirichlet velocity mask (velocity = 0).
# Note: At the N and S boundaries, only the normal (y) velocity component
#       will be zeroed out in Albany.  The x component can be nonzero,
#       supporting a free-slip boundary condition.

if 'dirichletVelocityMask' in gridfile.variables:
   print 'dirichletVelocityMask already in gridfile'
   kinbcmask = gridfile.variables['dirichletVelocityMask']
else:
   print 'dirichletVelocityMask not in gridfile; create new variable'
   datatype = gridfile.variables['xCell'].dtype  # Get the datatype for double precision float
   dirichletVelocityMask = gridfile.createVariable('dirichletVelocityMask', datatype, ('Time','nCells','nVertInterfaces'))

print "Defining velocity boundary conditions"
kinbcmask = np.zeros((nCells, nVertInterfaces))
kinbcmask[np.nonzero(yCell == yCell.min()), : ] = 1 # south row
kinbcmask[np.nonzero(yCell == yCell.max()), : ] = 1 # north row
kinbcmask[np.nonzero(xCell > xCell.max()-dx), : ] = 1          # east boundary
gridfile.variables['dirichletVelocityMask'][0,:] = kinbcmask

# Set the initial velocities to zero to enforce Dirichlet BC..
# May not be necessary, but doing this to be on the safe side.
if 'uReconstructX' in gridfile.variables:
   print 'uReconstructX already in gridfile'
else:
   print 'uReconstructX not in gridfile; create new variable'
   datatype = gridfile.variables['xCell'].dtype  # Get the datatype for double precision float
   uReconstructX = gridfile.createVariable('uReconstructX', datatype, ('Time','nCells'))

if 'uReconstructY' in gridfile.variables:
   print 'uReconstructY already in gridfile'
else:
   print 'uReconstructY not in gridfile; create new variable'
   datatype = gridfile.variables['xCell'].dtype  # Get the datatype for double precision float
   uReconstructY = gridfile.createVariable('uReconstructY', datatype, ('Time','nCells'))

Ux_slope = (50./seconds_per_year)/(xCell[:].max()-x0)  # x-velocity, linear ramp from 50 m /yr at the margin to 0 at the divide
Ux = Ux_slope*(xCell[:]-x0)-(50./seconds_per_year)
Ux[xCell[:]<=x0] = 0.0
for i in range(nVertInterfaces):
    gridfile.variables['uReconstructX'][0,:,i] = Ux

gridfile.variables['uReconstructY'][0,:] = 0.0

#############################################################################################

# Set basal traction coefficient, beta.
# For now, assume a Weertman-type power law, tau_b = C * u^(1/m), where C = beta.
# Asay-Davis et al. (2016) specify C = 3.160 x 10^6 Pa m^{-1/3} s^{1/3} for power-law friction,
#  with friction-law exponent m = 3.
# Later, we could support a Tsai friction law.

if 'beta' in gridfile.variables:
   print 'beta already in gridfile'
   kinbcmask = gridfile.variables['beta']
else:
   print 'beta not in gridfile; create new variable'
   datatype = gridfile.variables['xCell'].dtype  # Get the datatype for double precision float
   beta = gridfile.createVariable('beta', datatype, ('Time','nCells'))

print "Defining beta"
# For the Weertman power law, beta holds the 'C' coefficient.  The beta units in MPAS are a mess right now.
# Goeller et al. (2013) steals C_0 from MISMIP3D , C = 10^7 Pa m^-1/3 s^1/3 translates to beta = 31880.

C = 1.0e7   # Pa m^{-1/3} s^{1/3}
C = C / seconds_per_year**(1.0/3.0)  # convert to MPAS units
gridfile.variables['beta'][0,:] = C


# Set up layerThicknessFractions
gridfile.variables['layerThicknessFractions'][:] = 1.0 / float(nVertLevels)

#############################################################################################

### Setup water variables

# Check for all the hydrology variables and add them if they are not currently in the input file
for var in ['waterThickness','waterPressure','basalMeltInput','externalWaterInput']:
    if var in gridfile.variables:
        print var, ' already in gridfile'
    else:
        print var, ' not in gridfile; create new variable'
        datatype = gridfile.variables['xCell'].dtype    # Get the datatype for double precision float
        gridfile.createVariable(var,datatype, ('Time','nCells'))
        gridfile.variables[var][:] = 0.0

rhow = 1000.
bmelt_slope = -(rhow/seconds_per_year)*0.08/(xCell[:].max()-x0)  # SI mass rate, linear ramp from 8 cm /yr at the margin to 0 at the divide
bmelt = bmelt_slope*(xCell[:]-x0)+(rhow/seconds_per_year)*0.08
bmelt[xCell[:]<=x0] = 0.0
gridfile.variables['basalMeltInput'][0] = bmelt

# Approximate boundary conditions with a Dirichlet velocity mask (velocity = 0).
# Note: At the N and S boundaries, only the normal (y) velocity component
#       will be zeroed out in Albany.  The x component can be nonzero,
#       supporting a free-slip boundary condition.

if 'waterFluxMask' in gridfile.variables:
   print 'waterFluxMask already in gridfile'
   waterFluxMask = gridfile.variables['waterFluxMask'][0,:]
else:
   print 'waterFluxMask not in gridfile; create new variable'
   datatype = gridfile.variables['indexToEdgeID'].dtype
   gridfile.createVariable('waterFluxMask', datatype, ('Time','nEdges'))
   waterFluxMask = gridfile.variables['waterFluxMask'][0,:]

print "Defining water flux boundary conditions"
waterFluxMask[:] = 0
cONe = gridfile.variables['cellsOnEdge'][:]
for iEdge in range(nEdges):
    for j in [0,1]:
        if cONe[iEdge,j] == 0:
            waterFluxMask[iEdge] = 2
gridfile.variables['waterFluxMask'][0,:] = waterFluxMask





#WHL - sync ensures that values are written to the file before closing.
gridfile.sync()
gridfile.close()

print 'Successfully added initial conditions to: ', options.filename
