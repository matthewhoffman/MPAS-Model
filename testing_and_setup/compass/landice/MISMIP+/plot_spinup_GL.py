#!/usr/bin/env python

'''
Plot the total grounded area of spinup cases for MISMIP+
Separate figure for the sensitivity test on each variable
'''

import numpy as np
import pylab as plt
from netCDF4 import Dataset as NetCDFFile
import os
from matplotlib import cm

f_netcdf_in = 'globalStats.nc'

A_values = [1e-05,2e-05]
C0_values = [0.25,0.5,1.,2.]
cond_values = [1e-5,2.5e-5,5e-5,7.5e-5,0.0001,0.001]
rough_values = [0.05,0.125,0.25,1.0,2.0,5.0]
bump_values = [0.01,0.05,0.5,1.0]

############################################################################################################

### variable A with no hydrology

fig = plt.figure(1,figsize=(9,4.5))

plt.xlim(0,20000)
plt.xticks(np.arange(0,20000,2000))
plt.ylim(30,45)
plt.xlabel('Time (yrs)')
plt.ylabel(r'Grounded Area (1000 km$^2$)')
plt.title('Control (No Hydrology)',size=16)
# This is *approximately* the region of the retrograde slope
plt.fill_between(np.linspace(0,20000,10),np.ones(10)*36,np.ones(10)*41.5,color='k',alpha=0.1)

for A_value in A_values:
    # change line style based on X_al
    if A_value == 1e-5:
        ls = '-'
    else:
        ls = ':'

    # import data for the spinup run and plot
    directory = './A_' + str(A_value)+'_nohydrology'

    if os.path.exists(f_netcdf_in):
        print directory, ' exists, loading...'
        os.chdir(directory)
        f = NetCDFFile(f_netcdf_in,'r+')
        gArea = f.variables['groundedIceArea'][:]/(1e6*1e3)
        plt.plot(100.*np.arange(len(gArea)),gArea,ls=ls,c='k',label=A_value)
        os.chdir('../')
    else:
        print directory, ' does not exist.'


plt.legend(title='rate factor A')

############################################################################################################

### variable C0

fig = plt.figure(2,figsize=(9,4.5))

plt.xlim(0,20000)
plt.xticks(np.arange(0,20000,2000))
plt.ylim(30,45)
plt.xlabel('Time (yrs)')
plt.ylabel(r'Grounded Area (1000 km$^2$)')
plt.title('Variable C0',size=16)
# This is *approximately* the region of the retrograde slope
plt.fill_between(np.linspace(0,20000,10),np.ones(10)*36,np.ones(10)*41.5,color='k',alpha=0.1)

cond_value = 0.001
nl_alias='cond'
A_value = 1e-5

colors = [ cm.cool(x) for x in np.linspace(0.0, 1.0, 2*len(C0_values)) ]
mvcolor = int(.5*len(C0_values))

for C0_value in C0_values:
    # change line color based on C0_value
    c = colors[mvcolor+np.argwhere(np.array(C0_values)==C0_value)[0][0]]

    # import data for the spinup run and plot
    directory = './A_' + str(A_value)+'_'+nl_alias+'_'+ str(cond_value) + '_C0_' + str(C0_value)
    if os.path.exists(directory):
        print directory, ' exists, loading...'
        os.chdir(directory)
        if os.path.exists(f_netcdf_in):
            f = NetCDFFile(f_netcdf_in,'r+')
            gArea = f.variables['groundedIceArea'][:]/(1e6*1e3)
            plt.plot(100.*np.arange(len(gArea)),gArea,c=c,label=C0_value)
        else:
            print f_netcdf_in, 'does not exist.'
        os.chdir('../')
    else:
        print directory, ' does not exist.'

plt.legend(title='C0')

##########################################################################################

### variable conductivity

fig = plt.figure(3,figsize=(9,4.5))

plt.xlim(0,20000)
plt.xticks(np.arange(0,20000,2000))
plt.ylim(30,45)
plt.xlabel('Time (yrs)')
plt.ylabel(r'Grounded Area (1000 km$^2$)')
plt.title('Variable Conductivity',size=16)
# This is *approximately* the region of the retrograde slope
plt.fill_between(np.linspace(0,20000,10),np.ones(10)*36,np.ones(10)*41.5,color='k',alpha=0.1)

C0_value = 1.0
nl_alias='cond'
A_value = 1e-5

colors = [ cm.CMRmap(x) for x in np.linspace(0.0, 1.0, 2*len(cond_values)) ]
mvcolor = int(.5*len(cond_values))

for cond_value in cond_values:
    # change line color based on cond_value
    c = colors[mvcolor+np.argwhere(np.array(cond_values)==cond_value)[0][0]]

    # import data for the spinup run and plot
    directory = './A_' + str(A_value)+'_'+nl_alias+'_'+ str(cond_value) + '_C0_' + str(C0_value)
    if os.path.exists(directory):
        print directory, ' exists, loading...'
        os.chdir(directory)
        if os.path.exists(f_netcdf_in):
            f = NetCDFFile(f_netcdf_in,'r+')
            gArea = f.variables['groundedIceArea'][:]/(1e6*1e3)
            plt.plot(100.*np.arange(len(gArea)),gArea,c=c,label=cond_value)
        else:
            print f_netcdf_in, 'does not exist.'
        os.chdir('../')
    else:
        print directory, ' does not exist.'

plt.legend(title='k')

plt.draw()
plt.savefig('Spinup_Conductivity_TwoLines.png',dpi=300)

##########################################################################################

### variable roughness

fig = plt.figure(4,figsize=(9,4.5))

plt.xlim(0,20000)
plt.xticks(np.arange(0,20000,2000))
plt.ylim(30,45)
plt.xlabel('Time (yrs)')
plt.ylabel(r'Grounded Area (1000 km$^2$)')
plt.title('Variable Roughness',size=16)
# This is *approximately* the region of the retrograde slope
plt.fill_between(np.linspace(0,20000,10),np.ones(10)*36,np.ones(10)*41.5,color='k',alpha=0.1)

C0_value = 1.0
nl_alias='rough'
A_value = 1e-5

colors = [ cm.hot(x) for x in np.linspace(0.0, 1.0, 2*len(rough_values)) ]
mvcolor = int(.5*len(rough_values))

for rough_value in rough_values:
    # change line color based on X_nl
    c = colors[mvcolor+np.argwhere(np.array(rough_values)==rough_value)[0][0]]

    # import data for the spinup run and plot
    directory = './A_' + str(A_value)+'_'+nl_alias+'_'+ str(rough_value) + '_C0_' + str(C0_value)

    if os.path.exists(directory):
        print directory, ' exists, loading...'
        os.chdir(directory)
        if os.path.exists(f_netcdf_in):
            f = NetCDFFile(f_netcdf_in,'r+')
            gArea = f.variables['groundedIceArea'][:]/(1e6*1e3)
            plt.plot(100.*np.arange(len(gArea)),gArea,c=c,label=rough_value)
        else:
            print f_netcdf_in, 'does not exist.'
        os.chdir('../')
    else:
        print directory, ' does not exist.'

plt.legend(title='roughness')

############################################################################################################

### variable bump height

fig = plt.figure(5,figsize=(9,4.5))

plt.xlim(0,20000)
plt.xticks(np.arange(0,20000,2000))
plt.ylim(30,45)
plt.xlabel('Time (yrs)')
plt.ylabel(r'Grounded Area (1000 km$^2$)')
plt.title('Variable Bump Height',size=16)
# This is *approximately* the region of the retrograde slope
plt.fill_between(np.linspace(0,20000,10),np.ones(10)*36,np.ones(10)*41.5,color='k',alpha=0.1)

C0_value = 1.0
nl_alias='bump'
A_value = 1e-5

colors = [ cm.hot(x) for x in np.linspace(0.0, 1.0, 2*len(bump_values)) ]
mvcolor = int(.5*len(bump_values))

for bump_value in bump_values:
    # change line color based on X_nl
    c = colors[mvcolor+np.argwhere(np.array(bump_values)==bump_value)[0][0]]

    # import data for the control run and plot
    directory = './A_' + str(A_value)+'_'+nl_alias+'_'+ str(bump_value) + '_C0_' + str(C0_value)

    if os.path.exists(directory):
        print directory, ' exists, loading...'
        os.chdir(directory)
        if os.path.exists(f_netcdf_in):
            f = NetCDFFile(f_netcdf_in,'r+')
            gArea = f.variables['groundedIceArea'][:]/(1e6*1e3)
            plt.plot(100.*np.arange(len(gArea)),gArea,c=c,label=nl_value)
        else:
            print f_netcdf_in, 'does not exist.'
        os.chdir('../')
    else:
        print directory, ' does not exist.'

plt.legend(title='bump height')

############################################################################################################

### No Retrograde Bedslope

fig = plt.figure(6,figsize=(9,4.5))

plt.xlim(0,20000)
plt.xticks(np.arange(0,20000,2000))
plt.ylim(30,45)
plt.xlabel('Time (yrs)')
plt.ylabel(r'Grounded Area (1000 km$^2$)')
plt.title('No Retrograde Bedslope',size=16)
# This is *approximately* the region of the retrograde slope
plt.fill_between(np.linspace(0,20000,10),np.ones(10)*36,np.ones(10)*41.5,color='k',alpha=0.1)

# import data for the control run and plot
directory = './NoRetrograde/A_1e-05_conduc_1e-05_C0_1.0'

if os.path.exists(directory):
    os.chdir(directory)
    if os.path.exists(f_netcdf_in):
        print directory, ' exists, loading...'
        f = NetCDFFile(f_netcdf_in,'r+')
        gArea = f.variables['groundedIceArea'][:]/(1e6*1e3)
        plt.plot(100.*np.arange(len(gArea)),gArea,c='k')
    os.chdir('../../')
else:
    print directory, ' does not exist.'


plt.draw()
plt.show()

