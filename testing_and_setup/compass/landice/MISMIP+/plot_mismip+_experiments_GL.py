#!/usr/bin/env python

'''
Use this script to plot the grounded area through time from mismip+ experiments
'''

import numpy as np
import pylab as plt
from netCDF4 import Dataset as NetCDFFile
import os

stats_output = 'globalStats.nc'
ncf_output = 'output_00000.nc'

A_values = [1e-05]
cond_values = [0.0001,0.0005,0.001]
C0_values = [0.5,1.0,2.0]

for A_value in A_values:
    for cond_value in cond_values:
        for C0_value in C0_values:
            directory = 'A_' + str(A_value) + '_cond_' str(cond_value) + '_C0_' + str(C0_value)
            if os.path.exists(directory):
                print 'Directory, ', directory, ' exists, moving there...'
                os.chdir(directory)
                fig = plt.figure(figsize=(12,6.5))
            else:
                print 'Directory, ', directory, ' does not exist'
                continue

            # plot grounded area for both an evolving hydrology case and a constant basal traction case
            for trial in ['constant_beta','evolve']:
                if os.path.exists(trial):
                    print 'Directory, ', trial, ' exists, loading...'
                    os.chdir(trial)
                else:
                    print 'Directory, ', trial, ' does not exist'
                    continue

                # use different symbols for the evolving/constant cases
                if trial == 'evolve':
                    m = 'o'
                if trial == 'constant_beta':
                    m = 's'

                # plot lines for every experiment, see Asay-Davis et al. (2016) for details
                for exp in ['Ice0','Ice1r','Ice1ra','Ice1rr','Ice1rax','Ice1rrx',
                        'Ice2r','Ice2ra','Ice2rr','Ice2rax','Ice2rrx']:
                    print 'Move to directory ', exp
                    os.chdir(exp)

                    if exp == 'Ice0':
                        c = 'grey'
                    elif exp == 'Ice1r':
                        c = 'r'
                    elif exp in ['Ice1ra','Ice1rax']:
                        c = 'orange'
                    elif exp in ['Ice1rr','Ice1rrx']:
                        c = 'purple'
                    elif exp == 'Ice2r':
                        c = 'b'
                    elif exp in ['Ice2ra','Ice2rax']:
                        c = 'yellow'
                    elif exp in ['Ice2rr','Ice2rrx']:
                        c = 'pink'

                    # Import the results for this particular experiment
                    if os.path.exists(ncf_output):
                        f = NetCDFFile(ncf_output,'r+')
                        g = NetCDFFile(stats_output,'r+')
                        time = f.variables['daysSinceStart'][:]
                        time /= 365.
                        gArea = g.variables['groundedIceArea'][:]
                        gArea /= (1e6*1e3)
                        # the 'x' runs that go to 1000 years have 10 yr globalStats output and 100 yr for normal output. Grap only every 100 yr
                        if exp[-1] == 'x':
                            print np.shape(gArea)
                            gArea = gArea[::10]
                        f.close()
                        g.close()
                        plt.plot(time,gArea,marker=m,color='k',mfc=c)
                    else:
                        print 'Experiment ', exp, ' did not run, no globalStats.nc file'

                    os.chdir('../')
                os.chdir('../')

            # Create a legend with colors corresponding to the different experiments
            import matplotlib.patches as mpatches
            p1 = mpatches.Patch(color='grey')
            p2 = mpatches.Patch(color='r')
            p3 = mpatches.Patch(color='orange')
            p4 = mpatches.Patch(color='purple')
            p5 = mpatches.Patch(color='b')
            p6 = mpatches.Patch(color='yellow')
            p7 = mpatches.Patch(color='pink')
            plt.legend([p1,p2,p3,p4,p5,p6,p7],['Ice0','Ice1r','Ice1ra','Ice1rr','Ice2r','Ice2ra','Ice2rr'],loc=3,ncol=2,prop={'size':14},framealpha=1.)

            plt.ylabel(r'Grounded Area (1000 km$^2$)',size=20)
            plt.xlabel('Time (yrs)',size=20)
            plt.axvline(100,color='k',ls='--',alpha=0.5)
            plt.axvline(200,color='k',ls='--',alpha=0.5)
            plt.ylim(23,40)
            plt.xlim(-10,1000)
            plt.xticks(size=20)
            plt.yticks(size=20)

            # this is *approximately* the region of the retrograde slope
            plt.fill_between(np.linspace(-100,20000,10),np.ones(10)*35,np.ones(10)*42,color='k',alpha=0.1)

            plt.draw()
            plt.show()
            #plt.savefig('Mismip_experiments.png',dpi=300)
            plt.cla()
            os.chdir('../')
