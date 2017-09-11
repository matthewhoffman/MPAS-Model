#!/usr/bin/env python

'''
This script runs a set of simulations for the landice_model
with varying parameters/forcings for each simulation.
Each simulation writes its own set of output files in a new directory.

It was written to run a parameter sensitivity test on the MISMIP+ spinup.

Change the desired namelist/albany/netcdf parameters on lines 59-61 below.
Decide how to run the model (either on LANL HPC or on a desktop) on lines 140-143.

# Author: Ben Hills, Summer 2017
'''

import numpy as np
from netCDF4 import Dataset as NetCDFFile
import subprocess
import os,sys,glob
from shutil import copyfile

########################################################################################

### Setup

# input files
f_albany_in = '../albany_input.xml'
f_namelist_in = '../namelist.landice.Spinup'
f_streams_in = '../streams.landice.Spinup'
f_netcdf_in = '../landice_grid.nc'
f_slurm_in = '../slurm.run.Spinup'
for f in [f_albany_in,f_namelist_in,f_netcdf_in]:
    if not os.path.exists(f[1:]):
        print "File ", f[1:], "does not exist, exiting..."
        sys.exit()

### Change the parameters below depending on the desired variables to change in the parameter sensitivity test.

# parameters to be changed in the namelist.landice file, alias for naming the directory
nl_name = 'SGH_conduc_coeff'
nl_alias = 'cond'
nl_init_value = 0.001
#nl_name = 'SGH_bed_roughness'
#nl_alias = 'rough'
#nl_init_value = 0.5
#nl_name = 'SGH_bed_roughness_max'
#nl_alias = 'bump'
#nl_init_value = 0.1

# parameters to be changed in the albany_input.xml file
al_name = '"Glen\'s Law A"'
al_init_value = 2e-5

# forcings to be changed in the landice_grid.nc file
ncf_name = 'C0'
ncf_init_value = 1.0

### Change the arrays below depending on the desired values to start for each variable
# Be careful here because each value for a parameter iterates through values for the other parameters.
# For instance if there are 3 values for nl (namelist), al (albany), and ncf (netcdf) then 27 runs would start.
nl_values = [0.001]
al_values = [1e-05]
ncf_values = [1.0]

########################################################################################

### Define Functions for rewriting files

def writeFile(src,target,f_in,f_out):
    if not os.path.exists(f_out):
        copyfile(f_in,f_out)
    f = open(f_out,'r')
    filedata = f.read()
    f.close()

    newdata = filedata.replace(src,target)

    f = open(f_out,'w')
    f.write(newdata)
    f.close()

def writeNCF(C0,f_netcdf_out):
    # copy the original netCDF file into the new subdirectory
    if not os.path.exists(f_netcdf_out):
        subprocess.check_call(['ncks',f_netcdf_in,f_netcdf_out])

    # open netCDF file
    gridfile = NetCDFFile(f_netcdf_out,'r+')

    # Basal Friction Parameter
    C0 *= 1./((3600.*24.*365.)**(1/3.))
    gridfile.variables['beta'][:] = C0

    gridfile.close()

########################################################################################

### Run the model for a set of parameter/forcing values

# loop through the multiplication factors of conduc_coeff
for al_value in al_values:
    for ncf_value in ncf_values:
        for nl_value in nl_values:

            # make directory and copy all necessary files to it
            dir_name = './A_'+str(al_value)+'_'+nl_alias+'_'+str(nl_value)+'_'+ncf_name+'_'+str(ncf_value)
            if not os.path.isdir(dir_name):
                    os.mkdir(dir_name)
            os.chdir(dir_name)

            # edit files that change with every run
            src = nl_name + ' = ' + str(nl_init_value)
            target = nl_name + ' = ' + str(nl_value)
            writeFile(src,target,f_namelist_in,'./namelist.landice')

            src = 'output_interval="0010-00-00_00:00:00"'
            target = 'output_interval="0100-00-00_00:00:00"'
            writeFile(src,target,f_streams_in,'./streams.landice')

            src = al_name + ' type="double" value="' + '%f' %(al_init_value)
            target = al_name + ' type="double" value="' + '%f' %(al_value)
            writeFile(src,target,f_albany_in,'./albany_input.xml')

            src = '#SBATCH --job-name=spinup'
            target = '#SBATCH --job-name=' + nl_alias + str(nl_value)
            writeFile(src,target,f_slurm_in,'./slurm.run')

            writeNCF(ncf_value,'./landice_grid.nc')

            # symlink files that are the same for every run or those that did not get edited
            os.symlink('../landice_model','./landice_model')
            graph_files = glob.iglob(os.path.join("../graph.info.part.*"))
            for f in graph_files:
                os.symlink(f,f[3:])

            # Run the model on the LANL HPC system
            #subprocess.check_call(['sbatch','slurm.run'])
            # Run the model on a desktop
            #subprocess.check_call(['mpirun','np','-16','./landice_model','&'])

            os.chdir('..')
