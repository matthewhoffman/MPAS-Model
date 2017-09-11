#!/usr/bin/env python

'''
This script runs a set of mismip+ simulations with varying parameters/forcings for each simulation.
For each set of parameters this script will initiate a 'constant_beta' run which uses the betaSolve
field from the spinup and does not let hydrology evolve as well as an 'evolve' run with evolving hydrology.

Note that spinup runs need to be done first.

Change the desired namelist/albany/netcdf parameters on lines 48-50 below.
Decide how to run the simulations (either with LANL HPC or on a desktop) on lines 252-255.

Author: Ben Hills, Summer 2017
'''

import numpy as np
from netCDF4 import Dataset as NetCDFFile
import subprocess
import os,sys,glob
import shutil
import fileinput

########################################################################################

### Setup, change the parameters below depending on the desired variables to change in the parameter sensitivity test.

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

# MISMIP+ experiments to run
experiments = ['Ice0', 'Ice1r', 'Ice1ra', 'Ice1rr', 'Ice1rax', 'Ice1rrx', 'Ice2r', 'Ice2ra', 'Ice2rr', 'Ice2rax', 'Ice2rrx']


### Define Function for rewriting files

def writeFile(src,target,f_in,f_out):
    if not os.path.exists(f_out):
        shutil.copyfile(f_in,f_out)
    f = open(f_out,'r')
    filedata = f.read()
    f.close()

    newdata = filedata.replace(src,target)

    f = open(f_out,'w')
    f.write(newdata)
    f.close()

########################################################################################

### Run the model for a set of parameter/forcing values

# loop through the multiplication factors of conduc_coeff
for al_value in al_values:
    for ncf_value in ncf_values:
        for nl_value in nl_values:

            # make directory and copy all necessary files to it
            dir_name = 'A_'+str(al_value)+'_'+nl_alias+'_'+str(nl_value)+'_'+ncf_name+'_'+str(ncf_value)
            print "Creating new parameter directory, ", dir_name
            if not os.path.isdir(dir_name):
                    os.mkdir(dir_name)
            os.chdir(dir_name)

            # Check for the spinup run
            if not os.path.isdir('../../spinup/'+dir_name):
                print "These parameters have no spinup, run that first."
                continue

            # Copy the restart file from the spinup run with the same parameters
            with open('../../spinup/'+dir_name+'/restart_timestamp') as f:
                restart_time = f.readline()[1:6]
            subprocess.check_call(['ncks','-O','-x','-v','xtime',
                '../../spinup/'+dir_name+'/restart_'+restart_time+'.nc','landice_grid_evolve.nc'])
            subprocess.check_call(['ncks','-O','-d','Time,180','-v',
                'betaSolve','../../spinup/'+dir_name+'/output_00000.nc','betaSolve.nc'])
            subprocess.check_call(['ncks','-O','-x','-v',
                'beta','landice_grid_evolve.nc','landice_grid_constant_beta.nc'])
            subprocess.check_call(['ncrcat','-A','betaSolve.nc','landice_grid_constant_beta.nc'])
            subprocess.check_call(['ncrename','-v','betaSolve,beta','landice_grid_constant_beta.nc'])

            # Make a template directory with all of the namelist/streams files so that many sets MISMIP+ experiments can be run
            template_dir = 'mismip+_experiments_template/'

            # Create a separate run directory for evolving hydrology and for a constant beta case
            for run in ['evolve','constant_beta']:
                print "Run Mismip+ experiments in ", run, "directory"
                if not os.path.isdir(run):
                    os.mkdir(run)
                os.chdir(run)

                # edit files that change with every run
                for expt in experiments:
                    src = nl_name + ' = ' + str(nl_init_value)
                    target = nl_name + ' = ' + str(nl_value)
                    writeFile(src,target,'../../'+template_dir+'namelist.landice.'+expt,'./namelist.landice.'+expt)
                    # Turn off the SGH solver and use a constant beta in the 'constant_beta' runs
                    if run == 'constant_beta':
                        src = 'config_beta_use_effective_pressure = .true.'
                        target = 'config_beta_use_effective_pressure = .false.'
                        writeFile(src,target,'./namelist.landice.'+expt,'./namelist.landice.'+expt)
                        src = 'config_SGH = .true.'
                        target = 'config_SGH = .false.'
                        writeFile(src,target,'./namelist.landice.'+expt,'./namelist.landice.'+expt)
                    shutil.copyfile('../../'+template_dir+'streams.landice.'+expt,'./streams.landice.'+expt)

                # Write albany and slurm files
                src = al_name + ' type="double" value="' + '%f' %(al_init_value)
                target = al_name + ' type="double" value="' + '%f' %(al_value)
                writeFile(src,target,'../../albany_input.xml','./albany_input.xml')
                src = '#SBATCH --job-name=experiments'
                target = '#SBATCH --job-name=' + run
                writeFile(src,target,'../../slurm.run.experiments','./slurm.run')

                # symlink files that are the same for every run or those that did not get edited
                if not os.path.exists('./landice_model'):
                    os.symlink('../../landice_model','./landice_model')
                graph_files = glob.iglob(os.path.join("../../graph.info.part.*"))
                for f in graph_files:
                    if not os.path.exists(f):
                        os.symlink(f,f[3:])

                ##########################################################################################################

                # Set up each subdirectory, looping through experiments
                for expt in experiments:
                    print 'Setting up directory for experiment', expt

                    # Make the subdirectory if it does not exist already
                    try:
                        os.mkdir(expt)
                    except:
                        pass

                    # Go to the subdirectory
                    try:
                        os.chdir(expt)
                    except:
                        sys.exit('Error, could not change to subdirectory')
                    # Move the appropriate namelist and stream files from the parent directory.
                    # Note: In the subdirectory, the expt prefix (e.g., 'Ice0') is not included.
                    #       So there is no need for the -n and -s specifiers when launching a run.
                    namelistFile = '../namelist.landice.' + expt
                    shutil.move(namelistFile, './namelist.landice')

                    streamsFile = '../streams.landice.' + expt
                    shutil.move(streamsFile, './streams.landice')

                    # Link to the executable in the parent directory
                    executableName = 'landice_model'
                    if not os.path.exists('landice_model'):
                        os.symlink('../landice_model', 'landice_model')

                    # Link to any and all graph partition files in the parent directory
                    # Note: If a new file is needed, it can be created using metis.
                    #       For example to run on 128 cores on LANL IC:
                    #       > module load metis
                    #       > gpmetis graph.info 128
                    #       This creates a file called graph.info.part.128

                    for f in os.listdir('..'):
                        if f.startswith('graph.info.part'):
                            if not os.path.exists(f):
                                os.symlink('../' + f, f)

                    # Link to the albany input file in the parent directory
                    if not os.path.exists('albany_input.xml'):
                        os.symlink('../' + 'albany_input.xml', 'albany_input.xml')

                    # Link to the appropriate restart file and timestamp
                    # No restart file needed for the Spinup experiment
                    # Note: The symlinks will initially be empty (except for landice_grid.nc).
                    #       Ice0, Ice1r and Ice2r must follow the Spinup.
                    #       Ice1ra and Ice1rr must follow Ice1r; Ice2ra and Ice2rr must follow Ice2r.
                    #       Ice1rax must follow Ice1ra, and similiary for the other Ice*x.


                    if expt =='Ice0' or expt=='Ice1r' or expt=='Ice2r':
                        # Start from restart file at the end of Spinup, but call it landice_grid.nc,
                        #  so the run is treated as a cold start
                        # Note: This requires a one-line NCO command to rename the Spinup restart file
                        #       while removing the xtime variable
                        gridfile = 'landice_grid_'+run+'.nc'
                        griddir = '../../'
                        if not os.path.exists('landice_grid.nc'):
                            os.symlink(griddir + gridfile, 'landice_grid.nc')
                    else:
                        # Start from the appropriate restart file
                        if expt=='Ice1ra' or expt=='Ice1rr':
                            restartYear = 100
                            restartdir = '../Ice1r/'
                        elif expt=='Ice1rax':
                            restartYear = 200
                            restartdir = '../Ice1ra/'
                        elif expt=='Ice1rrx':
                            restartYear = 200
                            restartdir = '../Ice1rr/'
                        elif expt=='Ice2ra' or expt=='Ice2rr':
                            restartYear = 100
                            restartdir = '../Ice2r/'
                        elif expt=='Ice2rax':
                            restartYear = 200
                            restartdir = '../Ice2ra/'
                        elif expt=='Ice2rrx':
                            restartYear = 200
                            restartdir = '../Ice2rr/'

                        # Link to the restart file
                        restartfile = 'restart_00' + str(restartYear) + '.nc'
                        try:
                            os.symlink(restartdir + restartfile, restartfile)
                        except:
                            os.remove(restartfile)
                            os.symlink(restartdir + restartfile, restartfile)

                        # Create the restart_timestamp file
                        # Not using symbolic links because these allow files to be rewritten
                        #  from other directories
                        timestampFile = open('restart_timestamp', 'w')
                        restartTimestamp = ' ' + str(restartYear) + '-01-01_00:00:00'
                        timestampFile.write(restartTimestamp + '\n')
                        timestampFile.close()

                    # Go back to the main directory and continue
                    os.chdir('..')

                ########################################################################################################

                # Run the model on the LANL HPC system
                #subprocess.check_call(['sbatch','slurm.run'])
                # Run the model on a desktop
                #subprocess.check_call(['mpirun','np','-16','./landice_model','&'])

                os.chdir('..')
            os.chdir('..')


