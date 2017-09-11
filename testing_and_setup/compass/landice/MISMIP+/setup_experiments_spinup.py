#!/usr/bin/env python

import glob,sys, os
import shutil

# Create a spinup directory and move files there

print "Creating spinup directory"
try:
    os.mkdir('spinup')
except:
    pass
shutil.move('namelist.landice.Spinup','spinup')
shutil.move('streams.landice.Spinup','spinup')
shutil.move('run_spinup_parameter_tests.py','spinup')
shutil.move('slurm.run.Spinup','spinup')
shutil.move('plot_spinup_GL.py','spinup')
shutil.move('landice_grid.nc','spinup')
os.symlink('../plot_outputfile.py','spinup/plot_outputfile.py')
os.symlink('../setup_mismip+_initial_conditions.py','spinup/setup_mismip+_initial_conditions.py')
os.symlink('../landice_model', 'spinup/landice_model')
os.symlink('../albany_input.xml', 'spinup/albany_input.xml')
graph_files = glob.iglob(os.path.join("graph.info*"))
for file in graph_files:
    shutil.copyfile(file,'spinup/'+file)



# Create an experiments directory and move files there

print "Creating experiments directory"
try:
    os.mkdir('experiments')
except:
    pass
shutil.move('run_mismip+_experiments.py','experiments')
shutil.move('slurm.run.experiments','experiments')
shutil.move('plot_mismip+_experiments_GL.py','experiments')
os.symlink('../plot_outputfile.py','experiments/plot_outputfile.py')
os.symlink('../setup_mismip+_initial_conditions.py','experiments/setup_mismip+_initial_conditions.py')
os.symlink('../landice_model', 'experiments/landice_model')
os.symlink('../albany_input.xml', 'experiments/albany_input.xml')
for file in graph_files:
    shutil.move(file,'experiments/'+file)

# Within the esperiments directory make a template for the namelist/streams files that can be used over and over

print "Creating template for experiments"
experiments = ['Ice0', 'Ice1r', 'Ice1ra', 'Ice1rr', 'Ice1rax', 'Ice1rrx', 'Ice2r', 'Ice2ra', 'Ice2rr', 'Ice2rax', 'Ice2rrx']
template_dir = 'mismip+_experiments_template'
os.chdir('experiments')
try:
    os.mkdir(template_dir)
except:
    pass
# Loop through experiments
for expt in experiments:
    # Move the appropriate namelist and stream files from the parent directory.
    shutil.move('../namelist.landice.'+expt,template_dir)
    shutil.move('../streams.landice.'+expt,template_dir)
