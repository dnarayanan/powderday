#purpose: to set up slurm files and model *.py files from the
#positions written by caesar_cosmology_npzgen.py for a cosmological
#simulation.  This is written for the University of Florida's
#HiPerGator2 cluster.

import numpy as np
from subprocess import call
import pdb,ipdb
import caesar

#===============================================
#MODIFIABLE HEADER
#===============================================
#shell scripting
nnodes=1

snapshot = '/ufrc/narayanan/desika.narayanan/powderday_files/gizmo/fire2/snapshot_600.0.hdf5'
snap = 600
galaxy = 0
x_cent=29338.17963673
y_cent=30979.87927007
z_cent=32479.89844843
TCMB = 2.73



model_dir_base = '/ufrc/narayanan/desika.narayanan/paper/powderday_paper/codes/fire_images/high_fidelity_images/'
hydro_dir = '/ufrc/narayanan/desika.narayanan/powderday_files/gizmo/fire2//'


#if we want to write the files locally, but have the paths in the
#parameters files lead to differnet paths (for a different computer),
#put those paths here.  otherweise, set these equal to whatever is in
#model_dir and hydro_dir
hydro_dir_remote = hydro_dir

model_run_name='simba_m50'
COSMOFLAG=1 #flag for setting if the gadget snapshots are broken up into multiples or not and follow a nomenclature snapshot_000.0.hdf5
FILTERFLAG = 0 #flag for setting if the gadget snapshots are filtered or not, and follow a nomenclature snap305_galaxy1800_filtered.hdf5


SPHGR_COORDINATE_REWRITE = True

nprocessors = 100

#===============================================

model_dir = model_dir_base
model_dir_remote = model_dir
    
tcmb = TCMB

xpos = x_cent
ypos = y_cent
zpos = z_cent




for processor_number in range(nprocessors):
    
    cmd = "./high_fidelity_image_array.sh "+str(nnodes)+' '+model_dir+' '+hydro_dir+' '+model_run_name+' '+str(COSMOFLAG)+' '+str(FILTERFLAG)+' '+model_dir_remote+' '+hydro_dir_remote+' '+str(xpos)+' '+str(ypos)+' '+str(zpos)+' '+str(galaxy)+' '+str(snap)+' '+str(tcmb)+' '+str(processor_number)
    #print cmd
    call(cmd,shell=True)

    pdb
