#purpose: to set up slurm files and model *.py files from the
#positions written by caesar_cosmology_npzgen.py for a cosmological
#simulation.  This is written for the University of Florida's
#HiPerGator2 cluster.

import numpy as np
from subprocess import call
import pdb,ipdb
import caesar
from str_snap import str_snap

#===============================================
#MODIFIABLE HEADER
#===============================================
#shell scripting
nnodes=1

startsnap = 135
endsnap = 136
npzfile = '/Volumes/pegasus/pd_runs/mufasa/m25n512/fh_qr/attenuation/mufasa_m25n512.halos_pos_for_pd.npz'

model_dir = '/Volumes/pegasus/pd_runs/mufasa/m25n512/fh_qr/attenuation/'
hydro_dir = '/Volumes/pegasus2/gizmo_runs/mufasa/m25n512/fh_qr/output/'

#if we want to write the files locally, but have the paths in the
#parameters files lead to differnet paths (for a different computer),
#put those paths here.  otherweise, set these equal to whatever is in
#model_dir and hydro_dir
model_dir_remote = '/ufrc/narayanan/desika.narayanan/pd_runs/mufasa/m25n512/fh_qr/attenuation/'
hydro_dir_remote = '/ufrc/narayanan/desika.narayanan/gizmo_runs/mufasa/m25n512/fh_qr/output/'

model_run_name='mufasa_m25n512'
COSMOFLAG=0 #flag for setting if the gadget snapshots are broken up into multiples or not



SPHGR_COORDINATE_REWRITE = True


#===============================================

data = np.load(npzfile)
pos = data['pos'][()] #positions dictionary
NHALOS = data['NHALOS']


#first call the initial setup_all_cluster shell


for nh in range(NHALOS):
    for snap in range(startsnap,endsnap):
        
        xpos = pos['halo'+str(nh)]['snap'+str_snap(snap)][0]
        ypos = pos['halo'+str(nh)]['snap'+str_snap(snap)][1]
        zpos = pos['halo'+str(nh)]['snap'+str_snap(snap)][2]
        
        cmd = "./cosmology_setup_all_cluster.hipergator.sh "+str(nnodes)+' '+model_dir+' '+hydro_dir+' '+model_run_name+' '+str(COSMOFLAG)+' '+model_dir_remote+' '+hydro_dir_remote+' '+str(xpos)+' '+str(ypos)+' '+str(zpos)+' '+str(nh)+' '+str(snap)
        #print cmd
        call(cmd,shell=True)

