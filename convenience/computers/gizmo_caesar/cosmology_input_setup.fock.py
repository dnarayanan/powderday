#purpose: to set up torque files and model *.py files from the
#positions written by caesar_cosmology_npzgen.py for a cosmological
#simulation.  This is written for the Haverford College's Fock cluster

import numpy as np
from subprocess import call
import pdb,ipdb
import caesar

#===============================================
#MODIFIABLE HEADER
#===============================================
#shell scripting
nnodes=1

startsnap = 135
endsnap = 136
npzfile = '/Volumes/pegasus/pd_runs/mufasa/m25n512/fh_qr/irxbeta/mufasa_m25n512.halos_pos_for_pd.npz'

model_dir = '/Volumes/pegasus/pd_runs/mufasa/m25n512/fh_qr/smg_survey/'
hydro_dir = '/Volumes/pegasus2/gizmo_runs/mufasa/m25n512/fh_qr/output/'

#if we want to write the files locally, but have the paths in the
#parameters files lead to differnet paths (for a different computer),
#put those paths here.  otherweise, set these equal to whatever is in
#model_dir and hydro_dir
model_dir_remote = '/astro/desika/pd_runs/mufasa/m25n512/fh_qr/smg_survey/'
hydro_dir_remote = '/astro/desika/gizmo_runs/mufasa/m25n512/fh_qr/output/'

model_run_name='mufasa_m25n512'
COSMOFLAG=0 #flag for setting if the gadget snapshots are broken up into multiples or not



SPHGR_COORDINATE_REWRITE = True


#===============================================

data = np.load(npzfile)
pos = data['pos'][()] #positions dictionary
NHALOS = data['NHALOS']


#first call the initial setup_all_cluster shell

for snap in range(startsnap,endsnap):
    for nh in range(NHALOS):
        
        xpos = pos['halo'+str(nh)]['snap{:03d}'.format(snap)][0]
        ypos = pos['halo'+str(nh)]['snap{:03d}'.format(snap)][1]
        zpos = pos['halo'+str(nh)]['snap{:03d}'.format(snap)][2]
        
        cmd = "./cosmology_setup_all_cluster.fock.sh "+str(nnodes)+' '+model_dir+' '+hydro_dir+' '+model_run_name+' '+str(COSMOFLAG)+' '+model_dir_remote+' '+hydro_dir_remote+' '+str(xpos)+' '+str(ypos)+' '+str(zpos)+' '+str(nh)+' '+str(snap)+' '+str(NHALOS)
        #print cmd
        call(cmd,shell=True)

        
