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

startsnap = 305
endsnap = 306
npzfile = '/ufrc/narayanan/s.lower/simba_snap305_pos_test.npz'

model_dir_base = '/ufrc/narayanan/desika.narayanan/paper/powderday_paper/codes/sfr_lir/'
hydro_dir = '/orange/narayanan/s.lower/simba/filtered_snapshots/'

hydro_outputfile = '/ufrc/narayanan/pg3552/gizmo/output_time/output_m25.txt'

#if we want to write the files locally, but have the paths in the
#parameters files lead to differnet paths (for a different computer),
#put those paths here.  otherweise, set these equal to whatever is in
#model_dir and hydro_dir
hydro_dir_remote = hydro_dir

model_run_name='simba_m50'
COSMOFLAG=0 #flag for setting if the gadget snapshots are broken up into multiples or not and follow a nomenclature snapshot_000.0.hdf5
FILTERFLAG = 1 #flag for setting if the gadget snapshots are filtered or not, and follow a nomenclature snap305_galaxy1800_filtered.hdf5


SPHGR_COORDINATE_REWRITE = True


#===============================================

if (COSMOFLAG == 1) and (FILTERFLAG == 1):
    raise ValueError("COSMOFLAG AND FILTER FLAG CAN'T BOTH BE SET")


data = np.load(npzfile,allow_pickle=True)
pos = data['pos'][()] #positions dictionary

#ngalaxies is the dict that says how many galaxies each snapshot has, in case it's less than NGALAXIES_MAX
ngalaxies = data['ngalaxies'][()]


#get the scalefactor and redshifts
snaps_to_scalefactor = {}
snaps_to_redshift = {}
redshift_to_snaps = {}
scalefactor = np.loadtxt(hydro_outputfile)
for snap in range(len(scalefactor)):
    snaps_to_scalefactor[str(snap)] = scalefactor[snap]
    snaps_to_redshift[str(snap)] = (1./scalefactor[snap])-1.
    redshift_to_snaps[str((1./scalefactor[snap])-1.)] = snap


#first call the initial setup_all_cluster shell

for snap in range(startsnap,endsnap):

    model_dir = model_dir_base+'/snap{:03d}'.format(snap)
    model_dir_remote = model_dir
    
    redshift = snaps_to_redshift[str(snap)]
    tcmb = 2.73*(1.+redshift)

    NGALAXIES = ngalaxies['snap{:03d}'.format(snap)]
    
    for nh in range(NGALAXIES):
        
        xpos = pos['galaxy'+str(nh)]['snap{:03d}'.format(snap)][0]
        ypos = pos['galaxy'+str(nh)]['snap{:03d}'.format(snap)][1]
        zpos = pos['galaxy'+str(nh)]['snap{:03d}'.format(snap)][2]

        #figure out tcmb

        cmd = "./cosmology_setup_all_cluster.hipergator.sh "+str(nnodes)+' '+model_dir+' '+hydro_dir+' '+model_run_name+' '+str(COSMOFLAG)+' '+str(FILTERFLAG)+' '+model_dir_remote+' '+hydro_dir_remote+' '+str(xpos)+' '+str(ypos)+' '+str(zpos)+' '+str(nh)+' '+str(snap)+' '+str(tcmb)
        #print cmd
        call(cmd,shell=True)

        
