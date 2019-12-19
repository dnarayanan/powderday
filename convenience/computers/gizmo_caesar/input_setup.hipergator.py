from __future__ import print_function

#script intended to set up the slurm files and model*.py files
#assuming caesar have been written out to an npz file, and that the
#system we're running on is the University of Florida's HiPerGator2 cluster

import numpy as np
from subprocess import call
import pdb,ipdb
import caesar

#===============================================
#MODIFIABLE HEADER
#===============================================
#shell scripting
nnodes=2

model_dir = '/ufrc/narayanan/desika.narayanan/pd_runs/mufasa_zooms/m50n512/z0/halo401_dm2_ml12_nq/NSF_AAG_2018_BPT/lines_off/'
hydro_dir = '/ufrc/narayanan/desika.narayanan/gizmo_runs/mufasa_zooms/m50n512/z0/halo401_dm2_ml12_nq/output/'

#model_dir = '/ufrc/narayanan/desika.narayanan/pd_runs/MassiveFIRE/B100_N512_M3e13_TL00004_baryon_toz2_HR_9915/smg_survey/'
#hydro_dir = '/ufrc/narayanan/desika.narayanan/gizmo_runs/MassiveFIRE/SMGs/B100_N512_M3e13_TL00004_baryon_toz2_HR_9915/output/'

HALOS = False #default: if set to True then we center on the Halo center instead of the central galaxy center
localhalo = 0

#if we want to write the files locally, but have the paths in the
#parameters files lead to differnet paths (for a different computer),
#put those paths here.  otherweise, set these equal to whatever is in
#model_dir and hydro_dir
model_dir_remote = model_dir
hydro_dir_remote = hydro_dir

model_run_name='halo401'
COSMOFLAG=0 #flag for setting if the gadget snapshots are broken up into multiples or not



SPHGR_COORDINATE_REWRITE = True

#GAL = 14 #this is the galaxy from SPH_PROGEN we need to find the
         #progenitors for.  to do this you need to have run sphgr on
         #all the snaps. NOTE - sphgr has to be set up for the corret
         #galaxy for this to work!



#===============================================


#first call the initial setup_all_cluster shell

if HALOS == False:
    data = np.load(hydro_dir+'Groups/caesar_physical_properties.local.'+str(localhalo)+'.npz')
else:
    data = np.load(hydro_dir+'Groups/caesar_physical_properties.halos.local.'+str(localhalo)+'.npz')


for snap in data['snaps']:
    startsnap = snap
    endsnap = snap
    #startsnap = np.min(data['snaps'])
    #endsnap = np.max(data['snaps'])

    cmd = "./setup_all_cluster.hipergator.sh "+str(nnodes)+' '+str(startsnap)+' '+str(endsnap)+' '+model_dir+' '+hydro_dir+' '+model_run_name+' '+str(COSMOFLAG)+' '+model_dir_remote+' '+hydro_dir_remote
    print(cmd)
    call(cmd,shell=True)


    if SPHGR_COORDINATE_REWRITE == True:
        if HALOS == False:
            data = np.load(hydro_dir+'Groups/caesar_physical_properties.local.'+str(localhalo)+'.npz')
        else:
            data = np.load(hydro_dir+'Groups/caesar_physical_properties.halos.local.'+str(localhalo)+'.npz')

    
            
        sph_snap = data['snaps'][::-1]
        sph_cmx = data['xpos'][::-1]
        sph_cmy = data['ypos'][::-1]
        sph_cmz = data['zpos'][::-1]
        sph_redshifts = data['redshifts'][::-1]
        snaps = np.arange(startsnap,endsnap+1)
        
        for i in snaps:
    
            
            wsph = (np.where(sph_snap == i))[0][0]
            x_cent = sph_cmx[wsph]
            y_cent = sph_cmy[wsph]
            z_cent = sph_cmz[wsph]
            redshift = sph_redshifts[wsph]
            TCMB = 2.73*(1.+redshift)
        #append positions
            modelfile = model_dir+'/model_'+str(i)+'.py'
            print('appending coordinates to: ', modelfile)
            with open(modelfile,"a") as myfile:
                myfile.write("\n\n")
                myfile.write("#===============================================\n")
                myfile.write("#GRID POSITIONS\n")
                myfile.write("#===============================================\n")
                myfile.write("x_cent = %s\n" % x_cent)
                myfile.write("y_cent = %s\n" % y_cent)
                myfile.write("z_cent = %s\n" % z_cent)
               
                myfile.write("\n\n")
                myfile.write("#===============================================\n")
                myfile.write("#CMB\n")
                myfile.write("#===============================================\n")
                myfile.write("TCMB = %s\n" %TCMB)
    
                myfile.write("\n")
            
   

