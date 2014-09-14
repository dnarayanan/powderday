import numpy as np
from subprocess import call
import pdb

#===============================================
#MODIFIABLE HEADER
#===============================================

#shell scripting
nnodes=6
startsnap=50
endsnap=400 #set the same as startsnap if you just want to do one snapshot
model_dir='/data/desika/pd_runs/SIGS/G0/CF_off/8_gyr/dtm_1e-3/'
hydro_dir='/data/desika/gadgetruns/SIGS/G0/'
model_run_name='SIGS_dtm1e-3'
COSMOFLAG=0 #flag


#halo files
halo_dir = '/data/desika/halos/m13_mr_Dec16_2013/'


#===============================================


#first call the initial setup_all_cluster shell


cmd = "./setup_all_cluster.sh "+str(nnodes)+' '+str(startsnap)+' '+str(endsnap)+' '+model_dir+' '+hydro_dir+' '+model_run_name+' '+str(COSMOFLAG)
print cmd
call(cmd,shell=True)



'''
#now get the list of positions: read in the halos

snaps = np.arange(startsnap,endsnap)
for i in snaps:
    if i < 10:
        snapnum_str = '00'+str(i)
    elif i >= 10 and i <100:
        snapnum_str = '0'+str(i)
    else:
        snapnum_str = str(i)

    halo_file = halo_dir+'/halos_'+snapnum_str+'.hop.npz'
    print halo_file

    halo_data = np.load(halo_file)
    #find com of most massive halo
    mass = halo_data['total_halo_mass']
    w = np.where(mass == np.max(mass))[0]
    pos = [halo_data['x_com'][w],halo_data['y_com'][w],halo_data['z_com'][w]]
    pos = np.array(pos)
   

    modelfile = model_dir+'/model_'+str(i)+'.py'
    print 'appending coordinates to: ', modelfile
    with open(modelfile,"a") as myfile:
        myfile.write("\n\n")
        myfile.write("#===============================================\n")
        myfile.write("#GRID CENTERING\n")
        myfile.write("#===============================================\n")
        myfile.write("xcent = %s\n" % float(pos[0]))
        myfile.write("ycent = %s\n" % float(pos[1]))
        myfile.write("zcent = %s\n" % float(pos[2]))


        myfile.write("\n")
'''

    #append to the model file
   


    
