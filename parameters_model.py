#snapshot parameters

Gadget_snap_num = 239
hydro_dir = '/data/desika/gadgetruns/m13_mr_Dec16_2013/snapdir_'+str(Gadget_snap_num)+'/'

Gadget_snap_name = 'snapshot_'+str(Gadget_snap_num)+'.0.hdf5'

#where the files should go
PD_output_dir = '/home/desika/Dropbox/powderday/pd_runs/m13_mr_Dec16_2013/'
Auto_TF_file = 'snap'+str(Gadget_snap_num)+'.logical'
Auto_dustdens_file = 'snap'+str(Gadget_snap_num)+'.dustdens'


#===============================================
#FILE I/O
#===============================================
inputfile = PD_output_dir+'/example.'+str(Gadget_snap_num)+'.rtin'
outputfile = PD_output_dir+'/example.'+str(Gadget_snap_num)+'.rtout'

