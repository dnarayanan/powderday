#snapshot parameters

Gadget_snap_num = 20

if Gadget_snap_num < 10:
    snapnum_str = '00'+str(Gadget_snap_num)
elif Gadget_snap_num >= 10 and Gadget_snap_num <100:
    snapnum_str = '0'+str(Gadget_snap_num)
else:
    snapnum_str = str(Gadget_snap_num)

hydro_dir = '/data/desika/gadgetruns/m13_mr_Dec16_2013/snapdir_'+snapnum_str+'/'

Gadget_snap_name = 'snapshot_'+snapnum_str+'.0.hdf5'

#where the files should go
PD_output_dir = '/home/desika/Dropbox/powderday/pd_runs/m13_mr_Dec16_2013/'
Auto_TF_file = 'snap'+snapnum_str+'.logical'
Auto_dustdens_file = 'snap'+snapnum_str+'.dustdens'


#===============================================
#FILE I/O
#===============================================
inputfile = PD_output_dir+'/example.'+snapnum_str+'.rtin'
outputfile = PD_output_dir+'/example.'+snapnum_str+'.rtout'

