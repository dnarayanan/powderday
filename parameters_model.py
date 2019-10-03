#snapshot parameters

snapshot_num = 20
snapnum_str = '{:03d}'.format(snapshot_num)

hydro_dir = '/Volumes/pegasus/gadgetruns/SIGS/G2/'

snapshot_name = 'snapshot_'+snapnum_str+'.hdf5'

#where the files should go
PD_output_dir = '/Volumes/pegasus/pd_runs/test/'
Auto_TF_file = 'snap'+snapnum_str+'.logical'
Auto_dustdens_file = 'snap'+snapnum_str+'.dustdens'


#===============================================
#FILE I/O
#===============================================
inputfile = PD_output_dir+'/example.'+snapnum_str+'.rtin'
outputfile = PD_output_dir+'/example.'+snapnum_str+'.rtout'


#===============================================
#GRID POSITIONS
#===============================================
x_cent = 0
y_cent = 0
z_cent = 0

TCMB = 2.73
