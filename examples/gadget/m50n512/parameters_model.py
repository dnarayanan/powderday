#snapshot parameters

Gadget_snap_num = 80

if Gadget_snap_num < 10:
    snapnum_str = '00'+str(Gadget_snap_num)
elif Gadget_snap_num >= 10 and Gadget_snap_num <100:
    snapnum_str = '0'+str(Gadget_snap_num)
else:
    snapnum_str = str(Gadget_snap_num)

hydro_dir = '/ufrc/narayanan/desika.narayanan/pd/examples/gadget/m50n512/'

Gadget_snap_name = 'snapshot_'+snapnum_str+'.hdf5'

#where the files should go
PD_output_dir = '/ufrc/narayanan/desika.narayanan/pd/examples/gadget/m50n512/'
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
x_cent=46830.38634572
y_cent=33034.01080208
z_cent=19261.22960995

#===============================================
#CMB
#===============================================
TCMB = 2.73
