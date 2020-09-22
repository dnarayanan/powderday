#snapshot parameters

snapshot_num = 600
snapnum_str = '{:03d}'.format(snapshot_num)

hydro_dir = '/blue/narayanan/desika.narayanan/powderday_files/gizmo/fire1/FIRE_M12i_ref11/'

snapshot_name = 'snapshot_'+snapnum_str+'.hdf5'

#where the files should go
PD_output_dir = '/home/desika.narayanan/pd_git/tests/SKIRT/gizmo_fire'
Auto_TF_file = 'snap'+snapnum_str+'.logical'
Auto_dustdens_file = 'snap'+snapnum_str+'.dustdens'


#===============================================
#FILE I/O
#===============================================
inputfile = PD_output_dir+'/pd_skirt_comparison.'+snapnum_str+'.rtin'
outputfile = PD_output_dir+'/pd_skirt_comparison.'+snapnum_str+'.rtout'

#===============================================
#GRID POSITIONS
#===============================================
x_cent=29286.13137399
y_cent=31049.24813692
z_cent=32589.56154287

#===============================================
#CMB
#===============================================
TCMB = 2.73
