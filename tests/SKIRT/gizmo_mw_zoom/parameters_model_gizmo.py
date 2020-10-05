#snapshot parameters

snapshot_num = 134
snapnum_str = '{:03d}'.format(snapshot_num)

hydro_dir = '/blue/narayanan/desika.narayanan/powderday_files/gizmo/mufasa/'

snapshot_name = 'snapshot_'+snapnum_str+'.hdf5'

#where the files should go
PD_output_dir = '/home/desika.narayanan/pd_git/tests/SKIRT/gizmo_mw_zoom/'
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
x_cent=28077.51291590734
y_cent=22251.040921986722
z_cent=28531.439944942518

#===============================================
#CMB
#===============================================
TCMB = 2.73
