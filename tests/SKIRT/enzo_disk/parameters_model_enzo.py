#snapshot parameters

snapshot_num = 134
snapnum_str = '{:03d}'.format(snapshot_num)

hydro_dir = '/ufrc/narayanan/desika.narayanan/yt_datasets/enzo_iso_galaxy/galaxy0030/'

snapshot_name = 'galaxy0030'

#where the files should go
PD_output_dir = '/ufrc/narayanan/desika.narayanan/pd_git/tests/SKIRT/enzo_disk/'
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
x_cent=0.5
y_cent=0.5
z_cent=0.5

#===============================================
#CMB
#===============================================
TCMB = 2.73
