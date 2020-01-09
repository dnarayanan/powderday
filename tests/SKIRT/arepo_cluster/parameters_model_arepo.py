#snapshot parameters

snapnum_str = '59'
hydro_dir = '/ufrc/narayanan/desika.narayanan/yt_datasets/TNGHalo/'

snapshot_name = 'halo_59.hdf5'

#where the files should go
PD_output_dir = '/ufrc/narayanan/desika.narayanan/pd_git/tests/SKIRT/arepo_cluster/'
Auto_TF_file = 'halo_59.logical'
Auto_dustdens_file = 'halo_59.dustdens'


#===============================================
#FILE I/O
#===============================================
inputfile = PD_output_dir+'/pd_skirt_comparison.'+snapnum_str+'.rtin'
outputfile = PD_output_dir+'/pd_skirt_comparison.'+snapnum_str+'.rtout'

#===============================================
#GRID POSITIONS
#===============================================
x_cent=48669.34
y_cent= 53984.04
z_cent= 62114.9028077

#===============================================
#CMB
#===============================================
TCMB = 2.73
