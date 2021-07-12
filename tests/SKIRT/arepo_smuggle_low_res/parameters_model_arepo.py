#snapshot parameters

snapnum_str = '143'
hydro_dir = '/blue/narayanan/desika.narayanan/powderday_files/smuggle/low_res/'

galaxy_num = 0
galaxy_num_str = str(galaxy_num)

snapshot_name = 'smuggle_snapshot_143.low_res.hdf5'

#where the files should go
PD_output_dir = '/home/desika.narayanan/pd_git/tests/SKIRT/arepo_smuggle_low_res/'
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
x_cent=300
y_cent= 300
z_cent= 300

#===============================================
#CMB
#===============================================
TCMB = 2.73
