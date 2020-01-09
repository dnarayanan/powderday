#snapshot parameters

snapnum_str = '00300'
hydro_dir = '/ufrc/narayanan/desika.narayanan/yt_datasets/TipsyGalaxy/'

snapshot_name = 'galaxy.00300'

#where the files should go
PD_output_dir = '/ufrc/narayanan/desika.narayanan/pd_git/tests/SKIRT/changa_disk/'
Auto_TF_file = 'changa_test.logical'
Auto_dustdens_file = 'changa_test.dustdens'


#===============================================
#FILE I/O
#===============================================
inputfile = PD_output_dir+'/pd_skirt_comparison.changa.rtin'
outputfile = PD_output_dir+'/pd_skirt_comparison.changa.rtout'

#===============================================
#GRID POSITIONS
#===============================================
x_cent=0.
y_cent=0.
z_cent=0.

#===============================================
#CMB
#===============================================
TCMB = 2.73
