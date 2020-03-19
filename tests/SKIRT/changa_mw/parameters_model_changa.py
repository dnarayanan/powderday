#snapshot parameters

snapnum_str = '004096'
hydro_dir = '/ufrc/narayanan/desika.narayanan/powderday_files/tremmel/starform_example/'

snapshot_name = 'pioneer50h243.1536gst1bwK1BH.004096'

#where the files should go
PD_output_dir = '/home/desika.narayanan/pd_git/tests/SKIRT/changa_mw/'
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
x_cent=-0.41765088
y_cent= -0.31675988
z_cent= 0.12812596

#===============================================
#CMB
#===============================================
TCMB = 2.73
