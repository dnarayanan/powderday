#snapshot parameters

snapshot_num = 77
galaxy_num = 0


galaxy_num_str = str(galaxy_num)
snapnum_str = '{:03d}'.format(snapshot_num)

hydro_dir = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m12_n256_active_boxes/output_m12n256_7/'

snapshot_name = 'snapshot_'+snapnum_str+'.hdf5'

#where the files should go
PD_output_dir = '/home/desika.narayanan/pd_git/tests/SKIRT/m12_extinction/'
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
#x_cent=7791.04289063
#y_cent=11030.12925781
#z_cent=4567.8734375
 
x_cent = 5527.42349609  
y_cent = 4680.46623047 
z_cent = 10704.73

#x_cent = 5.00925088e+03 
#y_cent = 1.22012878e+04 
#z_cent = 8.88592731e+00


#x_cent = 10898.06515625 
#y_cent = 10868.67308594  
#z_cent = 4419.96082031

#===============================================
#CMB
#===============================================
TCMB = 2.73
