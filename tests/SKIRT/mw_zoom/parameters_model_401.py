#snapshot parameters

snapshot_num = 134

if snapshot_num < 10:
    snapnum_str = '00'+str(snapshot_num)
elif snapshot_num >= 10 and snapshot_num <100:
    snapnum_str = '0'+str(snapshot_num)
else:
    snapnum_str = str(snapshot_num)

hydro_dir = '/ufrc/narayanan/desika.narayanan/powderday_safe_packages/pd/examples/gadget/mw_zoom/'

snapshot_name = 'snapshot_'+snapnum_str+'.hdf5'

#where the files should go
PD_output_dir = '/ufrc/narayanan/desika.narayanan/pd/tests/SKIRT/mw_zoom/'
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
