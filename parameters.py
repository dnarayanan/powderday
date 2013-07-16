#GRIDDING PARAMETERS

Manual_TF = 0
Manual_TF_file = '/Users/desika/powderday/grid_data/mw_18_6.snap6.logical'
Manual_density_file = '/Users/desika/powderday/grid_data/mw_18_6.snap6.dens'

Gadget_dir = '/Users/desika/powderday/grid_data/'
Gadget_snap_num = 006

#file for writing the grid if it doesn't already exist

Auto_TF_file = '/Users/desika/powderday/grid_data/mw_18_6.snap6.pd.logical'
Auto_positions_file = '/Users/desika/powderday/grid_data/mw_18_6.snap6.pd.positions'





#DUST INFORMATION

dustfile = '/Users/desika/hyperion-dust-0.1.0/dust_files/kmh_lite.hdf5'


#GRID INFORMATION
#size in kpc: note - the parent grid corners are [-dx,dx; -dy,dy; -dz,dz]
dx = 100
dy = 100
dz = 100

#center cell position
x_cent = 0
y_cent = 0
z_cent = 0 



#SED INFORMATION

CALCULATE_SED  = 1 #if set to 0, no SED is calculated
N_viewing_angles = 10
n_wav = 250    # number of wavelengths in SED
wav_min = 0.01 # min wavelength for SED in micron
wav_max = 5000.# max wavelength for SED in micron



VERBOSE = False

