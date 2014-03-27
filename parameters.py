#GRIDDING PARAMETERS

Grid_Type = 'Octree'  #your choices here are 'Octree' and 'Cart'; we'll work more in in time.


Manual_TF = False
Manual_TF_file = '/Users/desika/Dropbox/powderday/grid_data/sunrise_logical.6.temp_interpolate.dat2'
Manual_density_file = '/Users/desika/Dropbox/powderday/grid_data/mw_18_6.snap6.dustdens'

GADGET_octree_gen = False
Gadget_dir = '/Users/desika/Dropbox/powderday/grid_data/'
Gadget_snap_num = 6
Gadget_snap_name = '/Users/desika/Dropbox/powderday/grid_data/snapshot_006.hdf5'

YT_octree_gen = True


#file for writing the grid if it doesn't already exist

Auto_TF_file = 'grid_data/mw_18_6.snap6.logical'
Auto_positions_file = 'grid_data/mw_18_6.snap6.positions'
Auto_dustdens_file = 'grid_data/mw_18_6.snap6.dustdens'




#===============================================
#PARTICLE SMOOTHING PARAMETERS
#===============================================
NCHUNK = 10. #number of particles per chunk 
NPARTICLES_DEBUG = -1 #ONLY FOR DEBUGGING. if set, we set this to be
                       #the number of particles that we smooth onto
                       #the grid. this needs to be set to -1 for all
                       #production runs.  if this happens to be > the
                       #the number of actual particles, the number of
                       #actual particles is what will be used.


#===============================================
#DUST INFORMATION
#===============================================
dustfile = 'dustfiles/kmh_lite.hdf5'


#===============================================
#HYDRO CODE UNITS
#===============================================

stellar_softening_length = 0.005 #kpc - typically the softening length in your SPH calculation
unit_mass = 1.e10 #msun
unit_length = 1. #kpc


#===============================================
#OLD STARS INFO
#===============================================
disk_stars_age = 1 #Gyr
bulge_stars_age = 1 #Gyr



#===============================================
#IMAGES AND SED
#===============================================

NTHETA = 2

#===============================================
#DEBUGGING
#===============================================

SOURCES_IN_CENTER = False
STELLAR_SED_WRITE = True
SUPER_SIMPLE_SED = False #just generate 1 oct of 100 pc on a side,
                         #centered on [0,0,0].  sources are added at
                         #random positions.
SKIP_GRID_READIN = False
NEW_STARS_ONLY = False # if set, we don't use oldstars and disk stars

CONSTANT_DUST_GRID = True #if set, then we don't create a dust grid by
                          #smoothing, but rather just make it the same
                          #size as the octree with a constant value of
                          #4e-20

#===============================================
#GRID INFORMATION
#===============================================

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

