#GRIDDING PARAMETERS

#snapshot parameters

GADGET_octree_gen = False
hydro_dir = '/Users/desika/gadgetruns/sbw_tests/mw_18_6_hr_hightimeres/'
Gadget_snap_num = 6
Gadget_snap_name = 'snapshot_006.hdf5'

YT_octree_gen = True



#where the files should go

PD_output_dir = '/Users/desika/Dropbox/powderday/pd_runs/sbw_tests/mw_18_6_hr_hightimeres/'

Auto_TF_file = 'snap6.logical'
Auto_positions_file = 'snap6.positions'
Auto_dustdens_file = 'snap6.dustdens'


#===============================================
#PARALLELIZATION
#===============================================

n_processes = 4

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
#dustfile = 'dustfiles/kmh_lite.hdf5'
dustfile = '/Users/desika/hyperion-dust-0.1.0/dust_files/d03_4.0_4.0_A.hdf5'
dusttometals_ratio = 0.4



#===============================================
#HYDRO CODE UNITS
#===============================================

stellar_softening_length = 0.005 #kpc - typically the softening length in your SPH calculation
unit_mass = 1.e10 #msun
unit_length = 1. #kpc


#===============================================
#STELLAR SEDS INFO
#===============================================

COSMOFLAG = False

disk_stars_age = 1 #Gyr ;meaningless if COSMOFLAG = True
bulge_stars_age = 1 #Gyr ; meaningless if COSMOFLAG = True

#bins for binning the stellar ages and metallicities for SED
#assignments in cases of many (where many ==
#>N_METALLICITY_BINS*N_STELLAR_AGE_BINS) stars; this is necessary for
#reduction of memory load; see manual for details.

N_METAL_BINS = 5
N_STELLAR_AGE_BINS = 500
N_MASS_BINS = 5  #these you want the most of for sure




#===============================================
#IMAGES AND SED
#===============================================

NTHETA = 10

#===============================================
#DEBUGGING
#===============================================

SOURCES_IN_CENTER = False
STELLAR_SED_WRITE = True
SUPER_SIMPLE_SED = False #just generate 1 oct of 100 pc on a side,
                         #centered on [0,0,0].  sources are added at
                         #random positions.
SKIP_GRID_READIN = False

CONSTANT_DUST_GRID = False #if set, then we don't create a dust grid by
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

N_viewing_angles = 10
n_wav = 250    # number of wavelengths in SED
wav_min = 0.01 # min wavelength for SED in micron
wav_max = 5000.# max wavelength for SED in micron



VERBOSE = False

