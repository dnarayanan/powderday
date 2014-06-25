#===============================================
#RESOLUTION KEYWORDS
#===============================================
oref = 0
n_ref = 64
zoom = False
zoom_box_len = 200 #kpc; so the box will be +/- zoom_box_len from the center

#===============================================
#PARALLELIZATION
#===============================================

n_processes = 3

#===============================================
#PARTICLE SMOOTHING PARAMETERS
#===============================================
NCHUNK = 10. #number of particles per chunk   #DEPRECATED; DON'T USE


#===============================================
#RT INFORMATION
#===============================================
n_photons_initial = 1.e7
n_photons_imaging = 1.e7
n_photons_raytracing_sources = 1.e7
n_photons_raytracing_dust = 1.e7


#===============================================
#DUST INFORMATION
#===============================================
dustfile = '/Users/desika/hyperion-dust-0.1.0/dust_files/d03_4.0_4.0_A.hdf5'
dusttometals_ratio = 0.4



#===============================================
#HYDRO CODE UNITS
#===============================================
unit_mass = 1.e10 #msun
unit_length = 1. #kpc


#===============================================
#STELLAR SEDS INFO
#===============================================

COSMOFLAG = False

disk_stars_age = 8 #Gyr ;meaningless if COSMOFLAG = True
bulge_stars_age = 8 #Gyr ; meaningless if COSMOFLAG = True

#bins for binning the stellar ages and metallicities for SED
#assignments in cases of many (where many ==
#>N_METALLICITY_BINS*N_STELLAR_AGE_BINS) stars; this is necessary for
#reduction of memory load; see manual for details.

N_STELLAR_AGE_BINS = 100
N_MASS_BINS = 100  #these you want the most of for sure

metallicity_legend= "/Users/desika/fsps/ISOCHRONES/Padova/Padova2007/zlegend_basel.dat"



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

