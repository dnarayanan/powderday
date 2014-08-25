#===============================================
#RESOLUTION KEYWORDS
#===============================================
oref = 0 #over refine factor - should typically be set to 0
n_ref = 64 #when n_particles > n_ref, octree refines further
zoom = True #False = use the entire grid; True = zoom in on the highest density peak
zoom_box_len = 200 #kpc; so the box will be +/- zoom_box_len from the center
bbox_lim = 1.e5 #kpc - this is the initial bounding box of the grid (+/- bbox_lim)
               #This *must* encompass all of the particles in the
               #simulation. if zoom is set, this is obviated; else, is
               #the simulated boxsize.

#===============================================
#PARALLELIZATION
#===============================================

n_processes = 3 #number of MPI processes to run


#===============================================
#RT INFORMATION
#===============================================
n_photons_initial = 1.e8
n_photons_imaging = 1.e8
n_photons_raytracing_sources = 1.e8
n_photons_raytracing_dust = 1.e8


#===============================================
#DUST INFORMATION
#===============================================
dustfile = '/Users/desika/hyperion-dust-0.1.0/dust_files/d03_4.0_4.0_A.hdf5'
dusttometals_ratio = 0.4



#===============================================
#HYDRO CODE UNITS
#===============================================
unit_mass = 1.e10 #msun/h
unit_length = 1. #kpc/h
unit_age = 1. #Gyr
unit_velocity = 1.e5 #cm/s

#===============================================
#STELLAR SEDS INFO
#===============================================
FORCE_BINNING = True #force SED binning
COSMOFLAG = False  #is this a cosmological simulation?

imf_type = 1

disk_stars_age = 8 #Gyr ;meaningless if COSMOFLAG = True; note, if this is <= 7, then these will live in birth clouds
bulge_stars_age = 8 #Gyr ; meaningless if COSMOFLAG = True; note, if this is <= 7, then these will live in birth clouds
disk_stars_metals = 19 #in fsps metallicity units
bulge_stars_metals = 19 #in fsps metallicity units


CF_on = True #if set to true, then we enable the Charlot & Fall birthcloud models 
birth_cloud_clearing_age = 0.01 #Gyr - stars with age <
                                #birth_cloud_clearing_age have
                                #charlot&fall birthclouds meaningless
                                #of CF_on  == False


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
#GRID INFORMATION  #DEPRECATED - no need to edit.
#===============================================

#center cell position #currently deprecated, though will be used for zooming later
x_cent = 0
y_cent = 0
z_cent = 0 

