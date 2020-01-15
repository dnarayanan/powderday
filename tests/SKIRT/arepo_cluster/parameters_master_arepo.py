# coding=utf-8

#===============================================
#HOME INFORMATION
#===============================================
pd_source_dir ='/ufrc/narayanan/desika.narayanan/pd_git/'

#===============================================
#RESOLUTION KEYWORDS
#===============================================
oref = 0 #over refine factor - should typically be set to 0
n_ref = 128 #when n_particles > n_ref, octree refines further
zoom_box_len = 250 #kpc; so the box will be +/- zoom_box_len from the center
bbox_lim = 1.e5 #kpc - this is the initial bounding box of the grid (+/- bbox_lim)
               #This *must* encompass all of the particles in the
               #simulation. 

#===============================================
#PARALLELIZATION
#===============================================

n_processes = 16 #number of pool processes to run
n_MPI_processes = 1 #number oF MPI processes to run

#===============================================
#RT INFORMATION
#===============================================
n_photons_initial = 1.e5
n_photons_imaging = 1.e5
n_photons_raytracing_sources = 1.e5
n_photons_raytracing_dust = 1.e5

FORCE_RANDOM_SEED = False
seed = -12345 #has to be an int, and negative.

#===============================================
#DUST INFORMATION 
#===============================================
dustdir = '/home/desika.narayanan/hyperion-dust-0.1.0/dust_files/' #location of your dust files
dustfile = 'd03_3.1_6.0_A.hdf5'
PAH = False
dust_grid_type = 'dtm' #needs to be in ['dtm','rr','manual','li_bestfit','li_ml']
dusttometals_ratio = 0.25 #skirt assumes 0.25: see http://www.skirt.ugent.be/tutorials/_tutorial_hydro_s_p_h.html ("dust system"subheading) 
enforce_energy_range = False #False is the default;  ensures energy conservation


SUBLIMATION = True #do we automatically kill dust grains above the
                    #sublimation temperature; right now is set to fast
                    #mode
SUBLIMATION_TEMPERATURE = 1600. #K -- meaningliess if SUBLIMATION == False

#===============================================
#STELLAR SEDS INFO
#===============================================
FORCE_BINNING = True #force SED binning

imf_type = 1 #FSPS imf types; 0 = salpeter, 1 = chabrier; 2 = kroupa; 3 and 4 (vandokkum/dave) not currently supported
pagb = 0 #weight given to post agb stars# 1 is the default
add_agb_dust_model=False #add circumstellar AGB dust model (100%); Villaume, Conroy & Jonson 2015


#===============================================
#NEBULAR EMISSION INFO
#===============================================
add_neb_emission = False #add nebular line emission from Cloudy Lookup tables (dev. by Nell Byler)

use_cloudy_tables = True    # If True, CLOUDY look up tables will be used to calculate nebular emission.
                            # Otherwise CLOUDY models are generated individually
                            # for each young star particle (under active development) (Default: True)

FORCE_gas_logu = False      # If set, then we force the ionization parameter (gas_logu) of HII regions to be
                            # gas_logu (next parameter) else, it is taken to be variable and dependent on ionizing
                            # radiation from star particles. (Default: False)

gas_logu = -2. #gas ionization parameter for HII regions; only relevant
              #if add_neb_emission = True default = -2

FORCE_gas_logz = False #if set, then we force the gas_logz of HII
                       #regions to be gas_logz (next parameter); else, it is taken to be the star particles metallicity.  default is False
gas_logz = 0. #units of log(Z/Z_sun); metallicity of the HII region
             #metallicity; only relevant if add_neb_emission = True;
             #default is 0

FORCE_logq = False          # If set, then we force the number of ionizing photons to be source_logq (next parameter)
                            # else, it is taken to be variable and dependent on ionizing radiation
                            # from star particles. (Default: False)

source_logq = 1.e47         # Number of ionizing photons emitted by the source in units of s^-1.
                            # Only relevant if add_neb_emission = True, use_cloudy_tables = True and
                            # FORCE_gas_logq = True (Default: 1.e47)


FORCE_inner_radius = False  # If set, then we force the inner radius of the cloud to be inner_radius (next parameter)
                            # else, it is taken to be the Stromgren sphere radius. (Default: False)

inner_radius = 1.e19        # This sets the inner radius of the cloud in cm. This is used only when add_neb_emission = True,
                            # use_cloudy_tables = True and FORCE_inner_radius = True (Default: 1.e19)


neb_abund = "dopita"        # This sets the HII region elemental abundances for generating CLOUDY models.
                            # Available abundaces are.
                            #    dopita:    Abundabces from Dopita (2001) with old solar abundances = 0.019 and ISM grains.
                            #    newdopita: Abundances from Dopita (2013). Solar Abundances from Grevasse 2010 - z= 0.013
                            #               includes smooth polynomial for N/O, C/O relationship functional form for He(z),
                            #               new depletion and factors in ISM grains.
                            #    gutkin:    Abundabces from Gutkin (2016) and PARSEC metallicity (Bressan+2012) based on Grevesse+Sauvel (1998)
                            #               and Caffau+2011
                            # This is used only when add_neb_emission = True and use_cloudy_tables = True. (Default: dopita)

use_Q = True                # If True, we run CLOUDY by specifying number of ionizing photons which are calculated
                            # based on the input sed and the inner radius which is set to the Strömgren radius.
                            # else, CLOUDY is run by specifying just the ionization parameter.Only relevant if
                            # add_neb_emission = True and use_cloudy_tables = True (Default: True)

HII_T = 1.e4                # Ionized gas temperature in K for calculating nebular emission.
                            # This is used only when add_neb_emission = True (Default = 1.e4)

HII_nh = 1.e2               # Gas hydrogen density for calcualting nebular emission in units if cm^-3.
                            # This is used only when add_neb_emission = True (Default = 1.e2)

HII_max_age = 2.e-3         # Sets the maximum age limit for calculating nebular emission in units of Gyr.
                            # This is used only when add_neb_emission = True (Default = 2.e-3)

neb_file_output = True      # If set to True creates an output file with ionization parameter (LogU),
                            # number of ionizing photons (LogQ), inner radius, stellar mass, age and
                            # metallicity(zmet) for each particle. (Default: True)


stellar_cluster_mass = 1.e4 # Mass of star clusters in Msun. This is used only when add_neb_emission = True (Default = 1.e4)

cloudy_cleanup = True       # If set to True, all the CLOUDY files will be deleted after the source addition is complete.
                            # Only relevant if add_neb_emission = True and use_cloudy_tables = True (Default: True)



CF_on = False #if set to true, then we enable the Charlot & Fall birthcloud models 
birth_cloud_clearing_age = 0.01 #Gyr - stars with age <
                                #birth_cloud_clearing_age have
                                #charlot&fall birthclouds meaningless
                                #of CF_on  == False

Z_init = 0 #force a metallicity increase in the newstar particles.
           #This is useful for idealized galaxies.  The units for this
           #are absolute (so enter 0.02 for solar).  Setting to 0
           #means you use the stellar metallicities as they come in
           #the simulation (more likely appropriate for cosmological
           #runs)

#Idealized Galaxy SED Parameters
disk_stars_age = 8 #Gyr ;meaningless if this is a cosmological simulation; note, if this is <= 7, then these will live in birth clouds
bulge_stars_age = 8 #Gyr ; meaningless if this is a cosmological simulation; note, if this is <= 7, then these will live in birth clouds
disk_stars_metals = 19 #in fsps metallicity units
bulge_stars_metals = 19 #in fsps metallicity units



#bins for binning the stellar ages and metallicities for SED
#assignments in cases of many (where many ==
#>N_METALLICITY_BINS*N_STELLAR_AGE_BINS) stars; this is necessary for
#reduction of memory load; see manual for details.

N_STELLAR_AGE_BINS = 100


metallicity_legend= "/home/desika.narayanan/fsps/ISOCHRONES/Padova/Padova2007/zlegend.dat"

#===============================================
#BLACK HOLE STUFF
#===============================================

BH_SED = False
BH_eta = 0.1 #bhluminosity = BH_eta * mdot * c**2.
BH_model = "Nenkova"
BH_modelfile = "/home/desika.narayanan/powderday/agn_models/clumpy_models_201410_tvavg.hdf5"
# The Nenkova BH_modelfile can be downloaded here:
# https://www.clumpy.org/downloads/clumpy_models_201410_tvavg.hdf5

nenkova_params = [5,30,0,1.5,30,40] #Nenkova+ (2008) model parameters

#===============================================
#IMAGES AND SED
#===============================================

NTHETA = 1
NPHI = 1
SED = True

SED_MONOCHROMATIC = False
FIX_SED_MONOCHROMATIC_WAVELENGTHS = False #if set, then we only use
                                         #nlam wavelengths in the
                                         #range between min_lam and
                                         #max_lam
SED_MONOCHROMATIC_min_lam = 0.3 #micron
SED_MONOCHROMATIC_max_lam = 0.4 #micron
SED_MONOCHROMATIC_nlam = 100 




IMAGING = False
filterdir = '/home/desika.narayanan/powderday/filters/'
filterfiles = [
    'arbitrary.filter',
#    'ACS_F475W.filter',
#    'ACS_F606W.filter',
#    'ACS_F814W.filter',
#    'B_subaru.filter',
]

# Insert additional filter files as above. In bash, the following command 
# formats the filenames for easy copying/pasting.
# $ shopt -s globstar; printf "#    '%s'\n" *.filter


npix_x = 128
npix_y = 128

#experimental and under development - not advised for use
IMAGING_TRANSMISSION_FILTER = False
filter_list = ['filters/irac_ch1.filter']
TRANSMISSION_FILTER_REDSHIFT = 0.001


#===============================================
#OTHER INFORMATION
#===============================================

solar = 0.013
PAH_frac = {'usg': 0.0586, 'vsg': 0.1351, 'big': 0.8063} # values will be normalized to 1

#===============================================
#DEBUGGING
#===============================================

SOURCES_IN_CENTER = False
STELLAR_SED_WRITE = True
SKIP_RT = False #skip radiative transfer (i.e. just read in the grids and maybe write some diagnostics)
SUPER_SIMPLE_SED = False #just generate 1 oct of 100 pc on a side,
                         #centered on [0,0,0].  sources are added at
                         #random positions.
SKIP_GRID_READIN = False

CONSTANT_DUST_GRID = False #if set, then we don't create a dust grid by
                          #smoothing, but rather just make it the same
                          #size as the octree with a constant value of
                          #4e-20
                          
N_MASS_BINS = 1 #this is really just a place holder that exists in
                #some loops to be able to insert some code downstream
                #for spatially varying IMFs.  right now for speed best
                #to set to 1 as it doesn't actually do anything.

FORCE_STELLAR_AGES = False
FORCE_STELLAR_AGES_VALUE = 0.05# Gyr

FORCE_STELLAR_METALLICITIES = False
FORCE_STELLAR_METALLICITIES_VALUE = 0.012 #absolute values (so 0.013 ~ solar)
