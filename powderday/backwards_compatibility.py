import powderday.config as cfg
import os.path

def variable_set():
    try:
        cfg.par.FORCE_RANDOM_SEED
    except:
        cfg.par.FORCE_RANDOM_SEED  = None

    try:
        cfg.par.BH_SED
    except:
        cfg.par.BH_SED  = None

    try:
        cfg.par.IMAGING
    except:
        cfg.par.IMAGING  = None

    try:
        cfg.par.SED
    except:
        cfg.par.SED = True

    try:
        cfg.par.IMAGING_TRANSMISSION_FILTER
    except:
        cfg.par.IMAGING_TRANSMISSION_FILTER = False

    try:
        cfg.par.SED_MONOCHROMATIC
    except:
        cfg.par.SED_MONOCHROMATIC = False

    try:
        cfg.par.SKIP_RT
    except:
        cfg.par.SKIP_RT = False


    try:
        cfg.par.FIX_SED_MONOCHROMATIC_WAVELENGTHS 
    except:
        cfg.par.FIX_SED_MONOCHROMATIC_WAVELENGTHS = True


    try:
        cfg.par.n_MPI_processes
    except:
        cfg.par.n_MPI_processes = cfg.par.n_processes #default is to make the same as the number of pool processes

        
    try:
        cfg.par.SOURCES_RANDOM_POSITIONS
    except:
        cfg.par.SOURCES_RANDOM_POSITIONS = False
        
    try:
        cfg.par.FORCE_gas_logu
    except:
        cfg.par.FORCE_gas_logu = False

    try:
        cfg.par.gas_logu
    except:
        cfg.par.gas_logu = -2

    try:
        cfg.par.gas_logu_init
    except:
        cfg.par.gas_logu_init = 0.0
        
    try:
        cfg.par.gas_logz
    except:
        cfg.par.gas_logz = 0

    try:
        cfg.par.FORCE_gas_logz
    except:
        cfg.par.FORCE_gas_logz = False

    try:
        cfg.par.source_logq
    except:
        cfg.par.source_logq = 1.e47

    try:
        cfg.par.FORCE_logq
    except:
        cfg.par.FORCE_logq = False

    try:
        cfg.par.FORCE_inner_radius
    except:
        cfg.par.FORCE_inner_radius = False
    
    try:
        cfg.par.inner_radius
    except:
        cfg.par.inner_radius = 1.e19

    try:
        cfg.par.use_Q
    except:
        cfg.par.use_Q = True

    try:
        cfg.par.neb_dust
    except:
        cfg.par.neb_dust = False

    try:
        cfg.par.cmdf_min_mass
    except:
        cfg.par.cmdf_min_mass = 3.5
    
    try:
        cfg.par.cmdf_max_mass
    except:
        cfg.par.cmdf_max_mass = 5.0

    try:
        cfg.par.cmdf_bins
    except:
        cfg.par.cmdf_bins = 6
    
    try:
        cfg.par.cmdf_beta

    except:
        cfg.par.cmdf_beta = -2.0

    try:
        cfg.par.SUBLIMATION
    except:
        cfg.par.SUBLIMATION = False

    try:
        cfg.par.SUBLIMATION_TEMPERATURE
    except:
        cfg.par.SUBLIMATION_TEMPERATURE = 1600

    try:
        cfg.model.TCMB
    except:
        cfg.model.TCMB = 2.73

    try:
        cfg.model.THETA
    except:
        cfg.model.THETA = 0

    try:
        cfg.model.PHI
    except:
        cfg.model.PHI = 0

    try:
        cfg.par.MANUAL_ORIENTATION
    except:
        cfg.par.MANUAL_ORIENTATION = False

    try:
        cfg.par.solar
    except:
        cfg.par.solar = 0.013

    try:
        cfg.par.dust_grid_type
    except:
        cfg.par.dust_grid_type = 'dtm'

    try:
        cfg.par.BH_model
    except:
        cfg.par.BH_model = None

    try:
        cfg.par.BH_modelfile
    except:
        cfg.par.BH_modelfile = None

    try:
        cfg.par.BH_var
    except:
        cfg.par.BH_var = False

    try:
        cfg.par.FORCE_STELLAR_AGES
    except:
        cfg.par.FORCE_STELLAR_AGES = False

    try:
        cfg.par.FORCE_STELLAR_AGES_VALUE
    except:
        cfg.par.FORCE_STELLAR_AGES_VALUE = 1

    try:
        cfg.par.FORCE_STELLAR_METALLICITIES
    except:
        cfg.par.FORCE_STELLAR_METALLICITIES = False

    try:
        cfg.par.FORCE_STELLAR_METALLICITIES_VALUE
    except:
        cfg.par.FORCE_STELLAR_METALLICITIES_VALUE = 0.013

    try:
        cfg.par.HII_T
    except:
        cfg.par.HII_T = 1.e4

    try: 
        cfg.par.HII_nh
    except:
        cfg.par.HII_nh = 1.e2

    try:
        cfg.par.HII_max_age
    except: 
        cfg.par.HII_max_age = 1.e-2
        
    try:
        cfg.par.HII_escape_fraction
    except:
        cfg.par.HII_escape_fraction = 0.0

    try:
        cfg.par.neb_abund
    except:
        cfg.par.neb_abund = "dopita"
    
    try:
        cfg.par.use_cloudy_tables
    except:
        cfg.par.use_cloudy_tables = False

    try:
        cfg.par.cloudy_cleanup
    except:
        cfg.par.cloudy_cleanup = True
    
    try:
        cfg.par.neb_file_output
    except:
        cfg.par.neb_file_output = True

    try:
        cfg.par.stellar_cluster_mass
    except:
        cfg.par.stellar_cluster_mass = 1.e4

    try:
        cfg.par.filterdir
    except:
        try:
            cfg.par.filter_file
            cfg.par.filterdir = os.path.dirname(cfg.par.filter_file)+'/'
        except:
            cfg.par.filterdir = ''

    try:
        cfg.par.filterfiles
    except:
        cfg.par.filterfiles = ['pdfilters.dat']

    try:
        cfg.par.PAH_frac
    except:
        cfg.par.PAH_frac = {'usg': 0.0586, 'vsg': 0.1351, 'big': 0.8063}

    return cfg.par.FORCE_RANDOM_SEED,cfg.par.BH_SED,cfg.par.IMAGING,cfg.par.SED,cfg.par.IMAGING_TRANSMISSION_FILTER,cfg.par.SED_MONOCHROMATIC,cfg.par.SKIP_RT,cfg.par.FIX_SED_MONOCHROMATIC_WAVELENGTHS,cfg.par.n_MPI_processes,cfg.par.SOURCES_RANDOM_POSITIONS,cfg.par.FORCE_gas_logu,cfg.par.gas_logu,cfg.par.gas_logu_init,cfg.par.gas_logz,cfg.par.FORCE_gas_logz,cfg.par.source_logq,cfg.par.FORCE_logq,cfg.par.FORCE_inner_radius,cfg.par.inner_radius,cfg.par.use_Q,cfg.par.neb_dust,cfg.par.cmdf_min_mass,cfg.par.cmdf_max_mass,cfg.par.cmdf_bins,cfg.par.cmdf_beta,cfg.par.SUBLIMATION,cfg.par.SUBLIMATION_TEMPERATURE,cfg.model.TCMB,cfg.model.THETA,cfg.model.PHI,cfg.par.MANUAL_ORIENTATION,cfg.par.solar,cfg.par.dust_grid_type,cfg.par.BH_model,cfg.par.BH_modelfile,cfg.par.BH_var,cfg.par.FORCE_STELLAR_AGES,cfg.par.FORCE_STELLAR_AGES_VALUE,cfg.par.FORCE_STELLAR_METALLICITIES,cfg.par.FORCE_STELLAR_METALLICITIES_VALUE,cfg.par.HII_T,cfg.par.HII_nh,cfg.par.HII_max_age,cfg.par.HII_escape_fraction,cfg.par.neb_abund,cfg.par.use_cloudy_tables,cfg.par.cloudy_cleanup,cfg.par.neb_file_output,cfg.par.stellar_cluster_mass, cfg.par.filterdir, cfg.par.filterfiles, cfg.par.PAH_frac
