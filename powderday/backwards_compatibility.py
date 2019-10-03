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
        cfg.par.gas_logz
    except:
        cfg.par.gas_logz = 0

    try:
        cfg.par.FORCE_gas_logz
    except:
        cfg.par.FORCE_gas_logz = False

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
        FORCE_STELLAR_AGES
    except:
        FORCE_STELLAR_AGES = False

    try:
        FORCE_STELLAR_AGES_VALUE
    except:
        FORCE_STELLAR_AGES_VALUE = 1

    try:
        FORCE_STELLAR_METALLICITIES
    except:
        FORCE_STELLAR_METALLICITIES = False

    try:
        FORCE_STELLAR_METALLICITIES_VALUE
    except:
        FORCE_STELLAR_METALLICITIES_VALUE = 0.013

    try:
        HII_T
    except:
        HII_T = 1.e4

    try: 
        HII_nh
    except:
        HII_nh = 1.e2

    try:
        HII_max_age
    except: 
        HII_max_age = 2.e-3
        
    try:
        neb_file_output
    except:
        neb_file_output = True

    try:
        stellar_cluster_mass
    except:
        stellar_cluster_mass = 1.e4

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

    return cfg.par.FORCE_RANDOM_SEED,cfg.par.BH_SED,cfg.par.IMAGING,cfg.par.SED,cfg.par.IMAGING_TRANSMISSION_FILTER,cfg.par.SED_MONOCHROMATIC,cfg.par.SKIP_RT,cfg.par.FIX_SED_MONOCHROMATIC_WAVELENGTHS,cfg.par.n_MPI_processes,cfg.par.SOURCES_RANDOM_POSITIONS,cfg.par.FORCE_gas_logu,cfg.par.gas_logu,cfg.par.gas_logz,cfg.par.FORCE_gas_logz,cfg.par.SUBLIMATION,cfg.par.SUBLIMATION_TEMPERATURE,cfg.model.TCMB,cfg.model.THETA,cfg.model.PHI,cfg.par.MANUAL_ORIENTATION,cfg.par.solar,cfg.par.dust_grid_type,cfg.par.BH_model,cfg.par.BH_modelfile,cfg.par.BH_var,cfg.par.FORCE_STELLAR_AGES,cfg.par.FORCE_STELLAR_AGES_VALUE,cfg.par.FORCE_STELLAR_METALLICITIES,cfg.par.FORCE_STELLAR_METALLICITIES_VALUE,HII_T,HII_nh,HII_max_age,neb_file_output,stellar_cluster_mass, cfg.par.filterdir, cfg.par.filterfiles, cfg.par.PAH_frac
