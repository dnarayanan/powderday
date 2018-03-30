import numpy as np
import config as cfg


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
        cfg.par.solar
    except:
        cfg.par.solar = 0.013

    try:
        cfg.par.dust_grid_type
    except:
        cfg.par.dust_grid_type = 'dtm'

    return cfg.par.FORCE_RANDOM_SEED,cfg.par.BH_SED,cfg.par.IMAGING,cfg.par.SED,cfg.par.IMAGING_TRANSMISSION_FILTER,cfg.par.SED_MONOCHROMATIC,cfg.par.SKIP_RT,cfg.par.FIX_SED_MONOCHROMATIC_WAVELENGTHS,cfg.par.n_MPI_processes,cfg.par.SOURCES_RANDOM_POSITIONS,cfg.par.gas_logu,cfg.par.gas_logz,cfg.par.FORCE_gas_logz,cfg.par.SUBLIMATION,cfg.par.SUBLIMATION_TEMPERATURE,cfg.model.TCMB,cfg.par.solar,cfg.par.dust_grid_type
