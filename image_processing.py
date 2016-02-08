import config as cfg
import numpy as np
from astropy import units as u
import ipdb

def add_transmission_filters(image):

    for i in range(len(cfg.par.filter_list)):
        lam,throughput = np.loadtxt(cfg.par.filter_list[i],unpack=True)
        f = image.add_filter()
        f.name = 'dum'
        lam /= (1.+cfg.par.TRANSMISSION_FILTER_REDSHIFT)
        f.spectral_coord = lam * u.micron
        f.transmission = throughput * u.percent
        f.detector_type = 'energy'
        f.alpha = 1
        f.central_spectral_coord = np.mean(lam)*u.micron

    return None
