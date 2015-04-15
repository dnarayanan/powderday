import numpy as np
import astropy.units as u
import yt
from yt import derived_field
import ipdb
import config as cfg
from astropy.cosmology import Planck13
import astropy.units as u


def ramses_field_add(fname,bounding_box = None,ds = None,starages = False):
    
    def _starmetals(field,data):
        return data[('io', 'particle_metallicity')]
    if fname != None:
        ds = yt.load(fname,bounding_box=bounding_box,over_refine_factor=cfg.par.oref,n_ref=cfg.par.n_ref)
        ds.index

        
        ad = ds.all_data()
    
        return ds

