from __future__ import print_function
import numpy as np
import h5py
import powderday.config as cfg
from powderday.agn_models.hopkins import agn_spectrum as underlying_agn_spectrum

'''
    agn_spectrum using Nenkova+ (2008) torus models. Model calculations from CLUMPY (https://www.clumpy.org). 
    Total spectrum is calculated as (torus flux ratio) x (underlying AGN spectra). 
    Underlying AGN spectra are assumed to be Hopkins+ (2007) templates. Default CLUMPY model parameters taken from Nenkova et al.

    CLUMPY returns spectra in lambda * Flambda (arbitrary units).

    - Ray Sharma
'''


class Nenkova2008:
    def __init__(self, N0=5, Y=30, i=0, q=1.5, sig=30, tv=40):
        self.N0 = N0
        self.Y = Y
        self.i = i
        self.q = q
        self.sig = sig
        self.tv = tv
        try:
            self.h = h5py.File(cfg.par.BH_modelfile, 'r')
        except IOError:
            raise IOError('Unable to find Nenkova BH model file. '
                          'Check the path in parameters master, or '
                          'download the file here: https://www.clump'
                          'y.org/downloads/clumpy_models_201410_tvav'
                          'g.hdf5')
        self.check_params()

    def agn_spectrum(self, log_L_bol):
        nu_vec = 3e14 / self.h['wave'][:]
        nu_vec = np.log10(nu_vec)
        nu_vec = np.concatenate((nu_vec, [-1, -2, -3, -4]))

        l_band_vec_torus = self.h['flux_tor'][:][self.ix][0]

        agn_nu, agn_l_band_vec = underlying_agn_spectrum(log_L_bol)
        l_band_vec = np.log10(l_band_vec_torus) + np.interp(
            nu_vec[:-4], agn_nu[:-4], agn_l_band_vec[:-4])
        l_band_vec = np.concatenate((l_band_vec, [0, 0, 0, 0]))

        return nu_vec, l_band_vec

    def check_params(self):
        self.ix = ((self.h['N0'][:] == self.N0) & (self.h['Y'][:] == self.Y) &
                   (self.h['i'][:] == self.i) & (self.h['q'][:] == self.q) &
                   (self.h['sig'][:] == self.sig) & (self.h['tv'][:] == self.tv))
        if self.ix.sum() == 0:
            raise IOError('Input Nenkova model parameters cannot be found in model file.')
