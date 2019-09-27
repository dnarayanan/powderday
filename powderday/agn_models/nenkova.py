from __future__ import print_function
import numpy as np
import h5py
import powderday.config as cfg

'''
    agn_spectrum using Nenkova+ (2008) torus models. Model calculations from CLUMPY (https://www.clumpy.org). Total spectrum is calculated as (Torus Flux) + (probability of AGN photon escape)x(AGN Flux). AGN spectra are assumed to be piecewise power-law (Rowan-Robinson 1995) with spectral breaks from Nenkova et al. Default CLUMPY model parameters taken from Nenkova et al.

    CLUMPY returns spectra in lambda * Flambda (arbitrary units). We scale such that the integrated Flambda gives the total IR luminosity.

    - Ray Sharma
'''


class Nenkova2008:
    def __init__(self, params=[5, 30, 0, 1.5, 30, 40]):
        N0, Y, i, q, sig, tv = params
        self.N0 = N0
        self.Y = Y
        self.i = i
        self.q = q
        self.sig = sig
        self.tv = tv

    def agn_spectrum(self, log_L_bol):
        try:
            h = h5py.File(cfg.par.BH_modelfile, 'r')
        except:
            raise IOError('Unable to find Nenkova BH model file. '
                                    'Check the path in parameters master, or '
                                    'download the file here: https://www.clump'
                                    'y.org/downloads/clumpy_models_201410_tvav'
                                    'g.hdf5')

        ix = ((h['N0'][:] == self.N0) &
              (h['Y'][:] == self.Y) &
              (h['i'][:] == self.i) &
              (h['q'][:] == self.q) &
              (h['sig'][:] == self.sig) &
              (h['tv'][:] == self.tv))

        nu_vec = 3e14 / h['wave'][:]

        frac_AGN_obsc = h['ptype1'][:][ix][0]
        l_band_vec_torus = h['flux_tor'][:][ix][0]
        l_band_vec_AGN = h['flux_toragn'][:][ix][0] - l_band_vec_torus
        l_band_vec = l_band_vec_torus + (frac_AGN_obsc * l_band_vec_AGN)

        l_band_vec = self.scale_spectrum(l_band_vec, nu_vec, log_L_bol)

        l_band_vec = np.log10(l_band_vec)
        l_band_vec = np.concatenate((l_band_vec, [0, 0, 0, 0]))
        nu_vec = np.log10(nu_vec)
        nu_vec = np.concatenate((nu_vec, [-1, -2, -3, -4]))

        to_cgs = np.log10(3.9) + 33
        return nu_vec, l_band_vec + to_cgs

    def scale_spectrum(self, l_band_vec, nu_vec, log_L_bol):
        ''' Scale the spectrum by (total IR luminosity) / (integrated spectrum in arb. units)
        '''
        L_IR = 10**self.bol_correct_IR(log_L_bol)
        integrated_spec = np.trapz(l_band_vec / nu_vec, nu_vec)
        norm = L_IR / abs(integrated_spec)
        return l_band_vec * norm

    def bol_correct_IR(self, log_L_bol, c1=17.87, k1=0.28, c2=10.03, k2=0.020):
        ''' Return log IR luminosity using bolometric corrections from Hopkins+ (2006). Defaults to 15micron band corrections.
        '''
        L_bol = 10**log_L_bol
        L_IR = L_bol / (c1 * pow(L_bol / 1e10, k1) +
                        c2 * pow(L_bol / 1e10, k2))
        return np.log10(L_IR)
