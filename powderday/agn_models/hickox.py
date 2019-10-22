from __future__ import print_function
import numpy as np
'''
    Vary the instantaneous BH luminosity by sampling the probability distribution from Hickox et al (2014) for short time-scale variations.
    - Ray Sharma
'''
def Hickox2014(L_cut=100, alpha=0.2):
    L = np.logspace(-5, 3, 100)
    p0 = pow(L / L_cut, -alpha) * np.exp(-L / L_cut)
    t0 = 0.00854
    return t0 * p0, L

def vary_bhluminosity(size):
    '''
        Returns a list of fractions by which to vary the luminosity.
    '''
    PDF, L_frac = Hickox2014()
    CDF = np.cumsum(PDF) / sum(PDF)
    choice = np.random.random(size)
    ix = [np.argmin(abs(CDF - x)) for x in choice]
    return L_frac[ix]
