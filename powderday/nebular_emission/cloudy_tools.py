from astropy import constants
import itertools
import numpy as np
from scipy import integrate

"""
----------------------------------------------------------------------------------------------------------------
From cloudyfsps written by Nell Byler.
(Source https://github.com/nell-byler/cloudyfsps/blob/master/cloudyfsps/generalTools.py
 retrieved in October 2019)
----------------------------------------------------------------------------------------------------------------
"""


def calc_LogU(nuin0, specin0, nh, T, efrac, mstar=1.0):
    '''
    Claculates the number of lyman ionizing photons for given a spectrum
    Input spectrum must be in ergs/s/Hz!!
    Q = int(Lnu/hnu dnu, nu_0, inf) , number of hydrogen ionizing photons
    mstar is in units of solar mass
    Rin is in units of cm-3
    nh is in units of cm-3
    '''

    c = constants.c.cgs.value  # cm/s
    h = constants.h.cgs.value  # erg/s
    alpha = 2.5e-13*((T/(10**4))**(-0.85)) # cm3/s
    lam_0 = 911.6 * 1e-8  # Halpha wavelength in cm

    nuin = np.asarray(nuin0)
    specin = np.asarray(specin0)
    nu_0 = c / lam_0
    inds, = np.where(nuin >= nu_0)
    hlam, hflu = nuin[inds], specin[inds]
    nu = hlam[::-1]
    f_nu = hflu[::-1]
    integrand = f_nu / (h * nu)
    logQ = np.log10(integrate.simps(integrand, x=nu)*mstar*(1-efrac)) 
    Rin = (3 * (10 ** logQ) / (4 * np.pi * nh * nh * alpha)) ** (1. / 3.)
    logU = np.log10((10**logQ)/(4*np.pi*Rin*Rin*nh*c))
    return logQ, Rin, logU


def air_to_vac(inpt, no_uv_conv=True):
    """
    from morton 1991
    preserves order of input array
    """
    if type(inpt) is float:
        wl = np.array([inpt])
    else:
        wl = np.asarray(inpt)
    to_vac = lambda lam: (6.4328e-5 + (2.94981e-2/(146.0-(1.0e4/lam)**2.0)) + (2.554e-4/(41.0-(1.0e4/lam)**2.0)))*lam + lam
    if no_uv_conv:
        outpt = np.array([to_vac(lam) if lam > 2000.0 else lam for lam in wl])
    else:
        outpt = to_vac(wl)
    return outpt


def sym_to_name(val=None):
    elem_keys = dict(He="helium",
                     C="carbon",
                     N="nitrogen",
                     O="oxygen",
                     Ne="neon",
                     Mg="magnesium",
                     Si="silicon",
                     S="sulphur",
                     Ar="argon",
                     Ca="calcium",
                     Fe="iron",
                     F="fluorine",
                     Na="sodium",
                     Al="aluminum",
                     Cl="chlorine",
                     Ni="nickel",
                     P="phosphorus",
                     Sc="scandium",
                     K="potassium",
                     Ti="titanium",
                     V="vanadium",
                     Cr="chromium",
                     Co="cobalt",
                     Cu="copper",
                     Mn="manganese",
                     Zn="zinc")
    if val is None:
        return elem_keys
    else:
        try:
            return elem_keys[val.title()]
        except KeyError:
            print("element not in ", elem_keys.keys())


def grouper(n, iterable):
    """
    Iterate through array in groups of n
    """
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, n))
        if not chunk:
            return
        yield chunk


def cmdf(stellar_mass, nbins, min_mass, max_mass, beta):
    """
    Calculates the number of clusters per mass interval assuming a cluster
    mass distribution function of the form dN/dM goes as M^(-beta)
    """
    interval = (max_mass-min_mass)/nbins
    num = []
    mass = []
    for i in range(nbins):
        m = min_mass + (i*interval)
        mass.append(m)

    denom = sum((10**q)**(beta + 2) for q in mass)
    A = (stellar_mass / denom)

    for i in range(nbins):
        N = A*((10**mass[i])**(1. + beta))
        num.append(round(N))
    return mass, num


def convert_metals(metals):
    """
    Converts metalicity from units of percentage by mass (SIMBA) 
    to atom per hydrogen atoms (CLOUDY)
    """
    # mass of elements in unified atomic mass units
    # [He, C, N, O, Ne, Mg, Si, S, Ca, Fe]
    per_H = 0.7314
    mass_H = 1.008
    mass = [4.002602, 12.001, 14.007, 15.999, 20.1797, 
            24.305, 28.085, 32.06, 40.078, 55.845]
    metals_conv = np.zeros(len(metals))
    for i in range(len(metals)):
        metals_conv[i] = np.log10((metals[i]/per_H)*(mass_H/mass[i]))

    return metals_conv

