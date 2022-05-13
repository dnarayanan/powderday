from astropy import constants
import itertools
import numpy as np
from scipy import integrate, spatial
import powderday.config as cfg

"""
----------------------------------------------------------------------------------------------------------------
From cloudyfsps written by Nell Byler.
(Source https://github.com/nell-byler/cloudyfsps/blob/master/cloudyfsps/generalTools.py
 retrieved in October 2019)
----------------------------------------------------------------------------------------------------------------
"""


def calc_LogQ(nuin0, specin0, efrac=0.0, mstar=1.0, mfrac=1.0):
    '''
    Claculates the number of lyman ionizing photons (Q) for given a spectrum
    Q = int(Lnu/hnu dnu, nu_0, inf)
    Here:
        specin0: Input spectrum must be in ergs/s/Hz
        nuin0: Freq of the input spectrum must be in Hz
        mstar: Mass of the particle in units of Msun
        efrac: The escape fraction of H-ionizing photons
    '''

    c = constants.c.cgs.value  # cm/s
    h = constants.h.cgs.value  # erg/s
    lam_0 = 911.6 * 1e-8  # Halpha wavelength in cm

    nuin = np.asarray(nuin0)
    specin = np.asarray(specin0)
    nu_0 = c / lam_0
    inds, = np.where(nuin >= nu_0)
    hlam, hflu = nuin[inds], specin[inds]
    nu = hlam[::-1]
    f_nu = hflu[::-1]
    integrand = f_nu / (h * nu)
    logQ = np.log10(integrate.simps(integrand, x=nu)*(mstar/mfrac)*(1-efrac))
    return logQ

   
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


def cmdf(stellar_mass, nbins, min_mass, max_mass, beta, rescale_masses=True):
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
        num.append(int(round(N)))

    ## rescale masses so that they sum to initial mass
    if rescale_masses:
        _mass = 10**np.array(mass)
        _scale = stellar_mass / np.sum(_mass * num)
        mass = np.log10(_mass * _scale)

    return mass, num


def convert_metals(metals):
    """
    Converts metalicity from units of percentage by mass (SIMBA) 
    to atom per hydrogen atoms (CLOUDY)
    """
    per_H = 0.7381 # Photospheric mass fraction from Asplund et al. 2009
    # mass of elements in unified atomic mass units
    # [He, C, N, O, Ne, Mg, Si, S, Ca, Fe]
    mass_H = 1.008
    mass = [4.002602, 12.001, 14.007, 15.999, 20.1797,
            24.305, 28.085, 32.06, 40.078, 55.845]
    metals_conv = np.zeros(len(metals))
    # Converting from mass fraction to atomic fraction
    for i in range(len(metals)):
        metals_conv[i] = np.log10((metals[i]/per_H)*(mass_H/mass[i]))

    return metals_conv


def get_nearest(particle_list, particle_central, num=32):
    
    all_particles = particle_list.copy()
    all_particles = np.array(all_particles)
    tree = spatial.KDTree(all_particles)
    arg = (tree.query(particle_central, k=num))
    # Removing
    mask = np.where(arg[1]<len(all_particles))[0]
    
    return np.array(arg[0][mask]), np.array(arg[1][mask])


def age_dist(Num, tavg, width=5, gamma=-0.65, bins=4, tries=100, tolerance = 0.1):
    # If the number of star particles is less than the number of bins or if the age is not within
    # the limits set in paramters_master then do not break the particle down

    if Num <= bins or tavg <= cfg.par.age_dist_min or tavg >= cfg.par.age_dist_max:
        return np.array([Num]), np.array([tavg])

    tavg = tavg*1e3

    # Solving for the limits of the age distribution (dN/dt is proportional to t^gamma) such that the
    # the average age of the distribution is as close to the age of the parent star particle. 
    # If the difference between the average age and the age of the parent star particle is above the threshold 
    # then the width (difference between start and end age of the distribution) is decreased at each step
    # until the threshold is met (This is done at most x times where is determined by the variable "tries". 
    # If the max number of tries is reached and the threshold is not met then the particle is not broken down.

    for i in range(tries):
        ti_temp = np.arange(1, cfg.par.HII_max_age*1e3, 0.1)
        x = (ti_temp + width) ** (gamma + 2) + ti_temp ** (gamma + 2) - 2 * (tavg ** (gamma + 2))
        index = np.argmin(np.abs(x))
        ti = ti_temp[index]
        tf = ti_temp[index] + width

        t = np.linspace(ti, tf, bins)
        A = Num / (np.sum(t ** (gamma + 1)))
        N = [int(k) for k in A * t ** (gamma + 1)]
        t_avg = np.sum(N * t) / np.sum(N)

        if np.round(width, 1) < 0:
            return np.array([Num]), np.array([tavg])

        elif np.abs(t_avg - tavg) <= tolerance:
            break

        else:
            width = width - 0.2

    return N, t*1e-3
