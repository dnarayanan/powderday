from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
from astropy.cosmology import Planck13
from astropy import units as u
from astropy import constants
import h5py

import fsps

#SKIRT STUFF
sedfile = '/home/desika.narayanan/SKIRT/run/pd_test_i90_sed.dat'
run = '/ufrc/narayanan/desika.narayanan/pd_dev_fork/examples/gadget/mw_zoom/example.134.rtout.sed'

data = np.loadtxt(sedfile)
skirt_lam = data[:,0]
skirt_flambda = data[:,1] #presumably in erg/s/cm^2/Angstrom

#PD stuff
m = ModelOutput(run)
pd_wav,pd_flux = m.get_sed(inclination='all',aperture=-1,units='ergs/s')
distance = 1.*u.Mpc.to(u.cm)
pd_flux *= u.erg/u.s
pd_flux/=(4.*np.pi*distance**2.)
pd_wav *= u.micron
pd_flux/=pd_wav.to(u.angstrom)

'''
#plot an fsps sed for 50 Myr and 0.012 metallicity
sp = fsps.StellarPopulation(sfh=0,logzsol=0.0)
wave,spec = sp.get_spectrum(tage = 0.05)
spec *= u.Lsun/u.Hz
wave *= u.angstrom
wave = wave.to(u.micron)
nu = const.c/wave
spec *= nu
spec/=wave
spec = spec.to(u.erg/u.s/u.angstrom)
spec/=(4.*np.pi*distance**2.)
'''

fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(skirt_lam,skirt_flambda,label='skirt')
ax.loglog(pd_wav,pd_flux[0,:],label='pd')
plt.legend()
ax.set_ylim([1.e-15,1.e0])
ax.set_xlabel(r'Wavelength ($\mu$m)')
ax.set_ylabel(r'f$_\lambda$ (cgs)')
plt.savefig('junk.png',dpi=300)

