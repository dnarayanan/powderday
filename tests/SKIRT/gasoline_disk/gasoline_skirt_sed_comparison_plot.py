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
sedfile = '/ufrc/narayanan/desika.narayanan/SKIRT/run/tipsy_test/changa_test.dust_i90_sed.dat'
run = '/ufrc/narayanan/desika.narayanan/pd_git/tests/SKIRT/changa_disk/pd_skirt_comparison.changa.rtout.sed'

data = np.loadtxt(sedfile)
skirt_lam = data[:,0]
skirt_flambda = data[:,1]*u.W/u.m**2/u.micron
skirt_flambda = skirt_flambda.to(u.erg/u.s/u.cm**2./u.angstrom)

#PD stuff
m = ModelOutput(run)
pd_wav,pd_flux = m.get_sed(inclination='all',aperture=-1,units='ergs/s')
distance = 1.*u.Mpc.to(u.cm)
pd_flux *= u.erg/u.s
pd_flux/=(4.*np.pi*distance**2.)  #SKIRT seems to say F = L/d^2 instead of F = L/(4*pi*d^2)
pd_wav *= u.micron
pd_flux/=pd_wav.to(u.angstrom)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(pd_wav,pd_flux[0,:],label='powderday',lw=3)
ax.loglog(skirt_lam,skirt_flambda,label='SKIRT')
plt.legend()
ax.set_xlim([0.05,1000])
ax.set_ylim([1.e-17,1.e-8])
ax.set_xlabel(r'Wavelength ($\mu$m)')
ax.set_ylabel(r'f$_\lambda$ (erg s$^{-1}$ cm$^{-1}$ $\AA^{-1}$)')
plt.savefig('powderday_skirt_comparison.png',dpi=300)

