from __future__ import print_function
import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../')
from agn_spectrum import *

from astropy import units as u
from astropy import constants as const

from hyperion.model import ModelOutput

fig = plt.figure()
ax = fig.add_subplot(111)

log_lum = np.linspace(9,13,100)*u.Lsun

norm = matplotlib.colors.Normalize(
    vmin = np.min(log_lum.value),
    vmax = np.max(log_lum.value))
c_m = matplotlib.cm.viridis_r
s_m  = matplotlib.cm.ScalarMappable(cmap = c_m,norm=norm)
s_m.set_array([])


for lum in log_lum:
    print('lum = %e'%lum.value)
    nu,bh_fnu = agn_spectrum(lum.value)

    #regularlizing units
    bh_fnu = bh_fnu[0:-4]
    bh_fnu = 10.**bh_fnu * u.erg/u.s
    bh_fnu = bh_fnu.to(u.Lsun)

    nu = nu[0:-4]
    nu = 10.**nu
    nu *= u.Hz

    lam = (const.c/nu).to(u.micron)
    ax.loglog(lam,bh_fnu,color=s_m.to_rgba(lum.value))
    


#load in the bh sed from a powderday run
data = np.load('/ufrc/narayanan/desika.narayanan/pd_runs/ena/bh_sed.npz')
nholes = data['luminosity'].shape[0]
for i in range(nholes):
    pd_nu =  data['nu']*u.Hz
    pd_lam = (const.c/pd_nu).to(u.micron)
    pd_fnu  = (data['fnu'][i][:]*u.erg/u.s).to(u.Lsun).value
    if data['luminosity'][i] > 0:
        ax.plot(pd_lam.value,pd_fnu)

#now plot the powderday SED
run = '/ufrc/narayanan/desika.narayanan/pd_runs/ena/example.094.rtout.bhon.sed'
m = ModelOutput(run)
wav,flux = m.get_sed(inclination='all',aperture=-1)
fullrun_wav = np.asarray(wav)*u.micron
fullrun_flux = np.asarray(flux)*u.erg/u.s
fullrun_nu = (const.c/fullrun_wav).to(u.Hz)
fullrun_fnu = fullrun_flux/fullrun_nu

ax.plot(fullrun_wav.value,fullrun_fnu[0,:].value/1.e20)

ax.set_xlabel(r'Wavelength ($\mu$m)')
ax.set_ylabel(r'F$_\nu$')
cb = fig.colorbar(s_m,orientation='vertical')
cb.set_label(r'Black Hole L$_\mathrm{bol}$ (L$_\odot$)')
cb.ax.tick_params(labelsize=8)




fig.savefig('junk.png',dpi=300)
