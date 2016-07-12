import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
from astropy.cosmology import Planck13
from astropy import units as u
from astropy import constants

#========================================================
#MODIFIABLE HEADER (make this a function later with argv)
z = 3
run = '/Users/desika/Dropbox/powderday/examples/gadget/example.200.rtout.sed'
#========================================================







fig = plt.figure()
ax = fig.add_subplot(1,1,1)


m = ModelOutput(run)
wav,flux = m.get_sed(inclination='all',aperture=-1)

wav  = np.asarray(wav)*u.micron #wav is in micron
wav *= (1.+z)

flux = np.asarray(flux)*u.erg/u.s
dl = Planck13.luminosity_distance(z)
dl = dl.to(u.cm)
    
flux /= (4.*3.14*dl**2.)
    
nu = constants.c.cgs/(wav.to(u.cm))
nu = nu.to(u.Hz)

flux /= nu
flux = flux.to(u.mJy)


for i in range(flux.shape[0]):
    ax.loglog(wav,flux[i,:])

ax.set_xlabel(r'$\lambda$ [$\mu$m]')
ax.set_ylabel('Flux (mJy)')
ax.set_ylim([1e-6,1.])
ax.set_xlim(0.05,15000)
ax.grid()

fig.savefig('/Users/desika/Dropbox/powderday/examples/gadget/sed.png')





