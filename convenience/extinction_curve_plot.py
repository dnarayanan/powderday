from hyperion.dust import SphericalDust
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from astropy import units as u
from astropy import constants as constants

dustfile = '/ufrc/narayanan/desika.narayanan/pd/hyperion-dust/dust_files/d03_3.1_6.0_A.hdf5'
dust = SphericalDust(dustfile)
mw_df_nu = dust.optical_properties.nu
mw_df_chi  = dust.optical_properties.chi

mw_df_lam = (constants.c/(mw_df_nu*u.Hz)).to(u.micron)
x = 1./mw_df_lam


fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x,mw_df_chi)
ax.set_xlim([1,10])
ax.set_xlabel('x')
ax.set_ylabel(r'$\chi$')
fig.savefig('extinction.png',dpi=300)
