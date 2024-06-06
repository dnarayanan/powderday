#plan:

#1. create an MRN distribution and plot its extinction law and kappa
#2. compare this to the WD01
#3. run a pd run and test it
#4. make albedos identical to one another and see what impact is
#4. now compare the mrn_test to the active dust model

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dust_file_writer import *
import astropy.units as u
import astropy.constants as constants
from hyperion.dust import IsotropicDust
import h5py


#pull in the WD01 kappa to compare to ours
f = h5py.File('/home/desika.narayanan/hyperion-dust-0.1.0/dust_files/d03_5.5_3.0_A.hdf5','r')
nu_wd = f['optical_properties']['nu']*u.Hz
chi_wd = f['optical_properties']['chi']
albedo_wd = f['optical_properties']['albedo']


loga = np.linspace(-4,0,16)
a = 10.**loga
mrn_dn_da = a**(-3.5)
mrn_dn_da /= np.max(mrn_dn_da)

da = [a[i+1]-a[i] for i in range(len(a)-1)]
da.append(da[-1])
mrn_dn = mrn_dn_da*da

np.savetxt('mrn_dn.txt',mrn_dn)


#wlen = 1. / np.logspace(-4,3,201)*u.micron
#nu = (constants.c/wlen).to(u.Hz)
wlen = (constants.c/nu_wd).to(u.micron)
nu = nu_wd

cfrac = 0.54

#compute quantities for exmaple dsf above
xtab, Qtab = Qext_tab_load()
t_Qext = Qext_get(loga,wlen.value,cfrac,xtab,Qtab)
t_Qext_V = Qext_get(loga,np.array([0.551]),cfrac,xtab,Qtab)
Aext, albedo, Qext = extinction_law(a,mrn_dn,wlen.value,cfrac,t_Qext,t_Qext_V)

#calculate kappa
#kappa (a) = 3 * Qext (a) / ( 4 * a * rho_grain)
ASSUMED_DENSITY_OF_DUST = 2.4*u.g/u.cm**3
kappa_a_lambda = np.zeros(Qext.shape)
for i in range(kappa_a_lambda.shape[1]):
    kappa_a_lambda[:,i] = 3.*Qext[:,i]/(a*u.micron.to(u.cm)*ASSUMED_DENSITY_OF_DUST)
kappa_lambda = np.sum(kappa_a_lambda,axis=0)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(1./wlen,Aext)
ax.set_xlim([1,10])
ax.set_ylim([1.e-8,1.e-5])
fig.savefig('extinction.png',dpi=300)



fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(nu_wd,chi_wd/np.max(chi_wd),label='wd01')
ax.loglog(nu,kappa_lambda/np.max(kappa_lambda),label='mrn')
ax.set_xlabel(r'Frequency (Hz)')
ax.set_ylim([1e-6,1])
ax.set_ylabel(r'$\chi$')
plt.legend(loc=2)
fig.savefig('chi.png',dpi=300)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.loglog(nu_wd,albedo_wd/np.max(albedo_wd),label='wd01')
ax.loglog(nu,albedo/np.max(albedo),label='mrn')
ax.set_xlabel(r'Frequency (Hz)')
ax.set_ylim([1e-6,1])
ax.set_ylabel(r'$\chi$')
plt.legend(loc=2)
fig.savefig('albedo.png',dpi=300)


#write dust files
from hyperion.dust import IsotropicDust
d = IsotropicDust(nu.value,albedo,kappa_lambda)
d.write('mrn_test.hdf5')
