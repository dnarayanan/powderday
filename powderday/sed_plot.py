import numpy as np
import astropy.constants as constants
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput

m = ModelOutput('example.rtout')

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
wav,nufnu = m.get_sed(inclination = 'all',aperture = -1)


for i in range(nufnu.shape[0]):
    ax.loglog(wav,nufnu[i,:])

ax.set_xlabel(r'$\lambda$ [$\mu$m]')
ax.set_ylabel(r'$\lambda L_\lambda$ [erg/s/cm$^2$]')

#set axis limits

ax.set_xlim(0.05,1000)
#ax.set_ylim(1e35,1e43)
#ax.set_ylim(1e-14,2e-6)
#ax.set_ylim(1e-10,1e-5)

fig.savefig('junk.png')



def stellar_seds_view():

    data = np.load('stellar_seds.npz')
    stellar_nu = data['arr_0']
    stellar_fnu = data['arr_1']
    disk_fnu = data['arr_2']
    bulge_fnu = data['arr_3']

    nstars = stellar_fnu.shape[0]
    
    #get lambda from nu
    lambda_cgs = constants.c.cgs.value/stellar_nu
    lambda_micron = lambda_cgs*1.e4
    

    stellar_lum_lsun = stellar_fnu
    for i in range(stellar_fnu.shape[0]):
        stellar_lum_lsun[i,:] = stellar_fnu[i,:]*stellar_nu


    for i in range(nstars):
       
        if i == 0: ax.plot(lambda_micron[:],stellar_fnu[i,:],color='grey')
        if i > 0:  ax.plot(lambda_micron[:],stellar_fnu[i,:],color='grey')
        
        
        
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel('micron')
        ax.set_ylabel('Lsun/Hz')
        
      
