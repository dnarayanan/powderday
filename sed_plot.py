import numpy as np
import pdb
import constants as const
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput

m = ModelOutput('example.rtout')

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
wav,nufnu = m.get_sed(inclination = 'all',aperture = -1, distance = 300.*const.pc)


for i in range(nufnu.shape[0]):
    ax.loglog(wav,nufnu[i,:])

ax.set_xlabel(r'$\lambda$ [$\mu$m]')
ax.set_ylabel(r'$\lambda F_\lambda$ [erg/s/cm$^2$]')

#set axis limits

ax.set_xlim(0.1,5000)


fig.savefig('/Users/desika/Dropbox/junk.png')

