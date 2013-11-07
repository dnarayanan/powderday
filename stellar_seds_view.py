import numpy as np
import pdb
import constants as const
import matplotlib.pyplot as plt

data = np.load('stellar_seds.npz')
stellar_nu = data['arr_0']
stellar_fnu = data['arr_1']
disk_fnu = data['arr_2']
bulge_fnu = data['arr_3']

nstars = stellar_fnu.shape[0]

#get lambda from nu
lambda_cgs = const.c/stellar_nu
lambda_micron = lambda_cgs*1.e4


fig = plt.figure()
ax = fig.add_subplot(111)

for i in range(nstars):
    if i == 0: ax.plot(lambda_micron[:],stellar_fnu[i,:],label='newstars')
    if i > 0:  ax.plot(lambda_micron[:],stellar_fnu[i,:])
ax.plot(lambda_micron[:],disk_fnu[:],color='red',lw=10,label='disk')
ax.plot(lambda_micron[:],bulge_fnu[:],color='blue',lw=10,label='bulge')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('micron')
ax.set_ylabel('Lsun/Hz')

plt.legend()
fig.savefig('stellar_seds.png')


