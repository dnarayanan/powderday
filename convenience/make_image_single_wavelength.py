import numpy as np
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
import astropy.units as u

# ------------------------
# modifiable header
# ------------------------

m = ModelOutput(
    '/Users/desika/Dropbox/powderday/verification/gadget/example.200.rtout.image')
wav = 200  # micron

# ------------------------

# Extract the image for the first inclination, and scale to 300pc. We
# have to specify group=1 as there is no image in group 0.
image = m.get_image(units='mJy')

# Open figure and create axes
fig = plt.figure()
ax = fig.add_subplot(111)

# Find the closest wavelength
iwav = np.argmin(np.abs(wav - image.wav))

# Calculate the image width in kpc
w = image.x_max * u.cm
w = w.to(u.kpc)

# plot the beast
cax = ax.imshow(np.log(image.val[0, :, :, iwav]), cmap=plt.cm.magma,
                origin='lower', extent=[-w.value, w.value, -w.value, w.value])


# Finalize the plot
ax.tick_params(axis='both', which='major', labelsize=10)
ax.set_xlabel('x (kpc)')
ax.set_xlabel('y (kpc)')

plt.colorbar(cax, label='log Luminosity (ergs/s)', format='%.0e')

fig.savefig('pd_image.png', bbox_inches='tight', dpi=150)
