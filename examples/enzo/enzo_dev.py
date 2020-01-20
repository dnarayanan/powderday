import yt
import numpy as np
import ipdb
from hyperion.model import Model
from hyperion.grid import AMRGrid
from hyperion.util.constants import pc, lsun

fname = "/Users/desika/Dropbox/yt_datasets/enzo_iso_galaxy/galaxy0030/galaxy0030"


ds = yt.load(fname)
ds.index

levels = ds.index.max_level

amr = AMRGrid()

'''
for ilevel in range(levels):
    level = amr.add_level()
    
    for igrid in ds.index.select_grids(ilevel):
        print igrid
        grid = level.add_grid()
        grid.xmin,grid.xmax = igrid.LeftEdge[0].in_units('cm'),igrid.RightEdge[0].in_units('cm')
        grid.ymin,grid.ymax = igrid.LeftEdge[1].in_units('cm'),igrid.RightEdge[1].in_units('cm')
        grid.zmin,grid.zmax = igrid.LeftEdge[2].in_units('cm'),igrid.RightEdge[2].in_units('cm')
        grid.quantities["density"] = np.transpose(np.array(igrid[("gas","metal_density")].in_units('g/cm**3')*0.4))
        
        grid.nx,grid.ny,grid.nz = igrid[("gas","metal_density")].shape

'''


ds = yt.load(fname)

def _dust_density(field, data):
    return data[('gas', 'metal_density')].in_units("g/cm**3")*0.4

ds.add_field(('gas', 'dust_density'), function=_dust_density, units='g/cm**3')

amr = AMRGrid.from_yt(ds, quantity_mapping={'density':('gas','dust_density')})

# Set up Hyperion model

import numpy as np

from hyperion.model import Model
from hyperion.util.constants import pc

m = Model()

m.set_amr_grid(amr)

m.add_density_grid(amr['density'], '/Users/desika//pd/hyperion-dust-0.1.0/dust_files/d03_3.1_6.0_A.hdf5')

m.set_minimum_temperature(10)

m.set_n_initial_iterations(0)

m.set_raytracing(True)

# Add a point source in the center
s = m.add_point_source()
s.position = (0, 0., 0.)
s.luminosity = 1000 * lsun
s.temperature = 6000.

i = m.add_peeled_images(sed=False, image=True)
i.set_image_limits(-1 * pc, 1 * pc, -1 * pc, 1 * pc)
i.set_image_size(256, 256)
i.set_viewing_angles(np.repeat(30, 8), np.linspace(0., 360, 9)[:-1])
i.set_wavelength_range(1, 1000, 1001)

m.set_n_photons(imaging=0, raytracing_sources=0, raytracing_dust=1e6)

m.set_copy_input(False)

m.write('test.rtin', overwrite=True)
m.run('test.rtout', overwrite=True)
