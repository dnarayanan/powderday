#Script showing how to extract some rtout files that were run on an
#octree format

from __future__ import print_function
from hyperion.model import ModelOutput
from hyperion.grid.yt3_wrappers import find_order
import astropy.units as u
import numpy as np

run = '/home/desika.narayanan/pd_git/tests/SKIRT/gizmo_mw_zoom/pd_skirt_comparison.134.rtout.sed'

m = ModelOutput(run)


oct = m.get_quantities()
#ds = oct.to_yt()

#ripped from hyperion/grid/yt3_wrappers.py -- we do this because
#something about load_octree in yt4.x is only returning the first cell
grid = oct
order = find_order(grid.refined)
refined = grid.refined[order]

quantities = {}
for field in grid.quantities:
    quantities[('gas', field)] = np.atleast_2d(grid.quantities[field][0][order][~refined]).transpose()

specific_energy = quantities['gas','specific_energy']*u.erg/u.s/u.g
dust_temp = quantities['gas','temperature']*u.K
dust_density = quantities['gas','density']*u.g/u.cm**3
