import hyperion
import yt
from hyperion.model import Model
from hyperion.grid import AMRGrid

snapshot = '//ufrc/narayanan/desika.narayanan/yt_datasets/enzo_cosmology_plus/RD0009/RD0009'

ds = yt.load(snapshot)

def _dust_density(field,data):
    return (0.4*data[('gas','metallicity')]*data[('gas','density')]).in_units('g/cm**3')

ds.add_field(('gas','dust_density'),function=_dust_density,units='g/cm**3')

m = Model()
#amr = AMRGrid.from_yt(ds,quantity_mapping={'density':('gas','density')})
amr = AMRGrid.from_yt(ds,quantity_mapping={'density':('gas','dust_density')})
m.set_amr_grid(amr)
