import yt
import numpy as np
import pdb

snapshot = '/blue/narayanan/desika.narayanan/yt_datasets/FIRE_M12i_ref11/snapshot_600.hdf5'
snapshot = '/orange/narayanan/pg3552/output_m12n256_boxes/output_m12n256_0/snapshot_064.hdf5'

def _size_with_units(field,data):
    return data.ds.parameters['size']

def _position_x(field,data):
    return data[('PartType3', 'Coordinates')][:,0]


ds = yt.load(snapshot)
ds.add_field('posx',function=_position_x,sampling_type='particle',units='cm',particle_type=True)
ad = ds.all_data()


#loop through the sizes now 
nsizes = ad['PartType3','Dust_Size'].shape[1]
#set up master array to hold the deposited sizes


for isize in range(2):

    ds.parameters['size'] = 0 #just to clear it out

    #we have to slice it before we add the field.  if you try to slice
    #in the field defintion (i.e., in _size_with_units), yt freaks
    ds.parameters['size'] = ad['PartType3','Dust_Size'][:,isize]

    #actually add the sliced field now.  we do this so that we can
    #prepare for depositing onto the octree. we call this a dummy size
    #since this is just there for a place holder to deposit into a dummy octree
    ds.add_field(('PartType3','dummy_size'),function=_size_with_units,sampling_type='particle',units='dimensionless',particle_type=True,force_override=True)

    pdb.set_trace()

    #deposit onto the octree.   this is what will get merged into the final octree
    ds.add_deposited_particle_field(('PartType3','dummy_size'),"sum")

    if isize == 0:
        octree_of_sizes = np.zeros((ad[('deposit','PartType3_sum_dummy_size')].shape[0],nsizes))
        
    octree_of_sizes[:,isize] = ad[('deposit','PartType3_sum_dummy_size')]
    



