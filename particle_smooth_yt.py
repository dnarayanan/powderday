import numpy as np
from yt.mods import *
from yt.fields.particle_fields import add_volume_weighted_smoothed_field
import pdb

def yt_smooth(pf):

    print 'starting yt_smooth'
    pf.index
    ad = pf.all_data()
    
    #Metallicity_00
    '''
    metal_fn = add_volume_weighted_smoothed_field("PartType0", "Coordinates", "Masses",
                                                  "SmoothingLength", "Density","Metallicity_00",
                                                  pf.field_info)
    '''
    metal_fn = add_volume_weighted_smoothed_field("PartType0", "Coordinates", "Masses",
                                                  "SmoothingLength", "Density","Metallicity",
                                                  pf.field_info)

    mass_fn = add_volume_weighted_smoothed_field("PartType0", "Coordinates", "Masses",
                                                  "SmoothingLength", "Density","Masses", 
                                                  pf.field_info)
    
    pf.field_info.alias(("gas", "metallicity"), metal_fn[0])
    pf.field_info.alias(("gas", "masses"), mass_fn[0])
    

    saved = pf.index.oct_handler.save_octree(always_descend=True)

    saved["metallicity"] = ad["gas", "metallicity"]
    saved["masses"] = ad["gas", "masses"]

    return saved["metallicity"],saved["masses"]

    for i in sorted(saved):
        if not hasattr(saved[i], 'shape'): continue
        print "% 20s => %s" % (i, saved[i].shape)

