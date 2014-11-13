import numpy as np
from yt.mods import *
from yt.fields.particle_fields import add_volume_weighted_smoothed_field
import pdb
import config as cfg
from gridstats import gridstats
from yt.geometry.selection_routines import AlwaysSelector
from yt.geometry.oct_container import OctreeContainer
import yt.units as units
import ipdb

def yt_smooth(pf):

    print 'starting yt_smooth'
    pf.index
    
    if  ('PartType4', 'Metallicity_00') in pf.derived_field_list:
        pf.field_info["PartType0","Metallicity_00"].particle_type=True
    else: 
        pf.field_info["PartType0","Metallicity"].particle_type=True

    

    ad = pf.all_data()
    saved = pf.index.oct_handler.save_octree()


    saved["density"] = ad["deposit","PartType0_smoothed_density"]
    saved["metallicity"] = ad["deposit","PartType0_smoothed_metallicity"]
    saved["masses"] = ad["deposit", "PartType0_smoothed_particle_mass"]

    
  

    #convert density to cgs
    saved["density"] = saved["density"].in_cgs()

    
    #direct calculation of the smoothed metal density (via the added
    #MetalDens field in grid_construction)
    metal_fn = add_volume_weighted_smoothed_field("PartType0", "Coordinates", "Masses",
                                                  "SmoothingLength", "Density","MetalDens",
                                                  pf.field_info)
    metal_fn = metal_fn[0]



    return saved["metallicity"],saved["density"],saved["masses"],ad[metal_fn]
   
    #return saved["metallicity"],saved["density"],saved["masses"]
   
    for i in sorted(saved):
        if not hasattr(saved[i], 'shape'): continue
        print "% 20s => %s" % (i, saved[i].shape)

