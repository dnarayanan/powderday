import numpy as np
from yt.mods import *
from yt.fields.particle_fields import add_volume_weighted_smoothed_field
import pdb
import config as cfg
from gridstats import gridstats
from yt.geometry.selection_routines import AlwaysSelector
from yt.geometry.oct_container import OctreeContainer
import yt.units as units

def yt_smooth(pf):

    print 'starting yt_smooth'
    pf.index
    

    pf.field_info["PartType0","metallicity"].particle_type=True
    ad = pf.all_data()
    saved = pf.index.oct_handler.save_octree()


    saved["density"] = ad["deposit","PartType0_smoothed_density"]
    saved["metallicity"] = ad["deposit","PartType0_smoothed_metallicity"]
    saved["masses"] = ad["deposit", "PartType0_smoothed_particle_mass"]

   


    #convert density to cgs
    saved["density"] = saved["density"].in_cgs()

    return saved["metallicity"],saved["density"],saved["masses"]

   
    for i in sorted(saved):
        if not hasattr(saved[i], 'shape'): continue
        print "% 20s => %s" % (i, saved[i].shape)

