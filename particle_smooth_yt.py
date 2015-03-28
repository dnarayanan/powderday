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

    #make sure the relevant fields make the transition
    pf.field_info["starmetals"].particle_type=True
    pf.field_info["gasmetals"].particle_type=True
    
    

    ad = pf.all_data()
    saved = pf.index.oct_handler.save_octree()
    
    saved["density"] = ad["gassmootheddensity"]
    saved["metallicity"] = ad["gassmoothedmetals"]
    saved["masses"] = ad["gassmoothedmasses"]
    
    
   
    
  

    #convert density to cgs
    saved["density"] = saved["density"].in_cgs()

    
 

    #return saved["metallicity"],saved["density"],saved["masses"],ad[metal_fn]
   
    return saved["metallicity"],saved["density"],saved["masses"]
   
    for i in sorted(saved):
        if not hasattr(saved[i], 'shape'): continue
        print "% 20s => %s" % (i, saved[i].shape)

