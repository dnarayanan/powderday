import numpy as np
from yt.mods import *
from yt.fields.particle_fields import add_volume_weighted_smoothed_field
import pdb
from octree_zoom import octree_zoom
import config as cfg


def yt_smooth(pf):

    print 'starting yt_smooth'
    pf.index
    ad = pf.all_data()
    
   
   




    if cfg.par.zoom == True:
        pf = octree_zoom(pf,cfg.par.zoom_box_len)

        #these have slightly different field names (yt internal names)
        #than if we don't call the octree_zoom function
        #[i.e. metallicity instead of Metallicity; particle_mass
        #instead of Masses).
        metal_fn = add_volume_weighted_smoothed_field("PartType0", "Coordinates", "Masses",
                                                      "SmoothingLength", "Density","metallicity",
                                                      pf.field_info)
        mass_fn = add_volume_weighted_smoothed_field("PartType0", "Coordinates", "Masses",
                                                     "SmoothingLength", "Density","particle_mass", 
                                                     pf.field_info)
        density_fn = add_volume_weighted_smoothed_field("PartType0", "Coordinates", "Masses",
                                                        "SmoothingLength", "Density","density", 
                                                        pf.field_info)

        pf.field_info.alias(("gas", "metallicity"), metal_fn[0])
        pf.field_info.alias(("gas", "particle_mass"), mass_fn[0])
        pf.field_info.alias(("gas","density"),density_fn[0])

        saved = pf.index.oct_handler.save_octree(always_descend=True)
        saved["metallicity"] = ad["gas", "metallicity"]
        saved["masses"] = ad["gas", "particle_mass"]
        saved["density"] = ad["density"]


    else:


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



        density_fn = add_volume_weighted_smoothed_field("PartType0", "Coordinates", "Masses",
                                                        "SmoothingLength", "Density","Density", 
                                                        pf.field_info)
    
        pf.field_info.alias(("gas", "metallicity"), metal_fn[0])
        pf.field_info.alias(("gas", "masses"), mass_fn[0])
        pf.field_info.alias(("gas","density"),density_fn[0])
        
        saved = pf.index.oct_handler.save_octree(always_descend=True)
        saved["metallicity"] = ad["gas", "metallicity"]
        saved["masses"] = ad["gas", "masses"]
        saved["density"] = ad["density"]





   
    return saved["metallicity"],saved["masses"],saved["density"]

    for i in sorted(saved):
        if not hasattr(saved[i], 'shape'): continue
        print "% 20s => %s" % (i, saved[i].shape)

