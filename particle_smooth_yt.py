import numpy as np
from yt.mods import *
from yt.fields.particle_fields import add_volume_weighted_smoothed_field
import pdb
import config as cfg
from gridstats import gridstats
from yt.geometry.selection_routines import AlwaysSelector
from yt.geometry.oct_container import OctreeContainer


def yt_smooth(pf):

    print 'starting yt_smooth'
    pf.index
    ad = pf.all_data()
    

   




    if cfg.par.zoom == True:



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

        pf.field_info.alias(("gas", "density"),('deposit', 'PartType0_smoothed_density'))
        saved = pf.index.oct_handler.save_octree(always_descend=True)
        saved["density"] = ad["gas","density"]



        print '\n\n\n'
        print "[particle_smooth_yt] ---------------------------------------------"
        print "[particle_smooth_yt] MANUALLY SETTING ALL METALLICITIES TO SOLAR"
        print "[particle_smooth_yt] THIS IS NOT LONG TERM CORRECT"
        print "[particle smooth yt] Also note the masses are not being passed around"
        print "[particle_smooth_yt] ---------------------------------------------"
        print '\n\n\n'


        saved["metallicity"] = np.ones(len(ad["density"]))*0.02
       

        ''' UNCOMMENT THIS LINE OUT WHEN THIS KLUGE IS REVERTED
        saved["metallicity"] = ad["gas", "metallicity"]
        '''

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





   
    return saved["metallicity"],saved["density"]

    for i in sorted(saved):
        if not hasattr(saved[i], 'shape'): continue
        print "% 20s => %s" % (i, saved[i].shape)

