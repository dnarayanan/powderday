import numpy as np
from yt.mods import *
from yt.fields.particle_fields import add_volume_weighted_smoothed_field
import pdb
from octree_zoom import octree_zoom
import config as cfg
from gridstats import gridstats
from yt.geometry.selection_routines import AlwaysSelector
from yt.geometry.oct_container import OctreeContainer



def yt_smooth(pf):

    print 'starting yt_smooth'
    pf.index
    ad = pf.all_data()
    

   




    if cfg.par.zoom == True:


      
        pf = octree_zoom(pf)

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



        print '----------------------------'
        print 'Grid was Zoomed in with box length (kpc): ',cfg.par.zoom_box_len
        print 'A new Octree was built with the stats: '
        print ' n_ref = ',pf.index.oct_handler.n_ref
        print ' max_level = ',pf.index.oct_handler.max_level
        print ' nocts = ',pf.index.oct_handler.nocts
        print ' The new Grid Statistics Follow:'
        print '----------------------------'

        always = AlwaysSelector(None)
        ir1 = pf.index.oct_handler.ires(always)  #refinement levels
        fc1 = pf.index.oct_handler.fcoords(always)  #coordinates in kpc
        fw1 = pf.index.oct_handler.fwidth(always)  #width of cell in kpc
        gridstats(ir1,fc1,fw1)







    else:


        #Metallicity_00
        
        
        metal_fn = add_volume_weighted_smoothed_field("PartType0", "Coordinates", "Masses",
        "SmoothingLength", "Density","Metallicity_00",
        pf.field_info)
        '''
    

        metal_fn = add_volume_weighted_smoothed_field("PartType0", "Coordinates", "Masses",
                                                      "SmoothingLength", "Density","Metallicity",
                                                      pf.field_info)
   
        '''
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

