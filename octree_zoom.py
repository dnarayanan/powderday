
import numpy as np
import pdb
from yt.mods import *
from yt.fields.particle_fields import add_volume_weighted_smoothed_field
import sys
import config as cfg



def octree_zoom(pf):


    pf.index
    ad = pf.all_data()

    print '----------------------------'
    print 'Entering Octree Zoom with parameters: '
    print "(...Calculating Center of Mass in octree_zoom)"
    com = ad.quantities.center_of_mass()
    print "Center of Mass is at coordinates (kpc): ",com
   



    minbox = np.array(com)-cfg.par.zoom_box_len
    maxbox = np.array(com)+cfg.par.zoom_box_len
    region = pf.region(com,minbox,maxbox)

    print 'minimum edges of the zoomed box are: (kpc)',minbox
    print 'maximum edges of the zoomed box are: (kpc)',maxbox

    data = {}
    skip = []
    vector_fields = [(r"particle_position_%s", "code_length"),
                #(r"particle_velocity_%s", "code_length/code_time")
                     ]
                     

    for ptype in pf.particle_types_raw:
        for fn, u in vector_fields:
            for ax in 'xyz':
                if region[ptype, fn % ax].size == 0:
                     skip.append(ptype)
                     continue
                data[ptype, fn % ax] = region[ptype, fn % ax].in_units(u)
                data[ptype, fn % ax] = (data[ptype, fn % ax].d, u)
                print ptype, fn % ax, data[ptype, fn % ax][0].min(), data[ptype, fn % ax][0].max()
            fn = fn[:-3]
            if ptype in skip: continue
            region.field_data.clear()
            data[ptype, fn] = region[ptype, fn].in_units(u)
            data[ptype, fn] = (data[ptype, fn].d, u)
            print ptype, fn, data[ptype, fn][0].min(), data[ptype, fn][0].max()

    for field, disk in sorted(pf.field_info.field_aliases.items()):
        if field in data: continue
        if field[0] not in pf.particle_types_raw: continue
        if field[1] == "smoothing_length":
            data[field] = region[field].in_units("code_length")
        else:
            data[field] = region[field].in_cgs()
        print "Setting", field, data[field].units
        data[field] = (data[field].d, str(data[field].units))


    import pprint
    pprint.pprint(sorted(data.keys()))

    #because minbox can be negative or positive (as can maxbox), and
    #we want to make sure we go *just* beyond those values for bbox to
    #encapsulate all of the particles in region, we have to have some
    #np.min and max arguments in the bbox definition.

    bbox = [[np.min([minbox[0]*1.01,minbox[0]*0.99]), 
             np.max([maxbox[0]*1.01,maxbox[0]*0.99])],
            [np.min([minbox[1]*1.01,minbox[1]*0.99]),
             np.max([maxbox[1]*1.01,maxbox[1]*0.99])],
            [np.min([minbox[2]*1.01,minbox[2]*0.99]),
             np.max([maxbox[2]*1.01,maxbox[2]*0.99])]]

 


    new_ds = load_particles(data,
                            length_unit = pf.length_unit,
                            mass_unit = pf.mass_unit,
                            velocity_unit = pf.velocity_unit,
                            bbox=np.array(bbox),
                            over_refine_factor=cfg.par.oref)

    new_ds.particle_types_raw = tuple(pt for pt in pf.particle_types_raw 
                                      if pt not in skip)
    new_ds.index
                


    '''
    DEBUG
    '''

    new_ds.field_info.alias(("gas", "density"),
    ('deposit', 'PartType0_smoothed_density'))
    saved = new_ds.index.oct_handler.save_octree(always_descend=True)
    saved["density"] = ad["gas","density"]
 
    new_ad = new_ds.all_data()


    print '[min,max] x particle position: ',np.min(new_ad["PartType0","particle_position_x"]),np.max(new_ad["PartType0","particle_position_x"])
    print '[min,max] y particle position: ',np.min(new_ad["PartType0","particle_position_y"]),np.max(new_ad["PartType0","particle_position_y"])
    print '[min,max] z particle position: ',np.min(new_ad["PartType0","particle_position_z"]),np.max(new_ad["PartType0","particle_position_z"])

  
    #p = ProjectionPlot(new_ds, "z", "density")
    #p.save()



    



    return new_ds
