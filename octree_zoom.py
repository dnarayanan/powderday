
import numpy as np
import pdb
from yt.mods import *
from yt.fields.particle_fields import add_volume_weighted_smoothed_field
import sys
import config as cfg


def octree_zoom(pf,box_len):


    pf.index


    
    region = pf.region([0,0,0],[-box_len,-box_len,-box_len],[box_len,box_len,box_len])
    
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
    bbox = [[-box_len*1.01, box_len*1.01],
            [-box_len*1.01, box_len*1.01],
            [-box_len*1.01, box_len*1.01]]




    new_ds = load_particles(data,
                            length_unit = pf.length_unit,
                            mass_unit = pf.mass_unit,
                            velocity_unit = pf.velocity_unit,
                            bbox=np.array(bbox),
                            over_refine_factor=cfg.par.oref)

    new_ds.particle_types_raw = tuple(pt for pt in pf.particle_types_raw 
                                      if pt not in skip)
    new_ds.index
                 


    



    return new_ds
