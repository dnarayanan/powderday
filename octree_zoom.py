import numpy as np
import pdb
from yt.mods import *
from yt.fields.particle_fields import add_volume_weighted_smoothed_field
import sys

bbox = [[-3220.0, 3220.0],
        [-3220.0, 3220.0],
        [-3220.0, 3220.0]]

unit_base = {'UnitLength_in_cm'         : 3.08568e+21,
             'UnitMass_in_g'            :   1.989e+43,
             'UnitVelocity_in_cm_per_s' :      100000}

pf = load('/Users/desika/gadgetruns/sbw_tests/mw_18_6_hr_hightimeres/snapshot_017.hdf5',
          unit_base = unit_base,
          bounding_box = bbox,
          over_refine_factor = 0)



pf.index


region = pf.region([0,0,0],[-10,-10,-10],[10,10,10])

data = {}
#data = dict((field, region[field].in_cgs().d) for field in pf.field_list if
#             field[0] != "all")
skip = []
vector_fields = [
                 (r"particle_position_%s", "code_length"),
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
bbox = [[-10.01, 10.01],
        [-10.01, 10.01],
        [-10.01, 10.01]]


new_ds = load_particles(data,
                        length_unit = unit_base['UnitLength_in_cm'],
                        mass_unit = unit_base['UnitMass_in_g'],
                        velocity_unit = unit_base['UnitVelocity_in_cm_per_s'],
                        bbox=np.array(bbox),
                        over_refine_factor=0)
new_ds.particle_types_raw = tuple(pt for pt in pf.particle_types_raw 
                                  if pt not in skip)
new_ds.index

ad = new_ds.all_data()

new_ds.field_info.alias(("gas", "metallicity"),
    ('deposit', 'PartType0_smoothed_metallicity'))
new_ds.field_info.alias(("gas", "density"),
    ('deposit', 'PartType0_smoothed_density'))


saved = new_ds.index.oct_handler.save_octree(always_descend=True)

saved["density"] = ad["gas","density"]

#saved["metallicity"] = ad["gas", "metallicity"]

p = ProjectionPlot(new_ds, "z", "density")
p.save()
