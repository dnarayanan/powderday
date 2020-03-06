import yt

#function takes a ds region as the input
def yt_field_map(obj):
    ds = obj.ds
    # First we figure out all the fields we want to cut out.
    field_list = ds.field_info.field_aliases.items()
    
    data = {}
    done = set([])
    skip = set([])
  
    pb = yt.get_pbar("Converting fields ", len(field_list))
    for i, ((field_type, field_name), (disk_type, disk_name)) in \
            enumerate(sorted(field_list)):
        pb.update(i)
        # We are going to skip any derived fields that slipped in.
        if (disk_type, disk_name) not in ds.field_list: continue
        if field_type in ds.particle_types and \
           field_type not in ds.particle_types_raw:
            continue
        d = obj[field_type, field_name]
        units = ds.field_info[disk_type, disk_name].units
        if d.size == 0:
            skip.add(field_type)
        if len(d.shape) == 2 and d.shape[1] == 3:
            # We split this up into multiple pieces as well as single.
            for i, ax in enumerate('xyz'):
                data[field_type, "%s_%s" % (field_name, ax)] = \
                    obj[field_type, field_name][:,i].in_units(units)
        data[field_type, field_name] = \
            obj[field_type, field_name].in_units(units)
        done.add((disk_type, disk_name))
    pb.finish()
    return data, skip
