from __future__ import print_function
import yt
import powderday.config as cfg



def stream(fname):

   
    
    def gadget():
        if ('PartType0', 'CS Temperature') in ds.derived_field_list:
            from powderday.front_ends.CSgadget2pd import gadget_field_add as field_add
        elif ('PartType4', 'TemperatureMax') in ds.derived_field_list:
            from powderday.front_ends.benopp_gadget2pd import gadget_field_add as field_add
        else:
            from powderday.front_ends.gadget2pd import gadget_field_add as field_add
        
        print ('[front_end_controller:] gadget data set detected')
        return field_add

    def tipsy():
        from powderday.front_ends.tipsy2pd import tipsy_field_add as field_add
        print ('[front_end_controller:] tipsy data set detected')
        return field_add

    def ramses():
        from powderday.front_ends.ramses2pd import ramses_field_add as field_add
        print ('[front_end_controller:] ramses data set detected')
        return field_add

    
    def enzo():
        from powderday.front_ends.enzo2pd import enzo_field_add as field_add
        print ('[front_end_controller:] enzo data set detected')
        return field_add


    def arepo():
        from powderday.front_ends.arepo2pd import arepo_field_add as field_add
        print('[front_end_controller:] arepo data set detected')
        return field_add


    bbox = [[-2.*cfg.par.bbox_lim,2.*cfg.par.bbox_lim],
            [-2.*cfg.par.bbox_lim,2.*cfg.par.bbox_lim],
            [-2.*cfg.par.bbox_lim,2.*cfg.par.bbox_lim]]
    
    try: 
        ds = yt.load(fname,bounding_box = bbox)
        ds.index
        print ('[front_end_controller:] bounding_box being used')
    except:
        ds = yt.load(fname)
        ds.index
        print ('[front_end_controller:] NO bounding_box being used')

    ds_type = ds.dataset_type 
    
  
    #define the options dictionary
    options = {'gadget_hdf5':gadget,
               'tipsy':tipsy,
               'ramses':ramses,
               'enzo_packed_3d':enzo,
               'arepo_hdf5':arepo}


    #grab the field from the right front end
    field_add = options[ds_type]()
   
        
    return field_add,ds
