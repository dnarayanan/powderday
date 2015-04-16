import numpy as np
import yt
import ipdb

def stream(fname):

   
    
    def gadget():
        if ('PartType0', 'CS Temperature') in ds.derived_field_list:
            from CSgadget2pd import gadget_field_add as field_add
        else:
            from gadget2pd import gadget_field_add as field_add
        
        print '[front_end_controller:] gadget data set detected'
        return field_add

    def tipsy():
        from tipsy2pd import tipsy_field_add as field_add
        print '[front_end_controller:] tipsy data set detected'
        return field_add

    def ramses():
        from ramses2pd import ramses_field_add as field_add
        print '[front_end_controller:] ramses data set detected'
        return field_add

  
    ds = yt.load(fname)
    ds_type = ds.dataset_type 
    
  
    #define the options dictionary
    options = {'gadget_hdf5':gadget,
               'tipsy':tipsy,
               'ramses':ramses}


    #grab the field from the right front end
    field_add = options[ds_type]()
    
        
    return field_add
