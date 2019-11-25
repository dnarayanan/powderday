from __future__ import print_function


def yt_smooth(pf):

    print ('starting yt_smooth')
    pf.index

    #make sure the relevant fields make the transition
    #pf.field_info["starmetals"].particle_type=True
    #pf.field_info["gasmetals"].particle_type=True
    
    
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
        print ("% 20s => %s" % (i, saved[i].shape))

