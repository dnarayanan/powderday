from __future__ import print_function
import numpy as np
import powderday.config as cfg
import yt
from yt.frontends.sph.data_structures import ParticleDataset



ParticleDataset.filter_bbox = True
ParticleDataset._skip_cache = True


def octree_zoom_bbox_filter(fname,pf,bbox0,field_add):

    ds0 = pf
    
    ds0.index
    ad = ds0.all_data()

    print ('\n\n')
    print ('----------------------------')
    print ("[octree zoom_bbox_filter:] Calculating Center of Mass")


    gas_com_x = np.sum(ad["gasdensity"] * ad["gascoordinates"][:,0])/np.sum(ad["gasdensity"])
    gas_com_y = np.sum(ad["gasdensity"] * ad["gascoordinates"][:,1])/np.sum(ad["gasdensity"])
    gas_com_z = np.sum(ad["gasdensity"] * ad["gascoordinates"][:,2])/np.sum(ad["gasdensity"])


    com = [gas_com_x,gas_com_y,gas_com_z]

    print ("[octree zoom_bbox_filter:] Center of Mass is at coordinates (kpc): ",com)


    center = [cfg.model.x_cent,cfg.model.y_cent,cfg.model.z_cent]
    print ('[octree zoom_bbox_filter:] using center: ',center)


    box_len = cfg.par.zoom_box_len
    #now begin the process of converting box_len to physical units in
    #case we're in a cosmological simulation.  We'll first give it
    #units of proper kpc, then convert to code length (which for
    #gadget is kpcm/h) for the bbox calculation (dropping the units of
    #course).  then when we re-convert to proper units, the box_len as
    #input in parameters_master will be in proper units.  if a
    #simulation isn't cosmological, then the only difference here will
    #be a 1/h
    
    box_len = ds0.quan(box_len,'kpc')
    box_len = box_len.convert_to_units('code_length').value
    bbox_lim = box_len

    
    bbox1 = [[center[0]-bbox_lim,center[0]+bbox_lim],
            [center[1]-bbox_lim,center[1]+bbox_lim],
            [center[2]-bbox_lim,center[2]+bbox_lim]]
    print ('[octree zoom] new zoomed bbox (comoving/h) in code units= ',bbox1)
    

    try: #particle
        ds1 = yt.load(fname,bounding_box=bbox1,n_ref = cfg.par.n_ref,over_refine_factor=cfg.par.oref)
    except: #amr
        ds1 = yt.load(fname,n_ref = cfg.par.n_ref,over_refine_factor=cfg.par.oref)
        bbox1 = None

    ds1.periodicity = (False,False,False)

    #re-add the new powderday convention fields; this time we need to
    #make sure to do the ages calculation since it hasn't been done
    #before.
    ds1 = field_add(None,bounding_box = bbox1,ds=ds1,starages=True)

    
    return ds1
