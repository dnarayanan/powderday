from __future__ import print_function
import numpy as np
import powderday.config as cfg
import yt
from yt.frontends.sph.data_structures import ParticleDataset
from yt.geometry.selection_routines import AlwaysSelector
import pdb


ParticleDataset.filter_bbox = True
ParticleDataset._skip_cache = True


def octree_zoom_bbox_filter(fname,ds,bbox0,field_add):

    
    ds.index
    ad = ds.all_data()

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
    #yt 3.x
    box_len = ds.quan(box_len,'kpc')
    #yt 4.x
    if yt.__version__ == '4.0.dev0':
        box_len = float(box_len.to('code_length').value)
        bbox_lim = box_len
    else:
        box_len = box_len.convert_to_units('code_length').value
        bbox_lim = box_len
    
    bbox1 = [[center[0]-bbox_lim,center[0]+bbox_lim],
            [center[1]-bbox_lim,center[1]+bbox_lim],
            [center[2]-bbox_lim,center[2]+bbox_lim]]
    print ('[octree zoom] new zoomed bbox (comoving/h) in code units= ',bbox1)
    

    #yt 3.x
    #ds1 = yt.load(fname,bounding_box=bbox1,n_ref = cfg.par.n_ref,over_refine_factor=cfg.par.oref)
    


    #What follows is tricky.  Broadly, the plan is to create a yt
    #region to cut out the dataset to our desired box size.  In yt4.x,
    #we will then pass around reg (which represents the cutout version
    #of the ds), as well as ds (which is the original ds).  the
    #original ds will contain all the original parameters of the
    #dataset.  We pass around the octree itself in a newly created
    #dictionary called reg.parameters

    if yt.__version__ == '4.0.dev0':
        
        #re load the field names, but now with the bounding box
        #set. this will allow us to map the field names to those
        #generated in the octree.  this represents a massive
        #inefficiency as we have to load the entire dataset a *second*
        #time.
        ds = field_add(fname,bounding_box = bbox1,ds=ds,add_smoothed_quantities=True)
        ds.periodicity = (False,False,False)
        reg = ds.region(center=center,left_edge = np.asarray(center)-bbox_lim,right_edge = np.asarray(center)+bbox_lim)

        #ds1 = reg.ds
        left = np.array([pos[0] for pos in bbox1])
        right = np.array([pos[1] for pos in bbox1])
        octree = ds.octree(left, right, n_ref=cfg.par.n_ref)#, force_build=True)

        reg.parameters={}
        reg.parameters['octree'] = octree

 

    else:
        #load up a cutout ds1 with a bounding_box so we can generate the octree on this dataset
        ds1 = yt.load(fname,bounding_box = bbox1,n_ref=cfg.par.n_ref,over_refine_factor=cfg.par.oref) 
        ds1.periodicity = (False,False,False)
        #now update the field names
        ds1 = field_add(None,bounding_box = bbox1,ds=ds1,add_smoothed_quantities=True)

        #now create the region so that we have the smoothed properties downstream correct
        reg = ds1.region(center=center,left_edge = np.asarray(center)-bbox_lim,right_edge = np.asarray(center)+bbox_lim)
        reg.parameters={}

        saved = ds1.index.oct_handler.save_octree()
        always = AlwaysSelector(None)
        #ir1 = ds.index.oct_handler.ires(always)  # refinement levels
        reg.parameters["fc1"] = ds1.index.oct_handler.fcoords(always)  # coordinates in code_length
        reg.parameters["fw1"] = ds1.index.oct_handler.fwidth(always)  # width of cell in code_length
        reg.parameters["refined"] = saved['octree'].astype('bool')
        
        reg.parameters["n_ref"] = ds1.index.oct_handler.n_ref
        reg.parameters["max_level"] = ds1.index.oct_handler.max_level
        reg.parameters["nocts"] = ds1.index.oct_handler.nocts



    #re-add the new powderday convention fields; this time we need to
    #make sure to do the ages calculation since it hasn't been done
    #before.
    #ds1 = field_add(None,bounding_box = bbox1,ds=ds1,starages=True)

    
    return reg



def arepo_zoom(fname,ds,bbox0,field_add):


    ds.index
    ad = ds.all_data()

    center = [cfg.model.x_cent,cfg.model.y_cent,cfg.model.z_cent]
    print ('[zoom/arepo_zoom:] using center: ',center)


    box_len = cfg.par.zoom_box_len

    #now begin the process of converting box_len to physical units in
    #case we're in a cosmological simulation.  We'll first give it
    #units of proper kpc, then convert to code length (which for
    #gadget is kpcm/h) for the bbox calculation (dropping the units of
    #course).  then when we re-convert to proper units, the box_len as
    #input in parameters_master will be in proper units.  if a
    #simulation isn't cosmological, then the only difference here will
    #be a 1/h

    #note: we dispense with the yt3.x options that are written in
    #octree_zoom_bbox_filter since you can't read in an arepo sim as
    #an arepo model without yt4.x (if it's yt3.x, it'll be read in via the sph_tributary

    box_len = ds.quan(box_len,'kpc')
    box_len = float(box_len.to('code_length').value)
    bbox_lim = box_len
   
    bbox1 = [[center[0]-bbox_lim,center[0]+bbox_lim],
            [center[1]-bbox_lim,center[1]+bbox_lim],
            [center[2]-bbox_lim,center[2]+bbox_lim]]
    print ('[zoom/arepo_zoom] new zoomed bbox (comoving/h) in code units= ',bbox1)

    #re load the field names, but now with the bounding box
    #set. this will allow us to map the field names to those
    #generated in the octree.  this represents a massive
    #inefficiency as we have to load the entire dataset a *second*
    #time.
    ds = field_add(fname,bounding_box = bbox1,ds=ds)
    ds.periodicity = (False,False,False)
    reg = ds.region(center=center,left_edge = np.asarray(center)-bbox_lim,right_edge = np.asarray(center)+bbox_lim)


    return reg
