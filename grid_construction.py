from __future__ import print_function
import random
import numpy as np
#import parameters as par
import config as cfg
from datetime import datetime
from astropy.table import Table
from astropy.io import ascii
import hyperion_octree_stats as hos
from gridstats import gridstats
from octree_zoom import octree_zoom_bbox_filter
from analytics import proj_plots,dust_histograms
from yt.mods import *
from yt.geometry.oct_container import OctreeContainer
from yt.geometry.selection_routines import AlwaysSelector
from yt.fields.particle_fields import add_volume_weighted_smoothed_field
from dust_grid_gen import dtm_grid,remy_ruyer,manual

import astropy.constants as constants



random.seed('octree-demo')
 
import pdb,ipdb
import os.path
import sys



def yt_octree_generate(fname,field_add):
    
   
    

    

    print ('[grid_construction]: bbox_lim = ',cfg.par.bbox_lim)

    
    bbox = [[-2.*cfg.par.bbox_lim,2.*cfg.par.bbox_lim],
            [-2.*cfg.par.bbox_lim,2.*cfg.par.bbox_lim],
            [-2.*cfg.par.bbox_lim,2.*cfg.par.bbox_lim]]
 
             
    
    '''
    unit_base = {'UnitLength_in_cm'         : cfg.par.unit_length*1.e3*const.pc,
                 'UnitMass_in_g'            : cfg.par.unit_mass*const.msun,
                 'UnitVelocity_in_cm_per_s' : cfg.par.unit_velocity}
    '''

    #load the DS and add pd fields; no need to put in stellar ages yet
    #as this will happen downstream in zoom
    pf = field_add(fname,bounding_box = bbox)
    

    
    #zoom if necessary
    # if cfg.par.zoom == True:
    pf = octree_zoom_bbox_filter(fname,pf,bbox,field_add)

    
    
 

    pf.index
    ad = pf.all_data()

    
   
    #---------------------------------------------------------------
    #PLOTTING DIAGNOSTIC PROJECTION PLOTS
    #---------------------------------------------------------------
    
    #proj_plots(pf)






    from yt.data_objects.particle_unions import ParticleUnion
    pu = ParticleUnion("all", list(pf.particle_types_raw))
    
    saved = pf.index.oct_handler.save_octree()
 
    always = AlwaysSelector(None)
    ir1 = pf.index.oct_handler.ires(always)  #refinement levels
    fc1 = pf.index.oct_handler.fcoords(always)  #coordinates in code_length
    fw1 = pf.index.oct_handler.fwidth(always)  #width of cell in code_length
    
    #convert fc1 and fw1 to YTArrays
    fc1 = pf.arr(fc1,'code_length')
    fw1 = pf.arr(fw1,'code_length')




    print ('----------------------------')
    print ('yt Octree Construction Stats')
    print ('----------------------------')
    print (' n_ref = ',pf.index.oct_handler.n_ref)
    print (' max_level = ',pf.index.oct_handler.max_level)
    print (' nocts = ',pf.index.oct_handler.nocts)
    print ('----------------------------')
        
    gridstats(ir1,fc1,fw1)
    
    



    #==================================



    refined = saved['octree']
    refined2 = []

 
    for i in range(len(refined)):
        if refined[i] == 1: refined2.append(True)
        if refined[i] == 0: refined2.append(False)

    refined = refined2



    #smooth the data on to the octree
    
    volume = np.zeros(len(refined))


    if cfg.par.CONSTANT_DUST_GRID == False: 
        
        #dust_histograms(refined,dust_smoothed_dtm,dust_smoothed_remy_ruyer)

        #crash the code if the parameter choice for dust grid type isn't in the hard coded valid list below
        dust_grid_type_list = ['dtm','rr','manual']
        try:
            dust_choice = dust_grid_type_list.index(cfg.par.dust_grid_type)
        except ValueError as e:
            print ('You chose a dust_choice that isnt a valid selection within the list: '+dust_grid_type_list+'....crashing now!')
            sys.exit()


        if cfg.par.dust_grid_type == 'dtm': 
            dust_smoothed_dtm = dtm_grid(pf,refined)
            dust_smoothed = dust_smoothed_dtm

        if cfg.par.dust_grid_type == 'rr': 
            dust_smoothed_remy_ruyer = remy_ruyer(pf,refined)
            dust_smoothed = dust_smoothed_remy_ruyer

        if cfg.par.dust_grid_type == 'manual': 
            dust_smoothed_manual = manual(pf,refined)
            dust_smoothed = dust_smoothed_manual


    else:
        print ('cfg.par.CONSTANT_DUST_GRID=True')
        print ('setting constant dust grid to 4.e-22')
        dust_smoothed = np.zeros(len(refined))+4.e-23




   

  

   


    #return refined,dust_smoothed,xmin,xmax,ymin,ymax,zmin,zmax,boost
    return refined,dust_smoothed,fc1,fw1,pf,ad


def grid_coordinate_boost(xmin,xmax,ymin,ymax,zmin,zmax):
    
    print ('\n boosting coordinates to [0,0,0] centering \n')
    xcent,ycent,zcent,dx,dy,dz = grid_center(xmin,xmax,ymin,ymax,zmin,zmax)
    xmin -= xcent
    xmax -= xcent
    ymin -= ycent
    ymax -= ycent
    zmin -= zcent
    zmax -= zcent
                
    return xmin,xmax,ymin,ymax,zmin,zmax


def stars_coordinate_boost(star_list,boost):
    
#center the stars
    nstars = len(star_list)
    for i in range(nstars):
        star_list[i].positions[0] -= boost[0]
        star_list[i].positions[1] -= boost[1]
        star_list[i].positions[2] -= boost[2]
    return star_list



def grid_center(xmin,xmax,ymin,ymax,zmin,zmax):
    
    xcent = np.mean([min(xmin),max(xmax)])
    ycent = np.mean([min(ymin),max(ymax)])
    zcent = np.mean([min(zmin),max(zmax)])
    
    
    #dx,dy,dz are the edges of the parent grid
    dx = (max(xmax)-min(xmin))/2
    dy = (max(ymax)-min(ymin))/2.
    dz = (max(zmax)-min(zmin))/2.

    return xcent,ycent,zcent,dx,dy,dz




