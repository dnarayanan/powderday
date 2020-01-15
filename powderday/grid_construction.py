from __future__ import print_function
import sys
import random
import numpy as np
import powderday.config as cfg
from powderday.gridstats import gridstats
from powderday.zoom import octree_zoom_bbox_filter,arepo_zoom
from yt.geometry.selection_routines import AlwaysSelector
#octree quantities for dust
from powderday.dust_grid_gen import dtm_grid_oct, remy_ruyer_oct, manual_oct,li_bestfit_oct,li_ml_oct

#particle and/or mesh quantities for dust
from powderday.dust_grid_gen import dtm_particle_mesh,remy_ruyer_particle_mesh,li_bestfit_particle_mesh,li_ml_particle_mesh

import yt
import pdb

random.seed('octree-demo')


def find_max_level(refined):
    #loop through refined and figure out how many trues we can have in a row
    master_max_level = 0
    max_level = 0.
    for i in range(len(refined)-1):
        if refined[i+1] == True:
            max_level+=1
        else:
            max_level = 0
        if max_level > master_max_level: master_max_level = max_level
        if max_level > 20: pdb.set_trace()
    return master_max_level


def yt_octree_generate(fname, field_add):

    print('[grid_construction]: bbox_lim = ', cfg.par.bbox_lim)

    bbox = [[-2.*cfg.par.bbox_lim, 2.*cfg.par.bbox_lim],
            [-2.*cfg.par.bbox_lim, 2.*cfg.par.bbox_lim],
            [-2.*cfg.par.bbox_lim, 2.*cfg.par.bbox_lim]]

    # load the DS and add pd fields.  this is the first field addition
    # of the simulation, so we don't yet need to add the smoothed
    # quantities (which can take some time in yt4.x).
    ds = field_add(fname, bounding_box=bbox,add_smoothed_quantities=False)

        

    #now zoom in.
    reg = octree_zoom_bbox_filter(fname, ds, bbox, field_add)

    
    # ---------------------------------------------------------------
    # PLOTTING DIAGNOSTIC PROJECTION PLOTS
    # ---------------------------------------------------------------
    # proj_plots(ds)


    #from yt.data_objects.particle_unions import ParticleUnion
    #pu = ParticleUnion("all", list(ds.particle_types_raw))

    
    if  yt.__version__ == '4.0.dev0':
        refined = reg.parameters['octree'][('index','refined')].astype('bool')

        #get central positions of cells as edges + width/2.
        #        xpos = (ds.parameters['octree'][('index', 'x')]+(ds.parameters['octree'][('index','dx')]/2.))[~refined]
        #        ypos = (ds.parameters['octree'][('index', 'y')]+(ds.parameters['octree'][('index','dy')]/2.))[~refined]
        #        zpos = (ds.parameters['octree'][('index', 'z')]+(ds.parameters['octree'][('index','dz')]/2.))[~refined]

        xpos = (reg.parameters['octree'][('index','x')])[~refined]
        ypos = (reg.parameters['octree'][('index','y')])[~refined]
        zpos = (reg.parameters['octree'][('index','z')])[~refined]

        #comebine these into the fc1 array with dimensions (nparticles,3)
        fc1 = np.array([xpos,ypos,zpos]).T

        dx = reg.parameters['octree'][('index','dx')][~refined]
        dy = reg.parameters['octree'][('index','dy')][~refined]
        dz = reg.parameters['octree'][('index','dz')][~refined]
        fw1 = np.array([dx,dy,dz]).T


        n_ref = reg.parameters['octree'].n_ref
        #max_level = find_max_level(refined)#'max_level not yet implemented in yt4.x'
        #note, we could figure this out from the max number of trues in a row
        nocts = len(refined)-np.sum(refined)


        
    else:
        #saved = ds.index.oct_handler.save_octree()
        #always = AlwaysSelector(None)
       
        #ir1 = ds.index.oct_handler.ires(always)  # refinement levels
        fc1 = reg.parameters["fc1"] #coordinates in code_length
        fw1 = reg.parameters["fw1"] #width of cell in code_length
        refined = reg.parameters["refined"]
        
        n_ref = reg.parameters["n_ref"]
        max_level = reg.parameters["max_level"]
        nocts = reg.parameters["nocts"]
    # convert fc1 and fw1 to YTArrays
    fc1 = ds.arr(fc1, 'code_length')
    fw1 = ds.arr(fw1, 'code_length')



    print('----------------------------')
    print('yt Octree Construction Stats')
    print('----------------------------')
    print(' n_ref = ', n_ref)
    #print(' max_level = ', max_level)
    print(' nocts = ', nocts)
    print('----------------------------')

    gridstats(fc1, fw1)

    # ==================================
    
    #try:
    #    refined.reshape(len(ir1))
    #except:
    refinements = 2**(3*cfg.par.oref)
    refined2 = []
    for r in refined:
        if r == 1:
            refined2.append(True)
        if r == 0:
            refined2.append(np.zeros(refinements).astype('bool'))
    refined = np.hstack(refined2)


    if cfg.par.CONSTANT_DUST_GRID == False:

        # crash the code if the parameter choice for dust grid type isn't in 
        # the hard coded valid list below
        dust_grid_type_list = ['dtm', 'rr', 'manual','li_bestfit','li_ml']
        try:
            dust_choice = dust_grid_type_list.index(cfg.par.dust_grid_type)
        except ValueError as e:
            print('You chose a dust_choice that isnt a valid selection within the list: dust_grid_type_list....crashing now!')
            sys.exit()

        if cfg.par.dust_grid_type == 'dtm':
            dust_smoothed_dtm = dtm_grid_oct(reg, refined)
            dust_smoothed = dust_smoothed_dtm

        if cfg.par.dust_grid_type == 'rr':
            dust_smoothed_remy_ruyer = remy_ruyer_oct(reg, refined)
            dust_smoothed = dust_smoothed_remy_ruyer

        if cfg.par.dust_grid_type == 'manual':
            dust_smoothed_manual = manual_oct(reg, refined)
            dust_smoothed = dust_smoothed_manual
        
        if cfg.par.dust_grid_type == 'li_bestfit':
            dust_smoothed_li_bestfit = li_bestfit_oct(reg,refined)
            dust_smoothed=dust_smoothed_li_bestfit

        if cfg.par.dust_grid_type == 'li_ml':
            dust_smoothed_li_ml = li_ml_oct(reg,refined)
            dust_smoothed = dust_smoothed_li_ml

    else:
        print('cfg.par.CONSTANT_DUST_GRID=True')
        print('setting constant dust grid to 4.e-22')
        dust_smoothed = np.zeros(len(refined))+4.e-22


    # return refined,dust_smoothed,xmin,xmax,ymin,ymax,zmin,zmax,boost
    return refined, dust_smoothed, fc1, fw1, reg, ds



def arepo_vornoi_grid_generate(fname, field_add):

    print('[grid_construction/arepo_vornoi_grid_generate]: bbox_lim = ', cfg.par.bbox_lim)

    bbox = [[-2.*cfg.par.bbox_lim, 2.*cfg.par.bbox_lim],
            [-2.*cfg.par.bbox_lim, 2.*cfg.par.bbox_lim],
            [-2.*cfg.par.bbox_lim, 2.*cfg.par.bbox_lim]]

    # load the DS and add pd fields.  this is the first field addition
    # of the simulation, so we don't yet need to add the smoothed
    # quantities (which can take some time in yt4.x).
    ds = field_add(fname, bounding_box=bbox)

    #now zoom in.
    reg = arepo_zoom(fname, ds, bbox, field_add)


    # ---------------------------------------------------------------
    # PLOTTING DIAGNOSTIC PROJECTION PLOTS
    # ---------------------------------------------------------------
    # proj_plots(ds)




    
    if cfg.par.CONSTANT_DUST_GRID == False:

        # crash the code if the parameter choice for dust grid type isn't in
        # the hard coded valid list below
        dust_grid_type_list = ['dtm', 'rr', 'manual','li_bestfit','li_ml']
        try:
            dust_choice = dust_grid_type_list.index(cfg.par.dust_grid_type)
        except ValueError as e:
            print('You chose a dust_choice that isnt a valid selection within the list: dust_grid_type_list....crashing now!')
            sys.exit()

        if cfg.par.dust_grid_type == 'dtm':
            dustdens = dtm_particle_mesh(reg)


        if cfg.par.dust_grid_type == 'rr':
            dustdens = remy_ruyer_particle_mesh(reg)


        if cfg.par.dust_grid_type == 'manual':
            raise ValueError(' "manual" dust grids not currently supported with Arepo simulations. Please try another choice amongst [dtm, rr, li_bestfit, li_ml]')
        #if cfg.par.dust_grid_type == 'manual':
        #    dust_smoothed_manual = manual(reg, refined)
        #    dust_smoothed = dust_smoothed_manual


        if cfg.par.dust_grid_type == 'li_bestfit':
            dustdens = li_bestfit_particle_mesh(reg)

        if cfg.par.dust_grid_type == 'li_ml':
            dustdens = li_ml_particle_mesh(reg)


    else:
        print('cfg.par.CONSTANT_DUST_GRID=True')
        print('setting constant dust grid to 4.e-22')
        dustdens = np.zeros(len(reg["gasmasses"]))+4.e-22



    #DEBUG -- this will eventually go into a refactored dust_grid_gen
#    metaldens = reg["metaldens"]
#    dustdens = (metaldens*cfg.par.dusttometals_ratio).to('g/cm**3').value

    return reg,ds,dustdens







def grid_coordinate_boost(xmin, xmax, ymin, ymax, zmin, zmax):

    print('\n boosting coordinates to [0,0,0] centering \n')
    xcent, ycent, zcent, dx, dy, dz = grid_center(
        xmin, xmax, ymin, ymax, zmin, zmax)
    xmin -= xcent
    xmax -= xcent
    ymin -= ycent
    ymax -= ycent
    zmin -= zcent
    zmax -= zcent

    return xmin, xmax, ymin, ymax, zmin, zmax


def stars_coordinate_boost(star_list, boost):

    # center the stars
    nstars = len(star_list)
    for i in range(nstars):
        star_list[i].positions[0] -= boost[0]
        star_list[i].positions[1] -= boost[1]
        star_list[i].positions[2] -= boost[2]
    return star_list


def grid_center(xmin, xmax, ymin, ymax, zmin, zmax):

    xcent = np.mean([min(xmin), max(xmax)])
    ycent = np.mean([min(ymin), max(ymax)])
    zcent = np.mean([min(zmin), max(zmax)])

    # dx,dy,dz are the edges of the parent grid
    dx = (max(xmax)-min(xmin))/2
    dy = (max(ymax)-min(ymin))/2.
    dz = (max(zmax)-min(zmin))/2.

    return xcent, ycent, zcent, dx, dy, dz
