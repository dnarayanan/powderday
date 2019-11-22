from __future__ import print_function
import sys
import random
import numpy as np
import powderday.config as cfg
from powderday.gridstats import gridstats
from powderday.octree_zoom import octree_zoom_bbox_filter
from yt.geometry.selection_routines import AlwaysSelector
from powderday.dust_grid_gen import dtm_grid, remy_ruyer, manual,li_bestfit
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

    # load the DS and add pd fields; no need to put in stellar ages yet
    # as this will happen downstream in zoom
    pf = field_add(fname, bounding_box=bbox)

    # zoom if necessary
    # if cfg.par.zoom == True:
    pf = octree_zoom_bbox_filter(fname, pf, bbox, field_add)

    pf.index
    ad = pf.all_data()

    # ---------------------------------------------------------------
    # PLOTTING DIAGNOSTIC PROJECTION PLOTS
    # ---------------------------------------------------------------

    # proj_plots(pf)

    from yt.data_objects.particle_unions import ParticleUnion
    pu = ParticleUnion("all", list(pf.particle_types_raw))

    
    if  yt.__version__ == '4.0.dev0':
        #get central positions of cells as edges + width/2.
        xpos = pf.parameters['octree'][('index', 'x')]+(pf.parameters['octree'][('index','dx')]/2.)
        ypos = pf.parameters['octree'][('index', 'y')]+(pf.parameters['octree'][('index','dy')]/2.)
        zpos = pf.parameters['octree'][('index', 'z')]+(pf.parameters['octree'][('index','dz')]/2.)
        #comebine these into the fc1 array with dimensions (nparticles,3)
        fc1 = np.array([xpos,ypos,zpos]).T

        dx = pf.parameters['octree'][('index','dx')]
        dy = pf.parameters['octree'][('index','dy')]
        dz = pf.parameters['octree'][('index','dz')]
        fw1 = np.array([dx,dy,dz]).T

        refined = pf.parameters['octree'][('index','refined')].astype('bool')

        n_ref = pf.parameters['octree'].n_ref
        max_level = find_max_level(refined)#'max_level not yet implemented in yt4.x'
        #note, we could figure this out from the max number of trues in a row
        nocts = len(refined)-np.sum(refined)
        pdb.set_trace()
    else:
        saved = pf.index.oct_handler.save_octree()
        always = AlwaysSelector(None)
        #ir1 = pf.index.oct_handler.ires(always)  # refinement levels
        fc1 = pf.index.oct_handler.fcoords(always)  # coordinates in code_length
        fw1 = pf.index.oct_handler.fwidth(always)  # width of cell in code_length        
        refined = saved['octree'].astype('bool')
        
        n_ref = pf.index.oct_handler.n_ref
        max_level = pf.index.oct_handler.max_level
        nocts = pf.index.oct_handler.nocts
    # convert fc1 and fw1 to YTArrays
    fc1 = pf.arr(fc1, 'code_length')
    fw1 = pf.arr(fw1, 'code_length')



    print('----------------------------')
    print('yt Octree Construction Stats')
    print('----------------------------')
    print(' n_ref = ', n_ref)
    print(' max_level = ', max_level)
    print(' nocts = ', nocts)
    print('----------------------------')

    gridstats(fc1, fw1)

    # ==================================
    
    if  yt.__version__ == '4.0.dev0':
        refined = pf.parameters['octree'][('index','refined')].astype('bool')
    else:
        refined = saved['octree'].astype('bool')

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


    # smooth the data on to the octree

    volume = np.zeros(len(refined))

    if cfg.par.CONSTANT_DUST_GRID == False:

        # crash the code if the parameter choice for dust grid type isn't in 
        # the hard coded valid list below
        dust_grid_type_list = ['dtm', 'rr', 'manual','li_bestfit']
        try:
            dust_choice = dust_grid_type_list.index(cfg.par.dust_grid_type)
        except ValueError as e:
            print('You chose a dust_choice that isnt a valid selection within the list: dust_grid_type_list....crashing now!')
            sys.exit()

        if cfg.par.dust_grid_type == 'dtm':
            dust_smoothed_dtm = dtm_grid(pf, refined)
            dust_smoothed = dust_smoothed_dtm

        if cfg.par.dust_grid_type == 'rr':
            dust_smoothed_remy_ruyer = remy_ruyer(pf, refined)
            dust_smoothed = dust_smoothed_remy_ruyer

        if cfg.par.dust_grid_type == 'manual':
            dust_smoothed_manual = manual(pf, refined)
            dust_smoothed = dust_smoothed_manual
        
        if cfg.par.dust_grid_type == 'li_bestfit':
            dust_smoothed_li_bestfit = li_bestfit(pf,refined)
            dust_smoothed=dust_smoothed_li_bestfit



    else:
        print('cfg.par.CONSTANT_DUST_GRID=True')
        print('setting constant dust grid to 4.e-22')
        dust_smoothed = np.zeros(len(refined))+4.e-23

    # return refined,dust_smoothed,xmin,xmax,ymin,ymax,zmin,zmax,boost
    return refined, dust_smoothed, fc1, fw1, pf, ad


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
