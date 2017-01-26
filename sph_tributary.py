import sys
script,pardir,parfile,modelfile = sys.argv
import numpy as np
import scipy.interpolate
import scipy.ndimage
import os.path
import copy
import pdb,ipdb

from hyperion.model import Model
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from hyperion.model import ModelOutput
import config as cfg
import yt
from yt.units.yt_array import YTQuantity
from astropy.table import Table
from astropy.io import ascii

import h5py
from grid_construction import yt_octree_generate,grid_coordinate_boost,grid_center
import SED_gen as sg
from find_order import *
import powderday_test_octree as pto
import hyperion_octree_stats as hos
import error_handling as eh
import backwards_compatibility as bc



def sph_m_gen(fname,field_add):
    
    refined,dustdens,fc1,fw1,pf,ad = yt_octree_generate(fname,field_add)
    xmin = (fc1[:,0]-fw1[:,0]/2.).convert_to_units('cm') #in proper cm 
    xmax = (fc1[:,0]+fw1[:,0]/2.).convert_to_units('cm')
    ymin = (fc1[:,1]-fw1[:,1]/2.).convert_to_units('cm')
    ymax = (fc1[:,1]+fw1[:,1]/2.).convert_to_units('cm')
    zmin = (fc1[:,2]-fw1[:,2]/2.).convert_to_units('cm')
    zmax = (fc1[:,2]+fw1[:,2]/2.).convert_to_units('cm')
    

    #dx,dy,dz are the edges of the parent grid
    dx = (np.max(xmax)-np.min(xmin)).value
    dy = (np.max(ymax)-np.min(ymin)).value
    dz = (np.max(zmax)-np.min(zmin)).value


    xcent = np.mean([np.min(xmin),np.max(xmax)]) #kpc
    ycent = np.mean([np.min(ymin),np.max(ymax)])
    zcent = np.mean([np.min(zmin),np.max(zmax)])
    
    boost = np.array([xcent,ycent,zcent])
    print '[pd_front end] boost = ',boost

    
    #Tom Robitaille's conversion from z-first ordering (yt's default) to
    #x-first ordering (the script should work both ways)

    refined_array = np.array(refined)
    refined_array = np.squeeze(refined_array)
    
    order = find_order(refined_array)
    refined_reordered = []
    dustdens_reordered = np.zeros(len(order))
    
    
    
    for i in range(len(order)): 
        refined_reordered.append(refined[order[i]])
        dustdens_reordered[i] = dustdens[order[i]]


    refined = refined_reordered
    dustdens=dustdens_reordered

    #hyperion octree stats
    max_level = hos.hyperion_octree_stats(refined)


    pto.test_octree(refined,max_level)

    
    np.save('refined.npy',refined)
    np.save('density.npy',dustdens)
    

    #========================================================================
    #Initialize Hyperion Model
    #========================================================================

    m = Model()
    
    if cfg.par.FORCE_RANDOM_SEED == True: m.set_seed(cfg.par.seed)

    print 'Setting Octree Grid with Parameters: '



    #m.set_octree_grid(xcent,ycent,zcent,
    #                  dx,dy,dz,refined)
    m.set_octree_grid(0,0,0,dx/2,dy/2,dz/2,refined)    


    if cfg.par.PAH == True:
        frac = {'usg': 0.0586, 'vsg': 0.1351, 'big': 0.8063}
        for size in ['usg', 'vsg', 'big']:
            m.add_density_grid(dustdens * frac[size], cfg.par.dustdir+'%s.hdf5' % size)
        
        m.set_enforce_energy_range(cfg.par.enforce_energy_range)
    else:
        m.add_density_grid(dustdens,cfg.par.dustdir+cfg.par.dustfile)



  





    return m,xcent,ycent,zcent,dx,dy,dz,pf,boost
