from common_analysis import *
import pylab as plt
import pdb,ipdb
import numpy as np
import cPickle
from progenner import *
from dn_galinfo_return import *
from progenner import *

def progen(SNAPNUM,p):

    ## specify snapshot number
    #SNAPNUM = 300

    ## what galaxy?
    #GAL = 14

    ## create a list of SPHGR objects using the
    ## loadSata function defined in common_analysis.py


    #get the galaxy ID from the last snap
    sims = cPickle.load(open(p,'rb'))
    galid,mstar,mgas,sfr,mhalo,haloid = galinfo_top(sims,10)
    GAL = galid[0]
    

    sims_loadsnap = iS.loadData(p)
    
    obj = sims_loadsnap

    #now get the progenitors for the galaxies and it's physical properties
    progenInit([sims_loadsnap], SNAPNUM, 0, [[GAL]], 'galaxy')
    
    #print whats available
    obj.galaxies[GAL].progen_galaxies[0].__dict__.keys()

    #assign the physical quantities

    z = np.asarray([s.redshift for s in obj.galaxies[GAL].progen_galaxies])
    
    

    cm = np.asarray([s.cm for s in obj.galaxies[GAL].progen_galaxies])
    cmx = cm[:,0]*0.7 #new progen has h divided out - need to put it back to put in code units
    cmy = cm[:,1]*0.7
    cmz = cm[:,2]*0.7
    
    snap = np.arange(len(z))
    snap += (SNAPNUM+1)-len(z)
    ipdb.set_trace()
    '''
    for i in range(0,len(gal.progen_z)):
        snap.append(SNAPNUM-i)
        print '%d %0.3f  %e %e  %f  %f  %f' % (snap[i], gal.progen_z[i],
                                               gal.progen_stellar_mass[i],
                                               gal.progen_halo_mass[i],
                                               gal.progen_cmx[i],
                                               gal.progen_cmy[i],
                                               gal.progen_cmz[i])
        
    '''
            
    return snap,cmx,cmy,cmz

