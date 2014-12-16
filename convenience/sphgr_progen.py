from common_analysis import *
import pylab as plt
import pdb,ipdb
import numpy as np

def progen(SNAPNUM,GAL):

    ## specify snapshot number
    #SNAPNUM = 300

    ## what galaxy?
    #GAL = 14

    ## create a list of SPHGR objects using the
    ## loadSata function defined in common_analysis.py
    sims = loadSnap(SNAPNUM)
    
    ## depending on your DIRPRES config, this list
    ## may have multiple objects, but for now lets
    ## focus on the 0th index
    obj = sims[0]

    ## import progen routine
    from progenner import *
    ## initialize progen data
    ## this basically fills the galaxy object with
    ## all of the .progen_xxxx data
    progenInit(sims, SNAPNUM, 0, [[GAL]], 'galaxy')

    snap = []
    gal = obj.galaxies[GAL]
    for i in range(0,len(gal.progen_z)):
        snap.append(SNAPNUM-i)
        print '%d %0.3f  %e %e  %f  %f  %f' % (snap[i], gal.progen_z[i],
                                               gal.progen_stellar_mass[i],
                                               gal.progen_halo_mass[i],
                                               gal.progen_cmx[i],
                                               gal.progen_cmy[i],
                                               gal.progen_cmz[i])
        

    return snap,gal.progen_cmx,gal.progen_cmy,gal.progen_cmz
