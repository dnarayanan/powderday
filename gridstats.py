import numpy as np
import config as cfg
import pdb
import constants as const

def gridstats(ir1,fc1,fw1):

    volume = (fw1 * const.pc* 1.e3)**3.
    
    xmin = fc1[:,0]-fw1[:,0]/2.
    xmax = fc1[:,0]+fw1[:,0]/2.
    ymin = fc1[:,1]-fw1[:,1]/2.
    ymax = fc1[:,1]+fw1[:,1]/2.
    zmin = fc1[:,2]-fw1[:,2]/2.
    zmax = fc1[:,2]+fw1[:,2]/2.


    print '----------------------------'
    print 'Grid Statistics'
    print '----------------------------'
    
    print 'Smallest Cell Edge (kpc): ',np.min(fw1)
    print 'Smallest Cell Volume (cm^3): ',np.min(volume)

    print 'Biggest Cell Edge (kpc): ',np.max(fw1)
    print 'Biggest Cell Volume (cm^3): ',np.max(volume)
   
    print 'Left Edge (kpc): ',np.min(xmin)
    print 'Right Edge (kpc): ',np.max(xmax)
    print 'Bottom Edge (kpc): ',np.min(ymin)
    print 'Top Edge (kpc): ',np.max(ymax)
    print 'Nearest Edge (kpc): ',np.min(zmin)
    print 'Farthest Edge (kpc): ',np.max(ymax)
    
