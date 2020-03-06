from __future__ import print_function
import numpy as np



def gridstats(fc1,fw1):

    volume = (fw1**3)
    xmin = fc1[:,0]-fw1[:,0]/2.
    xmax = fc1[:,0]+fw1[:,0]/2.
    ymin = fc1[:,1]-fw1[:,1]/2.
    ymax = fc1[:,1]+fw1[:,1]/2.
    zmin = fc1[:,2]-fw1[:,2]/2.
    zmax = fc1[:,2]+fw1[:,2]/2.


    print ('----------------------------')
    print ('Grid Statistics')
    print ('----------------------------')
    
    print ('Smallest Cell Edge: ',np.min(fw1.in_units('kpc')))
    print ('Smallest Cell Volume: ',np.min(volume.in_units('kpc**3')))

    print ('Biggest Cell Edge: ',np.max(fw1.in_units('kpc')))
    print ('Biggest Cell Volume: ',np.max(volume.in_units('kpc**3')))
   
    print ('Left Edge: ',np.min(xmin.in_units('kpc')))
    print ('Right Edge: ',np.max(xmax.in_units('kpc')))
    print ('Bottom Edge: ',np.min(ymin.in_units('kpc')))
    print ('Top Edge: ',np.max(ymax.in_units('kpc')))
    print ('Nearest Edge: ',np.min(zmin.in_units('kpc')))
    print ('Farthest Edge: ',np.max(ymax.in_units('kpc')))
    
