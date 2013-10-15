import random
import numpy as np
import pfh_readsnap
import parameters as par
from datetime import datetime
from astropy.table import Table
from astropy.io import ascii
import constants as const
import pdb
import math
import sys

def chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

def particle_smooth_linalg(x,y,z,hsml,coordinates,pos,m,refined):
 
    #define the grid in terms of refined_nomaster, to get the
    #dimensions right with respect to the coordinates array.  then,
    #later, we'll add a '0' back in the place of the base grid in
    #the final mass_grid.  
    
    nparticles = len(hsml)

    refined_nomaster = refined[1::]
    
    mass_grid = np.zeros(len(refined_nomaster))
    mass_gauss = np.zeros(len(refined_nomaster))

    wTrue = np.where(np.array(refined_nomaster) == True)[0]
    wFalse = np.where(np.array(refined_nomaster) == False)[0]
    

    #create chunks in sizes nparticles/nchunk
    master_chunk_array = chunks(range(long(nparticles)),long(par.NCHUNK))

   
    print 'in smooth_operator.particle_smooth_linalg: dividing the particles into '+str(nparticles/par.NCHUNK)+' chunks'

    chunk_counter=0
    t1 = datetime.now()
    for chunk_list in master_chunk_array:
        print 'chunk_counter = '+str(chunk_counter)

      
        print 'assigning distance arrays'
        dumpos=pos[chunk_list]
    

    
        dx = (dumpos[:,0,np.newaxis]-coordinates[np.newaxis,:,0])
        dy = (dumpos[:,1,np.newaxis]-coordinates[np.newaxis,:,1])
        dz = (dumpos[:,2,np.newaxis]-coordinates[np.newaxis,:,2])
        
       
  
        distance = np.sqrt(dx**2.+dy**2.+dz**2.)
 
        
       

        #p is the particle number
        for p in chunk_list:
            
        
          
            
            
            a_norm = m[p]
            #initialize all the values to be 0 just in case in mass_gauss - a temporary array holding the particle smoothed values
            mass_gauss[:] = 0
            
            #index distance by p-chunk_list[0] since distance is only
            #n_elements(chunk_list) big in its first dimension, so need to have it be in terms of 0's
            mass_gauss[wFalse] = a_norm * np.exp( (-1. * (distance[p-chunk_list[0],:]**2.))/(2. * (hsml[p]**2.)))
            
            #normalize to conserve mass in the smoothing
           
        
            mass_grid += mass_gauss
                
           
            
      
        chunk_counter += 1


            
        t2 = datetime.now()
        print 'in smooth_operator.particle_smooth_linalg: time taken for kernel smoothing = '+str(t2-t1)

    
    pdb.set_trace()

    #normalize to conserve mass in the smoothing
    mass_grid /= sum(mass_grid)/sum(m)





  

    assert( (sum(mass_grid)/sum(m[0:p+1]) < 1.01) & (sum(mass_grid)/sum(m[0:p+1]) > 0.99))
    assert( (sum(mass_gauss)/m[p] < 1.01) & (sum(mass_gauss)/m[p] > 0.99))

    #append a 0 to the beginning to account for the base grid
    mass_grid = np.insert(mass_grid,0,0,0)


        
 
    return mass_grid
    
    



