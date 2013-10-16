import numpy as np
cimport numpy as np
cimport cython
import parameters as par
#from libc.math cimport exp
from datetime import datetime
import pdb
from libc.stdlib cimport malloc, free

cdef extern from "math.h":
    double sqrt(double x) nogil
    double log10(double x) nogil
    double exp(double x) nogil


def particle_smooth_new(np.ndarray[double] x,
                        np.ndarray[double] y,
                        np.ndarray[double] z,
                        np.ndarray[double] hsml,
                        np.ndarray[double,ndim=2] coordinates,
                        np.ndarray[double] m,
                        refined,
                        np.ndarray[double] mass_grid):
 
    #define the grid in terms of refined_nomaster, to get the
    #dimensions right with respect to the coordinates array.  then,
    #later, we'll add a '0' back in the place of the base grid in
    #the final mass_grid

    
    cdef int p
   
    cdef double dist
    

    cdef int g
    cdef int counter=0

    cdef int i
    cdef int j

    cdef double a_norm
    cdef double hsml_var

    wTrue = np.where(np.array(refined) == True)[0]
    wFalse = np.where(np.array(refined) == False)[0]
    
    cdef int refined_len = len(refined)
    cdef int g_len = len(wFalse)
    cdef int nparticles = len(x)
    cdef double *kernel_sum = <double *>malloc(g_len *sizeof(double))
    cdef double total_kernel_sum = 0.


    #zero out the kernel_sum
    for i from 0<=i<g_len:
        kernel_sum[i] = 0
    

    #for each particle
    t1 = datetime.now()
    print 'particle_smooth_cython beginning at time: '+str(t1)


    if par.NPARTICLES_DEBUG != -1:
        print'================================================================='
        print 'WARNING WARNING WARNING WARNING WARNING'
        print'================================================================='
        print 'NPARTICLES_DEBUG IS SET TO A REAL NUMBER'
        print'================================================================='


        if par.NPARTICLES_DEBUG < nparticles: nparticles = par.NPARTICLES_DEBUG
        
    print 'particle_smooth_cython: using nparticles = '+str(nparticles)

   
    for p from 0<=p<nparticles:
 
        a_norm = m[p]
        hsml_var = hsml[p]
        
        for g from 0<=g<g_len:
            
            dist = ( (coordinates[g,0]-x[p])**2. +
                     (coordinates[g,1]-y[p])**2. +
                     (coordinates[g,2]-z[p])**2. )
            dist = dist**(0.5)/hsml_var
            
        
            
            #the line below is for gaussian smoothing - note, if we
            #were to use this, we'd have to not normalize the distance
            #by the hsml_var  as in the line above

            #kernel_sum[g] += a_norm * exp( (-1. * (dist**2.))/(2. * (hsml_var**2.)))

            kernel_sum[g] += sph_kernel(dist)

         
          


            

    


    t2 = datetime.now()
    print 'particle_smooth_cython: time for kernel smoothing  = '+str(t2-t1)

    print 'particle_smooth_cython: summing the kernel_sum'
   
    
    for i from 0<=i<g_len:
        total_kernel_sum += kernel_sum[i]
      


    print 'particle_smooth_cython: total_kernel_sum in code units = '+str(total_kernel_sum)

    print 'particle_smooth_cython: returning unnormalized kernel_sum to main'
    #copy the kernel sum over to the mass_grid
    for i from 0<=i<g_len:
        mass_grid[i] = kernel_sum[i]
        

     

    free(kernel_sum)
  

    
    return mass_grid
  
    

  
##############################################
#Standard SPH kernel for use with the Grid method
#Written by Bobby Thompson
cdef double sph_kernel(double x) nogil:
    cdef double kernel
    if x <= 0.5:
        kernel = 1.-6.*x*x*(1.-x)
    elif x>0.5 and x<=1.0:
        kernel = 2.*(1.-x)*(1.-x)*(1.-x)
    else:
        kernel = 0.
    return kernel
##############################################


     
