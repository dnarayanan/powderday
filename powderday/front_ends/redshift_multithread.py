from __future__ import print_function
from multiprocessing import Pool
import numpy as np
import powderday.config as cfg
from datetime import datetime
from astropy.cosmology import Planck13
from powderday.helpers import find_nearest

def redshift_multithread(formation_z):

    formation_z_list = formation_z
    #initialize the process pool and build the chunks
    p = Pool(processes = cfg.par.n_processes)
    nchunks = cfg.par.n_processes
    chunk_start_indices = []
    chunk_start_indices.append(0) #the start index is obviously 0
    delta_chunk_indices = int(len(formation_z_list) / nchunks)
    print ('delta_chunk_indices = ',delta_chunk_indices)
    for n in range(1,nchunks):
        chunk_start_indices.append(chunk_start_indices[n-1]+delta_chunk_indices)
    list_of_chunks = []
    for n in range(nchunks):
        formation_z_chunk = formation_z_list[chunk_start_indices[n]:chunk_start_indices[n]+delta_chunk_indices]
        if n == nchunks-1: 
            formation_z_chunk = formation_z_list[chunk_start_indices[n]::]
        list_of_chunks.append(formation_z_chunk)
    print ('Entering Pool.map multiprocessing for Stellar Age calculations')
    t1=datetime.now()
    chunk_sol = p.map(redshift_gen, [arg for arg in list_of_chunks])
    

    formation_time = []
    for i in range(len(chunk_sol)):
        sub_chunk_sol = chunk_sol[i].value
        for j in range(len(sub_chunk_sol)):
            formation_time.append(sub_chunk_sol[j])

    t2=datetime.now()
    print ('Execution time for Stellar Age calculations in Pool.map multiprocessing = '+str(t2-t1))
    
  
    
    return np.array(formation_time)
        

def redshift_gen(formation_z):
    age = Planck13.age(formation_z)
    
    return age


def redshift_vectorized(formation_z):
    t1=datetime.now()
    
    print ('Beginning redshift referencing')
    reference_z = np.arange(0,100,0.01)
    reference_t = Planck13.age(reference_z)
    t2=datetime.now()
    print ('Redshift referencing time: = '+str(t2-t1))
    
    formation_time = []

    t1=datetime.now()
    for i in range(len(formation_z)):
            
            idx = find_nearest(reference_z,formation_z[i])
            formation_time.append(reference_t[idx])
        

    t2=datetime.now()
    print ('redshift_vectorizing time: = '+str(t2-t1))

    return formation_time
