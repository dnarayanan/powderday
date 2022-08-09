import numpy as np
from powderday.pah.pah_file_read import read_draine_file
from powderday.helpers import find_nearest
from astropy import units as u
from astropy import constants as constants
import powderday.config as cfg
import pdb
from tqdm import tqdm
from powderday.pah.isrf_decompose import get_beta_nnls
import os,glob
from multiprocessing import Pool
from functools import partial
from datetime import datetime



#2a. have the code know if its neutral or ion and use that luminosity

#2b. manually add the logU>4 sources


def compute_grid_PAH_luminosity(cell_list,beta_nnls,grid_of_sizes,numgrains,draine_sizes,draine_lam,neutral_PAH_reference_objects,draine_bins_idx):

    #these are re-defined for each pool thread.  when they're
    #returned, they'll be packed into a master list with one extra
    #dimension for the thread number.
    grid_PAH_luminosity = np.zeros([len(cell_list),len(draine_lam)])


    for counter,cell in enumerate(cell_list):
        print(cell)
        beta_cell = beta_nnls[:,cell]
        beta_cell = beta_cell/np.max(beta_cell)
        
        #need to make a temporary (for this cell) PAH_list that is
        #just n_draine_sizes long that is convolved with beta_nnls
        pah_grid = np.zeros([len(draine_sizes),len(draine_lam)])
        for j in range(len(beta_cell)):
            PAH_list = neutral_PAH_reference_objects[j]
            temp_pah_grid = np.array([x.lum for x in PAH_list])
            temp_pah_grid *= beta_cell[j]
            pah_grid += temp_pah_grid #this is the running summation of the (n_sizes,n_lam) pah grid for the i-th cell
            
            
        #set the PAH luminosity of the cell to be the dot product of
        #the Draine luminosities (i.e., pah_grid[draine_bins_idx,:] which has
        #dimensions (simulation_sizes,wavelengths)) with the actual
        #grain size distribution in that cell (i.e.,
        #grid_of_sizes[i_cell,:]). note, we take the transpose of
        #grid_of_sizes to get the dimensions to match up correctly for the dot product
        

        grid_PAH_luminosity[counter,:] = np.dot(pah_grid[draine_bins_idx,:].T, grid_of_sizes.T[:,cell])
    #particle_PAH_luminosity = np.dot(pah_grid[draine_bins_idx,:].T,numgrains.T[:,cell])
    
    return grid_PAH_luminosity


def get_PAH_lum_cdf(nu_reverse,fnu,wpah_nu_reverse,grid_PAH_luminosity):
    lum_list = []

    for i in range(grid_PAH_luminosity.shape[0]): 
        fnu_reverse = fnu[i,:][::-1]
        lum = (np.absolute(np.trapz(nu_reverse[wpah_nu_reverse].cgs.value,fnu_reverse[wpah_nu_reverse])).item()*u.Lsun).to(u.erg/u.s).value
        lum_list.append(lum)

    lum_list = np.asarray(lum_list)
    lum_list[lum_list == 0] = np.min(lum_list[lum_list > 0])

    #omg this is so janky. but it works.
    percentile_list = []

    loglum_bins = np.linspace(np.min(np.log10(lum_list)),np.max(np.log10(lum_list)),100)
    for loglum in loglum_bins:
        percentile_list.append(np.sum(lum_list[lum_list > 10.**loglum])/np.sum(lum_list))
        
    percentile_idx = find_nearest(np.asarray(percentile_list),cfg.par.percentile_LPAH_to_include)
    lum_to_cut_below = 10.**(loglum_bins[percentile_idx])

    useful_idxs = np.where(lum_list >= lum_to_cut_below)[0]


    return useful_idxs


def pah_source_add(ds,reg,m,boost):
    
    LUM_FLOOR = 1.e20 #erg/s -- just some small value compared to the ~few Lsun we typically get in a cell

    
    #first - establish where we're working
    draine_directories = []
    print('powderday/pah/pah_source_create]: reading from the following Draine PAH directories')
    for it in os.scandir(cfg.par.draine_data_dir):
        if it.is_dir():
            print(it.path)
            draine_directories.append(it.path)


    #first establish the grain size distribution and sizes from the
    #hydro simulation
    grid_of_sizes = ds.parameters['reg_grid_of_sizes']

    simulation_sizes = (ds.parameters['grain_sizes_in_micron']*u.micron)


    #determine q_PAH for analysis and save it to parameters for
    #writing out 
    ad = ds.all_data()

    idx_pah = np.where(simulation_sizes.to(u.cm).value <= 3.e-7)[0]


    dN_pah = np.sum(reg['particle_dust','numgrains'][:,idx_pah],axis=1)
    dN_total = np.sum(reg['particle_dust','numgrains'],axis=1)

    q_pah = (dN_pah * reg['particle_dust','mass'])/(dN_total*reg['particle_dust','mass'])
    q_pah = q_pah * reg['particle_dust','carbon_fraction']

    reg.parameters['q_pah'] = q_pah
    

    #compute the mass weighted grain size distributions for comparison in analytics.py 
    #try: #for mesh based code
    #particle_mass_weighted_gsd = np.average(reg['dust','numgrains'],weights=reg['dust','mass'],axis=0)
    #grid_mass_weighted_gsd = np.average(grid_of_sizes,weights=reg['dust','mass'],axis=0)
    #except:
    particle_mass_weighted_gsd = np.average(reg['particle_dust','numgrains'],weights=reg['dust','mass'],axis=0)
    try: #for octree    
        grid_mass_weighted_gsd = np.average(grid_of_sizes,weights=reg['dust','smoothedmasses'],axis=0)
    except:
        grid_mass_weighted_gsd = np.average(grid_of_sizes,weights=reg['dust','mass'],axis=0)

    #second, read the information from the Draine files. We can do
    #this just for an arbitrary file in one of the draine_directories
    #since the size bins and wavelengths are all the same.
    temp_filename = glob.glob(draine_directories[0]+'/*iout_graD16*nb*_0.00')[0]
    temp_PAH_list = read_draine_file(temp_filename)
    draine_sizes = temp_PAH_list[0].size_list
    draine_lam = temp_PAH_list[0].lam*u.micron


    #now read in the full PAH_list for all logU = 0 files (note, we'll
    #fill in any logU>>4 PAH spectra on a case by case basis later in
    #this module. these are relatively rare, and not worth the effort
    #to carry around all that information otherwise. the logic here is
    #that the majority of spectra from logU=[0,4] are nearly
    #identical, and also the vast majority of cells are logU<4.  So we
    #can treat the logU>=4 on an indvidiaul basis.)
    
    neutral_logU_iout_files = []
    ion_logU_iout_files = []
    for directory in draine_directories:
        neutral_logU_iout_files.append(glob.glob(directory+'/*iout_graD16*nb*_0.00')[0])
        ion_logU_iout_files.append(glob.glob(directory+'/*iout_graD16*ib*_0.00')[0])

    neutral_PAH_reference_objects = np.zeros([len(neutral_logU_iout_files),len(temp_PAH_list)],dtype=object)
    ion_PAH_reference_objects = np.zeros([len(ion_logU_iout_files),len(temp_PAH_list)],dtype=object)

    print("[pah/pah_source_create:] building the reference PAH list for neutrals")
    for counter,neutral_file in tqdm(enumerate(neutral_logU_iout_files)):
        neutral_PAH_reference_objects[counter,:] = np.asarray(read_draine_file(neutral_file))

    print("[pah/pah_source_create:] building the reference PAH list for ions")
    for counter,ion_file in tqdm(enumerate(ion_logU_iout_files)):
        ion_PAH_reference_objects[counter,:] = np.asarray(read_draine_file(ion_file))

    #third, on a cell-by-cell basis, interpolate the luminosity for
    #each grain size bin, and multiply by the number of grains in that
    #bin
    ncells = grid_of_sizes.shape[0]

    total_PAH_luminosity = np.zeros(len(temp_PAH_list[0].lam))
    
    #get the logU and beta_nnls for the local ISRF
    beta_nnls,logU = get_beta_nnls(draine_directories,grid_of_sizes,simulation_sizes,reg)


    #find the indices of the Draine sizes that best match those that are in the simulation
    Draine_simulation_idx_left_edge_array = []
    for size in simulation_sizes.to(u.cm).value:
        idx0 = find_nearest(draine_sizes,size)
        #if draine_sizes[idx0] > size: idx0 -=1
        
        #this is really the nearest point in the Draine sizes to the
        #simulation_sizes. 
        Draine_simulation_idx_left_edge_array.append(idx0)



    #get the indices for where the Draine size bins match ours 
    size_arange = np.arange(len(simulation_sizes))
    draine_bins_idx = np.asarray(Draine_simulation_idx_left_edge_array)[size_arange]

    #pah_grid = np.array([x.lum for x in temp_PAH_list])


    #initialize the process pool and build the chunks
    t1 = datetime.now()
    nprocesses = np.min([cfg.par.n_processes,ncells]) #pool.map will barf in the corner case that we have less cells than cores


    p = Pool(processes = nprocesses)

    cell_list = np.arange(ncells)


    #chunking to speed up multiprocessing: since the processes are so
    #quick, we can lose a factor of 50% time in just spawning new
    #threads.  i saves a ton of time to chunk up the work and send it all off once.
    nchunks=nprocesses
    chunk_start_indices = []
    chunk_start_indices.append(0) #the start index is obviously 0
    #this should just be int(ncells/nchunks) but in case ncells < nchunks, we need to ensure that this is at least  1
    delta_chunk_indices = np.max([int(len(cell_list) / nchunks),1])
    print ('delta_chunk_indices = ',delta_chunk_indices)

    for n in range(1,nchunks):
        chunk_start_indices.append(chunk_start_indices[n-1]+delta_chunk_indices)

    list_of_chunks = []
    for n in range(nchunks):
        cells_list_chunk = cell_list[chunk_start_indices[n]:chunk_start_indices[n]+delta_chunk_indices]
        #if we're on the last chunk, we might not have the full list included, so need to make sure that we have that here
        if n == nchunks-1:
            cells_list_chunk = cell_list[chunk_start_indices[n]::]
        list_of_chunks.append(cells_list_chunk)



    print("Computing the PAH luminosities for every cell given its grain size distribution and logU. Entering Pool.map multiprocessing.")
    
    dum_numgrains = reg['particle_dust','numgrains'].value 


    
    #DEBUG DEBUG DEBUG UNCOMMENT THIS ENTIRE BLOCK TO GET BACK TO NORMAL PRODUCTION CODE WE JUST HAVE AN NPZ FILE FOR DEBUGGING RN
    '''
    temp_grid_PAH_luminosity = p.map(partial(compute_grid_PAH_luminosity,
                           beta_nnls = beta_nnls,
                           grid_of_sizes = grid_of_sizes.value,
                           numgrains = dum_numgrains,
                           draine_sizes = draine_sizes,
                           draine_lam = draine_lam.value,
                           neutral_PAH_reference_objects = neutral_PAH_reference_objects,
                           draine_bins_idx = draine_bins_idx),[arg for arg in list_of_chunks])
    
    #temp_grid_PAH_luminosity is now nprocesses long.  each index will
    #reveal a [ncells/nchunks,n_draine_lam] length array.  we hence
    #conctaenate them to get this to be ncells,ndraine_lam long.
    temp_grid_PAH_luminosity = np.asarray(temp_grid_PAH_luminosity)
    grid_PAH_luminosity = np.concatenate((temp_grid_PAH_luminosity),axis=0) 
    '''
    #DEBUG DEBUG DEBUG REMOVE THE NEXT TWO LINES THEY'RE CRAZY TALK AND JUST FOR DEBUGGING
    data = np.load('/blue/narayanan/desika.narayanan/pd_runs/powderday_testing/tests/SKIRT/MW_ultra_lowres_pah_sizes/grid_pah_luminosity.npz')
    grid_PAH_luminosity = data['grid_PAH_luminosity']


    t2 = datetime.now()
    print ('Execution time for PAH dot producting [is that a word?] across the grid = '+str(t2-t1))
 



   #z = zip(cell_list,beta_nnls,grid_of_sizes.value,dum_numgrains,draine_sizes,draine_lam.value,neutral_PAH_reference_objects,draine_bins_idx)
    #answer = p.starmap(compute_grid_PAH_luminosity,z)
    #p.close()
    #p.terminate()
    #p.join()

    #what i need back out is something like this with these dimensions
    #grid_PAH_luminosity = np.zeros([ncells,len(draine_lam)])
    #particle_PAH_luminosity = np.zeros([len(numgrains),len(draine_lam)])

    #grid_PAH_luminosity,particle_PAH_luminosity = compute_grid_PAH_luminosity(cell_list,beta_nnls,draine_sizes,draine_lam,neutral_PAH_reference_objects,grid_of_sizes,reg,draine_bins_idx)


    grid_PAH_luminosity[np.isnan(grid_PAH_luminosity)] = 0
    

    nu = (constants.c/draine_lam).to(u.Hz)
    #the units here are Lsun/Hz - this is to be consistent with our
    #stellar fnu addition later. The SEDs of individiual sources all
    #end up getting renormalized by the luminosity, so the exact units
    #don't matter as long as they're consistent across all the sources
    #(and types of sources) being added to the grid.
    fnu = np.divide((grid_PAH_luminosity*u.erg/u.s).to(u.Lsun).value,nu.to(u.Hz).value)

    #Because the Draine templates include re-emission, but we want to
    #add the PAHs as sources only, we restrict to the PAH range.
    nu_reverse = nu[::-1]

    #this is in here to set up the testing/debugging infrastructure
    #for a scenario where we only want to include certain wavelengths.
    #We set it to [0.1,1e3] as a default to catch all the emission.
    nu_pah2 = (constants.c/(0.1*u.micron)).to(u.Hz) #start of the pah range
    nu_pah1 = (constants.c/(1.e3*u.micron)).to(u.Hz) #end of pah range
    wpah_nu_reverse = np.where( (nu_reverse.value < nu_pah2.value) & (nu_reverse.value > nu_pah1.value))[0]




    #set a fnu floor since the 0's can propagate to NaNs
    fnu_floor = np.min(fnu[fnu>0])
    fnu[fnu==0]=fnu_floor

    #to reduce memory requirements, we can't really add a PAH source
    #for every single cell.  so we just do so for cells that are at a
    #luminosity such that the CDF (sum(L>L_threshold) > 99% of the
    #total luminosity).  
 
    only_important_PAH_idx = get_PAH_lum_cdf(nu_reverse,fnu,wpah_nu_reverse,grid_PAH_luminosity)
    
    for i in only_important_PAH_idx:#range(grid_PAH_luminosity.shape[0]): #np.arange(2500)
        #lum = np.trapz(draine_lam[wpah_lam].cgs.value,flam[i,wpah_lam]).value
        fnu_reverse = fnu[i,:][::-1]
        #if np.where(fnu_reverse == 0)[0] > 0:
        #    fnu_reverse[fnu_reverse ==0 ] = np.min(fnu_reverse[fnu_reverse > 0])

        lum = (np.absolute(np.trapz(nu_reverse[wpah_nu_reverse].cgs.value,fnu_reverse[wpah_nu_reverse])).item()*u.Lsun).to(u.erg/u.s).value


        if lum <= LUM_FLOOR: lum = LUM_FLOOR #just a jamky variable
                                            #defined at the top of
                                            #this function to define a
                                            #lowest luminosity so that
                                            #we don't add PAH cells
                                            #with 0 luminosity
        #reversing arrays to make nu increasing, and therefore correct for hyperion addition
        m.add_point_source(luminosity=lum,spectrum=(nu_reverse[wpah_nu_reverse].value,fnu_reverse[wpah_nu_reverse]),position=reg.parameters['cell_position'][i,:].in_units('cm').value-boost)


        #reg.parameters['cell_position'][i,:].in_units('cm').value-boost)
        #reg['particle_dust','coordinates'][i,:].in_units('cm').value-boost)


    if cfg.par.draine21_pah_grid_write: #else, the try/except in analytics.py will get caught and will just write a single -1 to the output npz file
        reg.parameters['grid_PAH_luminosity'] = grid_PAH_luminosity

    reg.parameters['PAH_lam'] = draine_lam.value

    total_PAH_luminosity =np.sum(grid_PAH_luminosity,axis=0)
    reg.parameters['total_PAH_luminosity'] = total_PAH_luminosity
    reg.parameters['only_important_PAH_idx'] = only_important_PAH_idx

    grid_PAH_L_lam = grid_PAH_luminosity/draine_lam.value
    integrated_grid_PAH_luminosity = np.trapz((grid_PAH_luminosity/draine_lam.value),draine_lam.value,axis=1)
    reg.parameters['integrated_grid_PAH_luminosity'] = integrated_grid_PAH_luminosity
    
    #save some information for dumping into analytics
    reg.parameters['q_pah'] = q_pah
    reg.parameters['particle_mass_weighted_gsd'] = particle_mass_weighted_gsd
    reg.parameters['grid_mass_weighted_gsd'] = grid_mass_weighted_gsd
    reg.parameters['simulation_sizes'] = simulation_sizes


#    for i_cell in tqdm(range(ncells)):
        
        
        #in principle this doesn't need to be done inside the loop for
        #a constant radiation field. However, for a radiation field
        #that varies cell by cell, then PAH_list will change as we
        #have different files that we read in, so we may as well keep it here for now.

 #       pah_grid = np.array([x.lum for x in PAH_list])
 #       idx = np.asarray(Draine_simulation_idx_left_edge_array)[size_arange]
        
        #set the PAH luminosity of the cell to be the dot product of
        #the Draine luminosities (i.e., pah_grid[idx,:] which has
        #dimensions (simulation_sizes,wavelengths)) with the actual
        #grain size distribution in that cell (i.e.,
        #grid_of_sizes[i_cell,:]). note, we take the transpose of
        #grid_of_sizes to get the dimensions to match up correctly for the dot product

  
  #      grid_PAH_luminosity[i_cell,:] = np.dot(pah_grid[idx,:].T,grid_of_sizes[i_cell,:].T)

  


    #import matplotlib.pyplot as plt
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.loglog(PAH_list[0].lam,total_PAH_luminosity[:]/PAH_list[0].lam)
    #ax.set_ylim([1e31,1e45])
    #ax.set_xlim([1,1000])
    #ax.set_xlabel(r'$\lambda (\mu $m)')
    #ax.set_ylabel(r'$L_\lambda$ (erg/s/$\mu$m)')
    #fig.savefig('/home/desika.narayanan/PAH_sed.png',dpi=300)
    
