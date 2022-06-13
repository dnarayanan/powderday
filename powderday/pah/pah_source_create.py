import numpy as np
from powderday.pah.pah_file_read import read_draine_file
from powderday.helpers import find_nearest
from astropy import units as u
from astropy import constants as constants
import powderday.config as cfg
import pdb
from tqdm import tqdm
from powderday.pah.isrf_decompose import get_beta_nnls
import os
#to do:

#0 DONE---- treat the PAH emission spectra as sources and add them as such to pd so that they attenuate accordingly

#1. speed up the convolution with grain size distribution -- currently
#takes ~2 hours per galaxy. Assume the slow part is the np.dot product.

#to do this:

#a. DONE --first, we need to have all the PAH filenames and SEDs on disk.

#b. DONE --then we need to figure out how to make those SEDs a basis
#function, and how to break up one of our ISRFs into a combination of
#those basis functions. -- in get_beta_nnls

#c. PAH_list needs to be read in for every PAH_file.  ideally we would
#have some way of ensuring that the PAH_files are in the same order as
#the ISRF directories -- perhaps the thing to do here is to read in
#the directories, and then pass that list into isrf_decompose to
#ensure that they're the same.

#c. isrf_decompose/get_beta_nnls takes in (once modified) the ISRF
#grid for the entire simulation, the GSD for every cell (once
#modified), and the dust mass, and returns the beta contribution of
#all the basis vectors, as well as the logU.  so now we need to shape
#that into place to return that properly here, and then for every
#cell, get the betas and log Us to combine the PAH emission files
#right. q

#4. put in the ionization fraction and appropriate file name per cell.


filename = '/blue/narayanan/desika.narayanan/powderday_files/PAHs/dataverse_files/mMMP/iout_graD16emtPAHib_mmpisrf_1.00'
#filename = '/blue/narayanan/desika.narayanan/powderday_files/PAHs/BC03_Z0.02_10Myr/iout_graD16emtPAHib_bc03_z0.02_1e7_1.00'
#filename = '/blue/narayanan/desika.narayanan/powderday_files/PAHs/BC03_Z0.0004_10Myr/iout_graD16emtPAHib_bc03_z0.0004_1e7_1.50'
#iout_DH20Ad_P0.20_0.00_bc03_z0.02_1e7_1.00'
#filename = '/home/desika.narayanan/PAHs/vsg_stat_therm.iout'

def pah_source_add(ds,reg,m,boost):

    
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
    #writing out DEBUG - WE SHOULD CHANGE THIS TO INCLUDE A FEW
    #DIFFRENT POSSIBILITIES, IUNCLUDING (A) COMPUTING QPAH AS IS, AND
    #(B) COMPUTING QPAH DIRECTLY FROM THE SIMULATION IN THE POSSIBLE
    #CASE THAT IT EXPLICITLY MODELS AROMATIC GRAPHITES

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

    #second, read the information from the Draine files
    PAH_list = read_draine_file(filename)
    print("reading Draine File",filename)
    draine_sizes = PAH_list[0].size_list
    draine_lam = PAH_list[0].lam*u.micron
    #third, on a cell-by-cell basis, interpolate the luminosity for
    #each grain size bin, and multiply by the number of grains in that
    #bin


    ncells = grid_of_sizes.shape[0]

    #define the grid that will store the PAH spectra for the entire mesh
    grid_PAH_luminosity = np.zeros([ncells,len(PAH_list[0].lam)])
    total_PAH_luminosity = np.zeros(len(PAH_list[0].lam))



    #get the logU and beta_nnls for the local ISRF
    beta_nnls,logU = get_beta_nnls(draine_directories,grid_of_sizes,simulation_sizes,reg)
    pdb.set_trace()


    #find the indices of the Draine sizes that best match those that are in the simulation
    Draine_simulation_idx_left_edge_array = []
    for size in simulation_sizes.to(u.cm).value:
        idx0 = find_nearest(draine_sizes,size)
        #if draine_sizes[idx0] > size: idx0 -=1
        
        #this is really the nearest point in the Draine sizes to the
        #simulation_sizes.  if we uncomment the above line, it will
        #become the left point.
        Draine_simulation_idx_left_edge_array.append(idx0)



    #this pah_grid is a n_draine_sizes length list of PAH SEDs
    #(i.e. the PAH SED for every draine size for a single draine
    #file).  

    #REMOVE THIS COMMENT WHEN THIS IS DONE: THIS NEEDS TO BE CHANGED
    #IN TWO WAYS: A. CURRENTLY THIS IS (167,2500) WHICH IS
    #N_DRAINE_SIZES,N_DRAINE_WAVELENGTHS.  This has to happen inside a
    #for(ncells) for loop. The reason for this is that the altnerative
    #is to dump it all inside an
    #(ncells,n_draine_sizes,n_draine_wavelengths) array, but this
    #would take even in small runs something like 50 GB memory.
    #better to lose some speed here.  a for loop for a small grid will
    #take 30 min, so if we can throw this into a pool.map that'll be
    #the business.

    #B. The key here is that this can't just work for a single Draine
    #file-- it has to have the beta_nnls from all the draine files
    #(where the correct draine file is known by the logU of the cell).
    #in essence, all the magic happens in that convolution.

    #C. how to actually do the dot product over all cells is not yet clear ot me.  it might be worth doing some dummy array testing.  

    size_arange = np.arange(len(simulation_sizes))
    pah_grid = np.array([x.lum for x in PAH_list])
    idx = np.asarray(Draine_simulation_idx_left_edge_array)[size_arange]

    pdb.set_trace()

    #set the PAH luminosity of the cell to be the dot product of
    #the Draine luminosities (i.e., pah_grid[idx,:] which has
    #dimensions (simulation_sizes,wavelengths)) with the actual
    #grain size distribution in that cell (i.e.,
    #grid_of_sizes[i_cell,:]). note, we take the transpose of
    #grid_of_sizes to get the dimensions to match up correctly for the dot product

    
    grid_PAH_luminosity = np.dot(pah_grid[idx,:].T, grid_of_sizes.T).T
    particle_PAH_luminosity = np.dot(pah_grid[idx,:].T,reg['particle_dust','numgrains'].T).T

    flam = (np.divide(particle_PAH_luminosity,draine_lam.cgs.value))
    fnu = draine_lam.cgs.value**2/constants.c.cgs.value*flam
    nu = (constants.c/draine_lam).to(u.Hz)

    #Because the Draine templates include re-emission, but we want to
    #add the PAHs as sources only, we restrict to the PAH range.
    nu_reverse = nu[::-1]
    nu_pah2 = (constants.c/(3*u.micron)).to(u.Hz) #start of the pah range
    nu_pah1 = (constants.c/(20.*u.micron)).to(u.Hz) #end of pah range
    wpah_nu_reverse = np.where( (nu_reverse.value < nu_pah2.value) & (nu_reverse.value > nu_pah1.value))[0]




    for i in range(particle_PAH_luminosity.shape[0]):
        #lum = np.trapz(draine_lam[wpah_lam].cgs.value,flam[i,wpah_lam]).value
        fnu_reverse = fnu[i,:][::-1]

        lum = np.absolute(np.trapz(nu_reverse[wpah_nu_reverse].cgs.value,fnu_reverse[wpah_nu_reverse])).value.item()
        #reversing arrays to make nu increasing, and therefore correct for hyperion addition
        m.add_point_source(luminosity = lum,spectrum=(nu_reverse[wpah_nu_reverse].value,fnu_reverse[wpah_nu_reverse].value), position = reg['particle_dust','coordinates'][i,:].in_units('cm').value-boost)
        

    if cfg.par.draine21_pah_grid_write: #else, the try/except in analytics.py will get caught and will just write a single -1 to the output npz file
        reg.parameters['grid_PAH_luminosity'] = grid_PAH_luminosity
    reg.parameters['PAH_lam'] = draine_lam.value

    total_PAH_luminosity =np.sum(grid_PAH_luminosity,axis=0)
    reg.parameters['total_PAH_luminosity'] = total_PAH_luminosity
    
    grid_PAH_L_lam = grid_PAH_luminosity/draine_lam.value
    integrated_grid_PAH_luminosity = np.trapz((grid_PAH_luminosity/draine_lam.value),draine_lam.value,axis=1)
    reg.parameters['integrated_grid_PAH_luminosity'] = integrated_grid_PAH_luminosity
    
    #save some information for dumping into analytics
    reg.parameters['q_pah'] = q_pah
    reg.parameters['particle_mass_weighted_gsd'] = particle_mass_weighted_gsd
    reg.parameters['grid_mass_weighted_gsd'] = grid_mass_weighted_gsd
    reg.parameters['simulation_sizes'] = simulation_sizes


#    for i_cell in tqdm(range(ncells)):
        
        #-----------------------------
        #THIS IS WHERE WE WILL NEED TO EVENTUALLY GET RID OF *FILENAME* UP TOP, AND DECIDE WHAT RADIATION FIELD WE'RE USING FOR THIS PARTICULAR CELL
        #-----------------------------
        
        
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
    


        
