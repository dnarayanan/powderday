import numpy as np
from powderday.pah.pah_file_read import read_draine_file
from powderday.helpers import find_nearest
from astropy import units as u
import powderday.config as cfg
import pdb
from tqdm import tqdm

#to do:

#1. speed up the convolution with grain size distribution -- currently
#takes ~2 hours per galaxy. Assume the slow part is the np.dot product.

#2. think about how to handle the thermal re-radiation in the FIR so
#that it is self consistent with absorption and re-radiation of the
#metagalactic field via hyperion.

#3. put in the radiation fields and appropriate filename per cell

#4. put in the ionization fraction and appropriate file name per cell.

filename = '/blue/narayanan/desika.narayanan/powderday_files/PAHs/iout_graD16emtPAHib_mmpisrf_1.00'
#filename = '/blue/narayanan/desika.narayanan/powderday_files/PAHs/BC03_Z0.02_10Myr/iout_graD16emtPAHib_bc03_z0.02_1e7_1.00'
#filename = '/blue/narayanan/desika.narayanan/powderday_files/PAHs/BC03_Z0.0004_10Myr/iout_graD16emtPAHib_bc03_z0.0004_1e7_1.50'
#iout_DH20Ad_P0.20_0.00_bc03_z0.02_1e7_1.00'
#filename = '/home/desika.narayanan/PAHs/vsg_stat_therm.iout'

def pah_source_add(ds,reg,m):

    #first establish the grain size distribution and sizes from the
    #hydro simulation
    grid_of_sizes = ds.parameters['reg_grid_of_sizes']
    simulation_sizes = (ds.parameters['grain_sizes_in_micron']*u.micron).to(u.cm).value

    #determine q_PAH for analysis and save it to parameters for writing out
    ad = ds.all_data()
    idx_pah = np.where(simulation_sizes <= 3.e-7)[0]
    dN_pah = np.sum(grid_of_sizes[:,idx_pah],axis=1)
    dN_total = np.sum(grid_of_sizes,axis=1)

    try: #mesh based codes or arepo 
        q_pah = (dN_pah * ad['dust','mass'])/(dN_total* ad['dust','mass'])
        ds.parameters['q_pah'] = q_pah
        #compute the mass weighted grain size distributions for comparison in analytics.py #DEBUG DEBUG DEBUG THIS IS TOO GIZMO CENTRIC
        particle_mass_weighted_gsd = np.average(reg['dust','numgrains'],weights=reg['dust','mass'],axis=0)
        grid_mass_weighted_gsd = np.average(grid_of_sizes,weights=reg['dust','mass'],axis=0)
    except:
        q_pah = (dN_pah * reg['dust','smoothedmasses'])/(dN_total* reg['dust','smoothedmasses'])
        q_pah[np.isnan(q_pah)] = 0 #since we often have octs with no dust
        ds.parameters['q_pah'] = q_pah
        #compute the mass weighted grain size distributions for comparison in analytics.py #DEBUG DEBUG DEBUG THIS IS TOO GIZMO CENTRIC
        particle_mass_weighted_gsd = np.average(reg['dust','numgrains'],weights=reg['dust','mass'],axis=0)
        grid_mass_weighted_gsd = np.average(grid_of_sizes,weights=reg['dust','smoothedmasses'],axis=0)
        

    #second, read the information from the Draine files
    PAH_list = read_draine_file(filename)
    print("reading Draine File",filename)
    draine_sizes = PAH_list[0].size_list
    draine_lam = PAH_list[0].lam
    #third, on a cell-by-cell basis, interpolate the luminosity for
    #each grain size bin, and multiply by the number of grains in that
    #bin

    ncells = grid_of_sizes.shape[0]

    #define the grid that will store the PAH spectra for the entire mesh
    grid_PAH_luminosity = np.zeros([ncells,len(PAH_list[0].lam)])
    total_PAH_luminosity = np.zeros(len(PAH_list[0].lam))

    #find the indices of the Draine sizes that best match those that are in the simulation
    Draine_simulation_idx_left_edge_array = []
    for size in simulation_sizes:
        idx0 = find_nearest(draine_sizes,size)
        #if draine_sizes[idx0] > size: idx0 -=1
        
        #this is really the nearest point in the Draine sizes to the
        #simulation_sizes.  if we uncomment the above line, it will
        #become the left point.
        Draine_simulation_idx_left_edge_array.append(idx0)


    size_arange = np.arange(len(simulation_sizes))

    #in principle, PAH_list will change cell by cell as the radiation
    #field changes cell by cell.  then, this should all go inside the
    #coming up commented loop.  however, this takes forever...like 2
    #hours vs a few sec.  for testing, for now, this is going to
    #remain out of the loop.


    pah_grid = np.array([x.lum for x in PAH_list])
    idx = np.asarray(Draine_simulation_idx_left_edge_array)[size_arange]

    #set the PAH luminosity of the cell to be the dot product of
    #the Draine luminosities (i.e., pah_grid[idx,:] which has
    #dimensions (simulation_sizes,wavelengths)) with the actual
    #grain size distribution in that cell (i.e.,
    #grid_of_sizes[i_cell,:]). note, we take the transpose of
    #grid_of_sizes to get the dimensions to match up correctly for the dot product

    grid_PAH_luminosity = np.dot(pah_grid[idx,:].T, grid_of_sizes.T).T
    reg.parameters = {}

    if cfg.par.draine21_pah_grid_write: #else, the try/except in analytics.py will get caught and will just write a single -1 to the output npz file
        reg.parameters['grid_PAH_luminosity'] = grid_PAH_luminosity
        reg.parameters['PAH_lam'] = draine_lam

    total_PAH_luminosity =np.sum(grid_PAH_luminosity,axis=0)
    reg.parameters['total_PAH_luminosity'] = total_PAH_luminosity
    
    grid_PAH_L_lam = grid_PAH_luminosity/draine_lam
    integrated_grid_PAH_luminosity = np.trapz((grid_PAH_luminosity/draine_lam),draine_lam,axis=1)
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
    


        
