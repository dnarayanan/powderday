import numpy as np
from powderday.pah.pah_file_read import read_draine_file
from powderday.helpers import find_nearest
from astropy import units as u
import pdb

#to do:

#1. speed up the convolution with grain size distribution

#2. test with a different filename to see if the resolution of the PAH spectra in the MIR changes

#3. think about how to handle the thermal re-radiation in the FIR so
#that it is self consistent with absorption and re-radiation of the
#metagalactic field via hyperion.

#4. put in the radiation fields and appropriate filename per cell

filename = '/home/desika.narayanan/PAHs/vsg_stat_therm.iout'

def pah_source_add(ds,reg,m):

    #first establish the grain size distribution and sizes from the
    #hydro simulation
    grid_of_sizes = ds.parameters['reg_grid_of_sizes']
    simulation_sizes = (ds.parameters['grain_sizes_in_micron']*u.micron).to(u.cm).value
    
    #second, read the information from the Draine files
    PAH_list = read_draine_file(filename)
    draine_sizes = PAH_list[0].size_list
    
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

    
    for i_cell in range(ncells):
        
        #-----------------------------
        #THIS IS WHERE WE WILL NEED TO EVENTUALLY GET RID IF *FILENAME* UP TOP, AND DECIDE WHAT RADIATION FIELD WE'RE USING FOR THIS PARTICULAR CELL
        #-----------------------------
        
        
        for i_size in range(len(simulation_sizes)):

            #while it would be optimal to interpolate, this is a very
            #slow process for every grain size and every wavelength.
            #Given that the Draine sizes are significantly oversampled
            #compared to the typical grain size array of a hydro sim,
            #it's easiest to just assign it to the nearest index
            #value.              

            idx0 = Draine_simulation_idx_left_edge_array[i_size]


            grid_PAH_luminosity[i_cell,:] += PAH_list[idx0].lum * grid_of_sizes[i_cell,i_size]

            total_PAH_luminosity += PAH_list[idx0].lum * grid_of_sizes[i_cell,i_size]
        
            
            
    pdb.set_trace()



    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.loglog(PAH_list[0].lam,total_PAH_luminosity[:]/PAH_list[0].lam)
    ax.set_ylim([1e38,1e43])
    ax.set_xlim([1,1000])
    ax.set_xlabel(r'$\lambda (\mu $m)')
    ax.set_ylabel(r'$L_\lambda$ (erg/s/micron)')
    fig.savefig('/home/desika.narayanan/junk2.png',dpi=300)
    
    pdb.set_trace()
        