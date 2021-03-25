#This file holds the function that adds the active dust grids.  This
#file is called from the tributaries, and contains all the common
#lines of code for adding dust grids once the grid_of_sizes has been
#established.

from __future__ import print_function
import numpy as np
import powderday.config as cfg
from hyperion.dust import SphericalDust
import pdb

def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx



def active_dust_add(m,grid_of_sizes,nsizes,dustdens,specific_energy,refined=[False]):
        #for empty cells, use the median size distribution
        for isize in range(nsizes):
            wzero = np.where(grid_of_sizes[:,isize] == 0)[0]
            wnonzero = np.where(grid_of_sizes[:,isize] != 0)[0]

            grid_of_sizes[wzero,isize] = np.median(grid_of_sizes[wnonzero,isize])

            print(len(wzero)/len(wnonzero))



        #now load the mapping between grain bin and filename for the lookup table
        data = np.load(cfg.par.pd_source_dir+'active_dust/dust_files/binned_dust_sizes.npz')
        grain_size_left_edge_array = data['grain_size_left_edge_array']
        grain_size_right_edge_array  = data['grain_size_right_edge_array']
        dust_filenames = data['outfile_filenames']

        nbins = len(grain_size_left_edge_array)




        #find which sizes in the hydro simulation correspond to the
        #pre-binned extinction law sizes from dust_file_writer.py

        dust_file_to_grain_size_mapping_idx = []
        x=np.linspace(cfg.par.otf_extinction_log_min_size,cfg.par.otf_extinction_log_max_size,nsizes)
        for i in range(nbins):
            dust_file_to_grain_size_mapping_idx.append(find_nearest(x,grain_size_left_edge_array[i]))


        #set up the frac array that is nbins big.  this is the
        #fractional contribution of each dust file bin which is based
        #on the total number of grains in the grid in that bin.

        #frac =np.zeros([dustdens.shape[0],nbins])

        dsf_grid = np.zeros([dustdens.shape[0],nbins])
        frac_grid = np.zeros([dustdens.shape[0],nbins])
        #------------------------
        #DEBUG BLOCK

        if cfg.par.OTF_EXTINCTION_MRN_FORCE == True:
            grid_sum = np.zeros(nbins)

            #how DNSF was set up.  not needed other than for testing
            x=np.linspace(-4,0,41)
            #load an example dust size function for testing against
            dsf = np.loadtxt(cfg.par.pd_source_dir+'active_dust/mrn_dn.txt')#DNSF_example.txt')

            nbins = len(grain_size_left_edge_array)


            for i in range(nbins):
                idx = find_nearest(x,grain_size_left_edge_array[i])
                #this sets the fraction of each bin size we need (for the
                #entire grid!)
                dsf_grid[:,i] = dsf[idx]
                grid_sum[i] = np.sum(dsf_grid[:,i])

            #set up the frac array that is nbins big.  this is the
            #fractional contribution of each dust file bin which is based
            #on the total number of grains in the grid in that bin.
            frac = grid_sum/np.sum(grid_sum)

            #now we need to set the localized extinction law. we do
            #this by comparing, fractionally, a given cell's number of
            #grains in that bin to the maximum number of grains that
            #the grid has in that bin.

            for i in range(nbins):
                frac_grid[:,i] = dsf_grid[:,i]/np.max(dsf_grid[:,i])*frac[i]



        #------------------------    

        else:


            grid_sum = np.zeros(nbins)


            #this sets the fraction of each bin size we need (for the
            #entire grid!)
            for i in range(nbins):
                grid_sum[i] = np.sum(grid_of_sizes[:,dust_file_to_grain_size_mapping_idx[i]])


            #set up the frac array that is nbins big.  this is the
            #fractional contribution of each dust file bin which is based
            #on the total number of grains in the grid in that bin.
            frac = grid_sum/np.sum(grid_sum)

            #now we need to set the localized extinction law. we do
            #this by comparing, fractionally, a given cell's number of
            #grains in that bin to the maximum number of grains that
            #the grid has in that bin.

            if np.sum(refined) > 0:
                wFalse = np.where(np.asarray(refined) == 0)[0]
                
                for i in range(nbins):
                        frac_grid[wFalse,i] = grid_of_sizes[:,dust_file_to_grain_size_mapping_idx[i]]/np.max(grid_of_sizes[:,dust_file_to_grain_size_mapping_idx[i]])*frac[i]
            else:
                frac_grid[:,i] = grid_of_sizes[:,dust_file_to_grain_size_mapping_idx[i]]/np.max(grid_of_sizes[:,dust_file_to_grain_size_mapping_idx[i]])*frac[i]

        #now add the dust grids to hyperion
        for bin in range(nbins):
            file = dust_filenames[bin]
            d = SphericalDust(cfg.par.pd_source_dir+'active_dust/'+file)
            m.add_density_grid(dustdens*frac_grid[:,bin],d,specific_energy=specific_energy)

