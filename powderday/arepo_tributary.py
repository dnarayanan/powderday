from __future__ import print_function
import numpy as np
import yt

from hyperion.model import Model
import matplotlib as mpl
mpl.use('Agg')
import powderday.config as cfg

from powderday.grid_construction import arepo_vornoi_grid_generate


from hyperion.dust import SphericalDust

from powderday.helpers import energy_density_absorbed_by_CMB

def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx


def arepo_m_gen(fname,field_add):
    
    reg,ds,dustdens = arepo_vornoi_grid_generate(fname,field_add)

    xcent = ds.quan(cfg.model.x_cent,'code_length').to('cm') #proper cm
    ycent = ds.quan(cfg.model.y_cent,'code_length').to('cm')
    zcent = ds.quan(cfg.model.z_cent,'code_length').to('cm')
    
    boost = np.array([xcent,ycent,zcent])
    print ('[arepo_tributary/vornoi_m_gen]:  boost = ',boost)

    #========================================================================
    #Initialize Hyperion Model
    #========================================================================

    m = Model()
    


    #because we boost the stars to a [0,0,0] coordinate center, we
    #want to make sure our vornoi tesslation is created in the same manner.
    
    particle_x = reg["gascoordinates"][:,0].to('cm')
    particle_y = reg["gascoordinates"][:,1].to('cm')
    particle_z = reg["gascoordinates"][:,2].to('cm')


    #just for the sake of symmetry, pass on a dx,dy,dz since it can be
    #used optionally downstream in other functions.
    dx = 2.* ds.quan(cfg.par.zoom_box_len,'kpc').to('cm')
    dy = 2.* ds.quan(cfg.par.zoom_box_len,'kpc').to('cm')
    dz = 2.* ds.quan(cfg.par.zoom_box_len,'kpc').to('cm')

    print ('[arepo_tributary] boost = ',boost)
    print ('[arepo_tributary] xmin (pc)= ',(xcent-dx/2.).to('pc'))
    print ('[arepo_tributary] xmax (pc)= ',(xcent+dx/2.).to('pc'))
    print ('[arepo_tributary] ymin (pc)= ',(ycent-dy/2.).to('pc'))
    print ('[arepo_tributary] ymax (pc)= ',(ycent+dy/2.).to('pc'))
    print ('[arepo_tributary] zmin (pc)= ',(zcent-dz/2.).to('pc'))
    print ('[arepo_tributary] zmax (pc)= ',(zcent+dz/2.).to('pc'))

    x_pos_boost = (particle_x-xcent).to('cm')
    y_pos_boost = (particle_y-ycent).to('cm')
    z_pos_boost = (particle_z-zcent).to('cm')
    
    m.set_voronoi_grid(x_pos_boost.value, y_pos_boost.value, z_pos_boost.value)

    #get CMB:
    
    energy_density_absorbed=energy_density_absorbed_by_CMB()
    specific_energy = np.repeat(energy_density_absorbed.value,dustdens.shape)

    if cfg.par.otf_extinction==False:

        if cfg.par.PAH == True:
        
            # load PAH fractions for usg, vsg, and big (grain sizes)
            frac = cfg.par.PAH_frac

            # Normalize to 1
            total = np.sum(list(frac.values()))
            frac = {k: v / total for k, v in frac.items()}

            for size in frac.keys():
                d = SphericalDust(cfg.par.dustdir+'%s.hdf5'%size)
                if cfg.par.SUBLIMATION == True:
                    d.set_sublimation_temperature('fast',temperature=cfg.par.SUBLIMATION_TEMPERATURE)
                    #m.add_density_grid(dustdens * frac[size], cfg.par.dustdir+'%s.hdf5' % size)
                m.add_density_grid(dustdens*frac[size],d,specific_energy=specific_energy)
            m.set_enforce_energy_range(cfg.par.enforce_energy_range)
        else:
            d = SphericalDust(cfg.par.dustdir+cfg.par.dustfile)
            if cfg.par.SUBLIMATION == True:
                d.set_sublimation_temperature('fast',temperature=cfg.par.SUBLIMATION_TEMPERATURE)
            m.add_density_grid(dustdens,d,specific_energy=specific_energy)
        #m.add_density_grid(dustdens,cfg.par.dustdir+cfg.par.dustfile)  

    else: #instead of using a constant extinction law across the
          #entire galaxy, we'll compute it on a cell-by-cell bassis by
          #using information about the grain size distribution from
          #the simulation itself.

        ad = ds.all_data()
        nsizes = reg['PartType0','NumGrains'].shape[1]
        try:
            assert(np.sum(ad['PartType0','NumGrains']) > 0)
        except AssertionError:
            raise AssertionError("[arepo_tributary:] There are no dust grains in this simulation.  This can sometimes happen in an early snapshot of a simulation where the dust has not yet had time to form.")

        #for empty cells, use the median size distribution
        for isize in range(nsizes):
            wzero = np.where(reg['PartType0','NumGrains'][:,isize] == 0)[0]
            wnonzero = np.where(reg['PartType0','NumGrains'][:,isize] != 0)[0]
            
            reg['PartType0','NumGrains'][wzero,isize] = np.median(reg['PartType0','NumGrains'][wnonzero,isize])

            
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
                grid_sum[i] = np.sum(reg['PartType0','NumGrains'][:,dust_file_to_grain_size_mapping_idx[i]])


            #set up the frac array that is nbins big.  this is the
            #fractional contribution of each dust file bin which is based
            #on the total number of grains in the grid in that bin.
            frac = grid_sum/np.sum(grid_sum)

            #now we need to set the localized extinction law. we do
            #this by comparing, fractionally, a given cell's number of
            #grains in that bin to the maximum number of grains that
            #the grid has in that bin.

            

            for i in range(nbins):

                frac_grid[:,i] = reg['PartType0','NumGrains'][:,dust_file_to_grain_size_mapping_idx[i]]/np.max(reg['PartType0','NumGrains'][:,dust_file_to_grain_size_mapping_idx[i]])*frac[i]


        #now add the dust grids to hyperion
        for bin in range(nbins):
            file = dust_filenames[bin]
            d = SphericalDust(cfg.par.pd_source_dir+'active_dust/'+file)
            m.add_density_grid(dustdens*frac_grid[:,bin],d,specific_energy=specific_energy)


    m.set_specific_energy_type('additional')



    
    return m,xcent,ycent,zcent,dx.value,dy.value,dz.value,reg,ds,boost
