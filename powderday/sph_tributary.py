from __future__ import print_function
import numpy as np
import yt

from hyperion.model import Model
import matplotlib as mpl
mpl.use('Agg')
import powderday.config as cfg

from powderday.grid_construction import yt_octree_generate
from powderday.find_order import find_order
import powderday.powderday_test_octree as pto
import powderday.hyperion_octree_stats as hos

from hyperion.dust import SphericalDust

from powderday.helpers import energy_density_absorbed_by_CMB
from powderday.analytics import dump_cell_info

def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx


#only comes to be necessary if cfg.par.otf_extinction is set
def _dust_density(field,data):
    return data.ds.arr(data[('PartType3','Dust_Density')],'code_mass/code_length**3')

def _size_with_units(field,data):
    return data.ds.parameters['size']

def sph_m_gen(fname,field_add):
    
    refined,dustdens,fc1,fw1,reg,ds = yt_octree_generate(fname,field_add)
    
    if yt.__version__ == '4.0.dev0':
        xmin = (fc1[:,0]-fw1[:,0]/2.).to('cm') #in proper cm 
        xmax = (fc1[:,0]+fw1[:,0]/2.).to('cm')
        ymin = (fc1[:,1]-fw1[:,1]/2.).to('cm')
        ymax = (fc1[:,1]+fw1[:,1]/2.).to('cm')
        zmin = (fc1[:,2]-fw1[:,2]/2.).to('cm')
        zmax = (fc1[:,2]+fw1[:,2]/2.).to('cm')
    else:
        xmin = (fc1[:,0]-fw1[:,0]/2.).convert_to_units('cm') #in proper cm
        xmax = (fc1[:,0]+fw1[:,0]/2.).convert_to_units('cm')
        ymin = (fc1[:,1]-fw1[:,1]/2.).convert_to_units('cm')
        ymax = (fc1[:,1]+fw1[:,1]/2.).convert_to_units('cm')
        zmin = (fc1[:,2]-fw1[:,2]/2.).convert_to_units('cm')
        zmax = (fc1[:,2]+fw1[:,2]/2.).convert_to_units('cm')

    #dx,dy,dz are the edges of the parent grid
    dx = (np.max(xmax)-np.min(xmin)).value
    dy = (np.max(ymax)-np.min(ymin)).value
    dz = (np.max(zmax)-np.min(zmin)).value


    xcent = float(ds.quan(cfg.model.x_cent,"code_length").to('cm').value)
    ycent = float(ds.quan(cfg.model.y_cent,"code_length").to('cm').value)
    zcent = float(ds.quan(cfg.model.z_cent,"code_length").to('cm').value)

    boost = np.array([xcent,ycent,zcent])
    print ('[sph_tributary] boost = ',boost)
    print ('[sph_tributary] xmin (pc)= ',np.min(xmin.to('pc')))
    print ('[sph_tributary] xmax (pc)= ',np.max(xmax.to('pc')))
    print ('[sph_tributary] ymin (pc)= ',np.min(ymin.to('pc')))
    print ('[sph_tributary] ymax (pc)= ',np.max(ymax.to('pc')))
    print ('[sph_tributary] zmin (pc)= ',np.min(zmin.to('pc')))
    print ('[sph_tributary] zmax (pc)= ',np.max(zmax.to('pc')))
    #Tom Robitaille's conversion from z-first ordering (yt's default) to
    #x-first ordering (the script should work both ways)



    refined_array = np.array(refined)
    refined_array = np.squeeze(refined_array)
    
    order = find_order(refined_array)
    refined_reordered = []
    dustdens_reordered = np.zeros(len(order))
    


    
    for i in range(len(order)): 
        refined_reordered.append(refined[order[i]])
        dustdens_reordered[i] = dustdens[order[i]]


    refined = refined_reordered
    dustdens=dustdens_reordered

    #hyperion octree stats
    max_level = hos.hyperion_octree_stats(refined)


    pto.test_octree(refined,max_level)

    dump_cell_info(refined,fc1,fw1,xmin,xmax,ymin,ymax,zmin,zmax)
    np.save('refined.npy',refined)
    np.save('density.npy',dustdens)
    

    #========================================================================
    #Initialize Hyperion Model
    #========================================================================

    m = Model()
    


    print ('Setting Octree Grid with Parameters: ')



    #m.set_octree_grid(xcent,ycent,zcent,
    #                  dx,dy,dz,refined)
    m.set_octree_grid(0,0,0,dx/2,dy/2,dz/2,refined)    

    #get CMB:
    
    energy_density_absorbed=energy_density_absorbed_by_CMB()
    specific_energy = np.repeat(energy_density_absorbed.value,dustdens.shape)


    ''' #debug and remove this after testing is done'''
    if cfg.par.otf_extinction == False:
        
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
          #entire galaxy, we'll compute it on a cell-by-cell basis by
          #using information about the grain size distribution from
          #the simulation itself.



        print("==============================================\n")
        print("Entering OTF Extinction Calculation\n")
        print("Note: For very high-resolution grids, this may cause memory issues due to adding ncells dust grids")
        print("==============================================\n")



        '''
        #THIS HAS BEEN MOVED TO GADGET2PD - hopefully can stay there and this can all be deleted.


        #We're assuming that the particle type that the active dust is
        #in in PartType3 and adding it to the sph types so that it can
        #be deposited onto the octree
        ds._sph_ptypes = ('PartType0','PartType3')
        
        #add the dust density field
        ds.add_field(('PartType3','density'),function=_dust_density,units='code_mass/code_length**3',sampling_type='particle',particle_type=True)
        ad = ds.all_data()
        nsizes = ad['PartType3','Dust_Size'].shape[1] #number of dust size bins

        #now loop through th esize bins and project them each onto
        #their own octree.  we'll then, after projecting them into
        #individual grids, collate those grids back together into a
        #master size grid

        for isize in range(nsizes):
            
            ds.parameters['size'] = 0 #just to clear it out
            
            #we have to slice it before we add the field.  if you try to slice
            #in the field defintion (i.e., in _size_with_units), yt freaks
            ds.parameters['size'] = ad['PartType3','Dust_Size'][:,isize]
            
            #actually add the sliced field now.  we do this so that we can
            #prepare for depositing onto the octree. we call this a dummy size
            #since this is just there for a place holder to deposit into a dummy octree
            
            print('adding and depositing fields for dust size bin '+str(isize))
            ds.add_field(('PartType3','dummy_size_bin'+str(isize)),function=_size_with_units,sampling_type='particle',units='dimensionless',particle_type=True,force_override=True)
            ad = ds.all_data()

            #deposit onto the octree.   this is what will get merged into the final octree
            if yt.__version__ != '4.0.dev0':
                ds.add_deposited_particle_field(('PartType3','dummy_size_bin'+str(isize)),"sum")
        '''

        '''
        EVERYTHING BEFORE THIS MAY NOT BE NEEDED IF ALL GOES WELL BECAUSE
        THIS IS INCLUDED IN THE FRONT END. THIS SAID, AFTER THIS WE
        NEED TO NOW ACCESS THESE DEPOSITED OCTREES FROM REG AND SUM THEM INTO THE FINAL OCTREE_OF_SIZES
        '''
        #WHERE WE'RE LEAVING OFF: AT THE END OF GADGET2PD, DOING
        #LEN(DUM_SIZE_OCTREE) gives 512 (i.e., the size of the refined
        #grid coming out of the reg region).  what we need to do is:
        #(a) ensure that the dust octree HERE is the same size, anda
        #(b) once we do, see if it has any information in it.  if it
        #doesn't, then we need to put in a try/except like that
        #commented below to crash the code.

        ad = ds.all_data()
        nsizes = ad['PartType3','Dust_Size'].shape[1] #number of dust size bins
        for isize in range(nsizes):

            if isize == 0:
                if  yt.__version__ != '4.0.dev0':
                    octree_of_sizes = np.zeros((ad[('deposit','PartType3_sum_dummy_size_bin'+str(isize))].shape[0],nsizes))
                else:
                    octree = ds.octree()
                    octree_of_sizes = np.zeros((octree[('PartType3', 'dummy_size_bin'+str(isize))].shape[0],nsizes))

            
            if  yt.__version__ != '4.0.dev0':
                octree_of_sizes[:,isize] = ad[('deposit','PartType3_sum_dummy_size_bin'+str(isize))]
            else:
                octree_of_sizes[:,isize] = octree[('PartType3', 'dummy_size_bin'+str(isize))]

        octree_of_sizes[np.isnan(octree_of_sizes)] = 0 #just because the density can be zero in some extreme cases which screws up the particle deposition
        import pdb
        pdb.set_trace()
        
        #                try:
        #            np.sum(dum_size_octree) > 0
         #       except ValueError:
         #           raise ValueError("[front_ends/gadget2pd:] The grain size distribution smoothed onto the octree has deposited no particles.  Try either increasing your box size, or decreasing n_ref in parameters_master")


        #now load the mapping between grain bin and filename for the lookup table
        data = np.load(cfg.par.pd_source_dir+'active_dust/dust_files/binned_dust_sizes.npz')
        grain_size_left_edge_array = data['grain_size_left_edge_array']
        grain_size_right_edge_array  = data['grain_size_right_edge_array']
        dust_filenames = data['outfile_filenames']

        nbins = len(grain_size_left_edge_array)
        frac = np.empty(nbins)

        #find which sizes in the hydro simulation correspond to the
        #pre-binned extinction law sizes from dust_file_writer.py
        x=np.linspace(cfg.par.otf_extinction_log_min_size,cfg.par.otf_extinction_log_max_size,nsizes)
        for i in range(nbins):
            idx = find_nearest(x,grain_size_left_edge_array[i])

        
        #now loop through the cells in the octree, and create a frac *
        #dustdens for each cell density grid, and add that density
        #grid.  this means we add ncells density grids NOTE: THIS
        #COULD END UP BEING A MEMORY PROBLEM FOR VERY HIGH-RES SIMS.
        
        for cell in range(len(dustdens)):
            frac = octree_of_sizes[cell]
            import pdb
            if np.sum(frac > 0):
                print(cell)
                print(frac)

        

        #this sets the fraction of each bin size we need (for the
        #entire grid! this eventually needs to be cell by cell)

        


        import pdb
        pdb.set_trace()
        #STOPPING HERE: NEED TO FIGURE OUT HOW TO SET THIS RELATIVE RATIO OF GRAINS IN THIS SIZE FOR EVERY CELL
        frac[i] = dsf[idx]

        

    '''
    #TESTING BLOCK ONLY - MUCH OF THIS NEEDS TO BE REMOVED
    '''

    '''
    #load the mapping between grain bin and filename
    data = np.load(cfg.par.pd_source_dir+'active_dust/dust_files/binned_dust_sizes.npz')
    grain_size_left_edge_array = data['grain_size_left_edge_array']
    grain_size_right_edge_array  = data['grain_size_right_edge_array']
    dust_filenames = data['outfile_filenames']


    #how DNSF was set up.  not needed other than for testing
    x=np.linspace(-4,0,41)
    #load an example dust size function for testing against
    dsf = np.loadtxt(cfg.par.pd_source_dir+'active_dust/mrn_dn.txt')#DNSF_example.txt')

    nbins = len(grain_size_left_edge_array)
    frac = np.empty(nbins)
    for i in range(nbins):
        idx = find_nearest(x,grain_size_left_edge_array[i])
        
        #this sets the fraction of each bin size we need (for the
        #entire grid! this eventually needs to be cell by cell)
        frac[i] = dsf[idx]

    # Normalize to 1
    total = np.sum(frac)
    frac/=total
    
    for counter,file in enumerate(dust_filenames):
        d = SphericalDust(cfg.par.pd_source_dir+'active_dust/'+file)
        m.add_density_grid(dustdens*frac[counter],d,specific_energy=specific_energy)

    '''
    
    





    m.set_specific_energy_type('additional')








    return m,xcent,ycent,zcent,dx,dy,dz,reg,ds,boost
