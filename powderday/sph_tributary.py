#where we're leaving off: make the normal part of the code that does
#not assume MRN for the fractional distribution of grain sizes match
#the mrn part (i.e, put in the frac_grid stuff etc.)

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
from powderday.tributary_dust_add import active_dust_add


#only comes to be necessary if cfg.par.otf_extinction is set
def _dust_density(field,data):
    return data.ds.arr(data[('PartType3','Dust_Density')],'code_mass/code_length**3')

def _size_with_units(field,data):
    return data.ds.parameters['size']

def sph_m_gen(fname,field_add):

    refined,dustdens,fc1,fw1,reg,ds = yt_octree_generate(fname,field_add)
    
    if float(yt.__version__[0:3]) >= 4:
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

    #save some information that can be used in the PAH model compute
    #an effective 'size' of a cell by density = mass/volume and assume
    #spherical geometry.  similarly, saving the particle location information
    if float(yt.__version__[0:3]) >= 4:
        dump_cell_info(refined,fc1.to('cm'),fw1.to('cm'),xmin,xmax,ymin,ymax,zmin,zmax)
    else:
        dump_cell_info(refined,fc1.convert_to_units('cm'),fw1.convert_to_units('cm'),xmin,xmax,ymin,ymax,zmin,zmax)
    reg.parameters['cell_size']=fw1.convert_to_units('cm') #so that we can have a uniform naming scheme for different front ends for saving in analytics/dump_data(
    reg.parameters['cell_position'] = fc1

    np.save('refined.npy',refined)
    np.save('density.npy',dustdens)
    

    #save some information that can be used in the PAH model
    reg.parameters['fw1'] = fw1.convert_to_units('cm')

    #========================================================================
    #Initialize Hyperion Model
    #========================================================================

    m = Model()
    
    #save in the m__dict__ that we're in an oct geometry
    m.__dict__['grid_type']='oct'

    print ('Setting Octree Grid with Parameters: ')



    #m.set_octree_grid(xcent,ycent,zcent,
    #                  dx,dy,dz,refined)
    m.set_octree_grid(0,0,0,dx/2,dy/2,dz/2,refined)    

    #get CMB:
    
    energy_density_absorbed=energy_density_absorbed_by_CMB()
    specific_energy = np.repeat(energy_density_absorbed.value,dustdens.shape)


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
        
        ad = ds.all_data()
        nsizes = ad['PartType3','Dust_Size'].shape[1]
        ncells = reg.parameters["octree_of_sizes"].shape[0]
        #ensure that the grid has particles
        for isize in range(nsizes):
            try:
                assert (np.sum(reg.parameters["octree_of_sizes"][:,isize]) > 0)
            except AssertionError:
                raise AssertionError("[sph_tributary:] The grain size distribution smoothed onto the octree has deposited no particles.  Try either increasing your box size, or decreasing n_ref in parameters_master.  Alternatively, run the simulation with otf_extinction=False")

        #define the grid of sizes that will be used in tributary_dust_add
        grid_of_sizes = reg.parameters["octree_of_sizes"]
        
        #for now, we just include these grids as -1 because there's no
        #gizmo model [that we know of] that delineates different
        #species.  when this occurs, we should fix this section,
        #modeled after arepo_tributary.
        grid_of_sizes_graphite = [-1]
        grid_of_sizes_silicates = [-1]
        grid_of_sizes_aromatic_fraction = [-1]
        


        active_dust_add(ds,m,grid_of_sizes,nsizes,dustdens,specific_energy,refined,grid_of_sizes_graphite,grid_of_sizes_silicates,grid_of_sizes_aromatic_fraction)

    m.set_specific_energy_type('additional')

    return m,xcent,ycent,zcent,dx,dy,dz,reg,ds,boost
