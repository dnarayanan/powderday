from __future__ import print_function
import numpy as np
from hyperion.model import Model
from hyperion.grid import AMRGrid
import powderday.config as cfg
from powderday.analytics import proj_plots
from powderday.helpers import energy_density_absorbed_by_CMB
from hyperion.dust import SphericalDust
import yt


def enzo_m_gen(fname,field_add):
    
    
    #add the fields in pd format
    ds = field_add(fname)
    
    #def. dust density
    def _dust_density(field, data):
        return data[('gas', 'metal_density')].in_units("g/cm**3")*cfg.par.dusttometals_ratio
    ds.add_field(('gas', 'dust_density'), function=_dust_density, units = 'g/cm**3')


    ad = ds.all_data()
   
    #now set up a region DEBUG DEBUG DEBUG THIS IS TOTALLY HACKED AND NOT THE RIGHT REGION
    #box_len = ds.quan(0.1,'code_length')
    #center = ds.domain_center
    #min_region = [center[0]-box_len,center[1]-box_len,center[2]-box_len]
    #max_region = [center[0]+box_len,center[1]+box_len,center[2]+box_len]
    
    #reg = ds.region(center,min_region,max_region)


    #set up the region again
    center = ds.arr([cfg.model.x_cent,cfg.model.y_cent,cfg.model.z_cent],'code_length')
    box_len = ds.quan(cfg.par.zoom_box_len,'kpc').in_units('code_length')
    min_region = [center[0]-box_len,center[1]-box_len,center[2]-box_len]
    max_region = [center[0]+box_len,center[1]+box_len,center[2]+box_len]
    reg = ds.region(center,min_region,max_region)


    reg.save_as_dataset('temp_enzo.h5',fields=[('all','creation_time'),('gas','metal_density'),('gas','density'),('newstars','metallicity_fraction'),('newstars','particle_mass'),('all', 'particle_index'),('index', 'grid_level'),('gas','dust_density')])
    ds1 = yt.load('temp_enzo.h5')
    ad1 = ds1.all_data()
    
    
    ds1.index.get_levels = reg.index.get_levels
    ds1.index.get_smallest_ds = reg.index.get_smallest_dx
    ds1.index.grid = reg.index.grid
    ds1.index.grid_corners = reg.index.grid_corners
    ds1.index.grid_dimensions = reg.index.grid_dimensions
    ds1.index.grid_levels = reg.index.grid_levels
    ds1.index.grid_left_edge = reg.index.grid_left_edge
    ds1.index.grid_right_edge = reg.index.grid_right_edge
    ds1.index.grid_particle_count = reg.index.grid_particle_count
    ds1.index.grids = reg.index.grids
    ds1.index.index_filename = reg.index.index_filename
    ds1.index.max_level = reg.index.max_level
    ds1.index_num_grds = reg.index.num_grids
    ds1.index.num_stars = reg.index.num_stars
    ds1.index.parameters = reg.index.parameters
    
    #ds1 = field_add('temp_enzo.h5',ds=ds1,starages=True)

    amr = AMRGrid.from_yt(ds1, quantity_mapping={'density':('gas','dust_density')})
    


    '''
    levels = ds.index.max_level
    
    amr = AMRGrid()
    for ilevel in range(levels):
        level = amr.add_level()
        
    for igrid in ds.index.select_grids(ilevel):
        print igrid
        grid = level.add_grid()
        grid.xmin,grid.xmax = igrid.LeftEdge[0].in_units('cm'),igrid.RightEdge[0].in_units('cm')
        grid.ymin,grid.ymax = igrid.LeftEdge[1].in_units('cm'),igrid.RightEdge[1].in_units('cm')
        grid.zmin,grid.zmax = igrid.LeftEdge[2].in_units('cm'),igrid.RightEdge[2].in_units('cm')
        grid.quantities["density"] = np.transpose(np.array(igrid[("gas","metal_density")].in_units('g/cm**3')*cfg.par.dusttometals_ratio))
        grid.nx,grid.ny,grid.nz = igrid[("gas","metal_density")].shape
    '''


    m = Model()

    m.set_amr_grid(amr)
    
    #DEBUG CMB DISABLED -- UNCOMMENT THIS TO FIX THIS
    #energy_density_absorbed=energy_density_absorbed_by_CMB()
    #energy_density_absorbed = np.repeat(energy_density_absorbed.value,amr['density'].shape)


    d = SphericalDust(cfg.par.dustdir+cfg.par.dustfile)
    if cfg.par.SUBLIMATION == True:
        d.set_sublimation_temperature('fast',temperature=cfg.par.SUBLIMATION_TEMPERATURE)
    m.add_density_grid(amr['density'],d)#DEBUG CMB DISABLED UNCOMMENT TO FIX THIS,specific_energy=energy_density_absorbed)
    #DEBUG CMB DISABLED UNCOMMENT TO FIX THIS m.set_specific_energy_type('additional')
#m.add_density_grid(amr['density'], cfg.par.dustdir+cfg.par.dustfile)
    

    #define the random things needed for parsing out the output args
    #center = ds.domain_center
    [xcent,ycent,zcent] = center.in_units('cm') #boost needs to be in cm since that's what the 
   
    boost = np.array([xcent,ycent,zcent])

    dx = ds.domain_width.in_units('cm')
    dy = ds.domain_width.in_units('cm')
    dz = ds.domain_width.in_units('cm')
    
    return m,xcent,ycent,zcent,dx,dy,dz,reg,ds,boost

       
