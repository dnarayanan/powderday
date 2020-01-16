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
    

    '''
    #add the fields in pd format
    ds = field_add(fname)
    
    #def. dust density
    def _dust_density(field, data):
        return data[('gas', 'metal_density')].in_units("g/cm**3")*cfg.par.dusttometals_ratio

    ds.add_field(('gas', 'dust_density'), function=_dust_density, units = 'g/cm**3')


    ad = ds.all_data()
   
 

    #cutout
    center = ds.arr([cfg.model.x_cent,cfg.model.y_cent,cfg.model.z_cent],'code_length')
    
    box_len = ds.quan(cfg.par.zoom_box_len,'kpc').in_units('code_length')
   
    min_region = [center[0]-box_len,center[1]-box_len,center[2]-box_len]
    max_region = [center[0]+box_len,center[1]+box_len,center[2]+box_len]
    reg = ds.region(center,min_region,max_region)
    
    reg.save_as_dataset('temp_enzo.h5')
    ds1 = yt.load('temp_enzo.h5')
    ad1 = ds1.all_data()
    
    print(ad1[('all', 'creation_time')])
    '''

    ds = yt.load(fname)
    ds = field_add(fname)

    #DEBUG DEBUG DEBUG DEBUG ALL THIS NEEDS TO BE MERGED IN WITH WHATS ABOVE AS ITS SUPER HACKY AND PRETENDS E.G. LIKE WE KNOW WHAT THE BOX LENGTH IS ETC. ACK!!!!

    #def. dust density
    def _dust_density(field, data):
        return data[('gas', 'metal_density')].in_units("g/cm**3")*cfg.par.dusttometals_ratio

    ds.add_field(('gas', 'dust_density'), function=_dust_density, units = 'g/cm**3')


    ad = ds.all_data()
    print(ad[('all', 'creation_time')])
    
    #set up a region
    box_len = ds.quan(0.1,'code_length')
    center = ds.domain_center
    min_region = [center[0]-box_len,center[1]-box_len,center[2]-box_len]
    max_region = [center[0]+box_len,center[1]+box_len,center[2]+box_len]
    
    reg = ds.region(center,min_region,max_region)
    
    import pdb
    pdb.set_trace()
    
    #now save the region to disk as a ds force the addition of some
    #fields that don't otherwise get saved in the save_as_dataset
    #method, because they are derived fields
    reg.save_as_dataset('temp_enzo.h5',fields=[('all','creation_time'),('gas','metal_density')])
    
    #load the region back up
    ds1 = yt.load('temp_enzo.h5')
    ad1 = ds1.all_data()
    print(ad1[('all', 'creation_time')])

    import pdb
    pdb.set_trace()

    ds1 = field_add('temp_enzo.h5',ds=ds1,starages=True)
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

    energy_density_absorbed=energy_density_absorbed_by_CMB()
    energy_density_absorbed = np.repeat(energy_density_absorbed.value,amr['density'].shape)


    d = SphericalDust(cfg.par.dustdir+cfg.par.dustfile)
    if cfg.par.SUBLIMATION == True:
        d.set_sublimation_temperature('fast',temperature=cfg.par.SUBLIMATION_TEMPERATURE)
    m.add_density_grid(amr['density'],d,specific_energy=energy_density_absorbed)
    m.set_specific_energy_type('additional')
 #m.add_density_grid(amr['density'], cfg.par.dustdir+cfg.par.dustfile)
    

    #define the random things needed for parsing out the output args
    #center = ds.domain_center
    [xcent,ycent,zcent] = center
   
    boost = np.array([xcent,ycent,zcent])
    dx = ds.domain_width.in_units('cm')
    dy = ds.domain_width.in_units('cm')
    dz = ds.domain_width.in_units('cm')
    
    
    return m,xcent,ycent,zcent,dx,dy,dz,ds,boost

       
