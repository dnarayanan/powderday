from __future__ import print_function
import numpy as np
from hyperion.model import Model
from hyperion.grid import AMRGrid
import powderday.config as cfg
from powderday.analytics import proj_plots
from powderday.helpers import energy_density_absorbed_by_CMB
from hyperion.dust import SphericalDust
from powderday.grid_construction import enzo_grid_generate
import yt
import os

def enzo_m_gen(fname,field_add):
    
    reg,ds1 = enzo_grid_generate(fname,field_add)

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
    
    #CMB DISABLED -- UNCOMMENT THIS TO FIX THIS.  The main issue is
    #that I'm not sure what shape to give to the np.repeat
    #array of energy_density_absorbed; I think it needs to be the ARM Grid shape but i'm not quite sure if it needs to be an AMRGrid()
    #energy_density_absorbed=energy_density_absorbed_by_CMB()
    #energy_density_absorbed =np.repeat(energy_density_absorbed.value,reg.index.num_grids)#amr['density'].shape)

    d = SphericalDust(cfg.par.dustdir+cfg.par.dustfile)
    if cfg.par.SUBLIMATION == True:
        d.set_sublimation_temperature('fast',temperature=cfg.par.SUBLIMATION_TEMPERATURE)
        
    m.add_density_grid(amr["density"],d)
    #uncomment when we're ready to put CMB in (and comment out previous line)
    #m.add_density_grid(amr['density'],d,specific_energy=energy_density_absorbed)
    #m.set_specific_energy_type('additional')
    
    center = ds1.arr([cfg.model.x_cent,cfg.model.y_cent,cfg.model.z_cent],'code_length')
    [xcent,ycent,zcent] = center.in_units('cm') #boost needs to be in cm since that's what the 
   
    boost = np.array([xcent,ycent,zcent])

    dx = ds1.domain_width.in_units('cm')
    dy = ds1.domain_width.in_units('cm')
    dz = ds1.domain_width.in_units('cm')
    
    return m,xcent,ycent,zcent,dx,dy,dz,reg,ds1,boost

       
