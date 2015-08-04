import yt
import numpy as np
import ipdb
from hyperion.model import Model
from hyperion.grid import AMRGrid
import config as cfg
from plot_generate import proj_plots


    


def enzo_m_gen(fname,field_add):
    

    
    #add the fields in pd format
    pf = field_add(fname)
    ad = pf.all_data()
   
 

    #cutout
    if cfg.par.MANUAL_CENTERING == True:
        center = pf.arr([cfg.model.x_cent,cfg.model.y_cent,cfg.model.z_cent],'code_length')
    else:
        center = pf.domain_center
    
    box_len = pf.quan(cfg.par.zoom_box_len,'kpc').in_units('code_length')
   
    min_region = [center[0]-box_len,center[1]-box_len,center[2]-box_len]
    max_region = [center[0]+box_len,center[1]+box_len,center[2]+box_len]
    region = pf.region(center,min_region,max_region)
  
    pf = region.ds
  
    proj_plots(pf)
    ipdb.set_trace()
    #def. dust density
    def _dust_density(field, data):
        return data[('gas', 'metal_density')].in_units("g/cm**3")*cfg.par.dusttometals_ratio
    
    pf.add_field(('gas', 'dust_density'), function=_dust_density, units = 'g/cm**3')
       
    amr = AMRGrid.from_yt(pf, quantity_mapping={'density':('gas','dust_density')})
    


    '''
    levels = pf.index.max_level
    
    amr = AMRGrid()
    for ilevel in range(levels):
        level = amr.add_level()
        
    for igrid in pf.index.select_grids(ilevel):
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
                
    m.add_density_grid(amr['density'], cfg.par.dustdir+cfg.par.dustfile)
   

    #define the random things needed for parsing out the output args
    #center = pf.domain_center
    [xcent,ycent,zcent] = center
   
    boost = np.array([xcent,ycent,zcent])
    dx = pf.domain_width.in_units('cm')
    dy = pf.domain_width.in_units('cm')
    dz = pf.domain_width.in_units('cm')
    
    
    return m,xcent,ycent,zcent,dx,dy,dz,pf,boost

       
