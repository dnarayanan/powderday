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
import six

def yt_dataset_to_amr_grid_xyz(ds, quantity_mapping={}):
    """
    Version of Hyperion's yt_dataset_to_amr_grid (from_yt) function for xyz grids
    Convert a yt dataset to a Hyperion AMRGrid object
    .. note:: This function requires yt 3.0 or later
    Parameters
    ----------
    ds : yt Dataset
        The yt dataset
    quantity_mapping : dict
        A dictionary mapping the name of the quantity to use in Hyperion (the
        key) to the name of the field to extract in yt (the value). An example
        is provided below.
    Notes
    -----
    The domain is always re-centered so that the position at
    ds.domain_center in yt becomes the origin in Hyperion.
    """

    field_list = "\n    ".join([str(x) for x in ds.derived_field_list])

    if len(quantity_mapping) == 0:
        raise ValueError("quantity_mapping needs to specified with key:value "
                         "pairs where the key is the name to give the quantity "
                         "in Hyperion and value is the name of the field in the "
                         "yt dataset. Available quantities are: \n\n    {0}".format(field_list))

    for output_quantity, input_field in six.iteritems(quantity_mapping):
        if not isinstance(output_quantity, six.string_types):
            raise ValueError("quantity_mapping keys should be strings")
        if input_field not in ds.derived_field_list:
            raise ValueError("yt field {0} does not exist. Available fields "
                             "are: \n\n    {1}".format(input_field, field_list))

    x0, y0, z0 = ds.domain_center.in_units('cm').ndarray_view() # reordered
    dx, dy, dz = ds.domain_width.in_units('cm').ndarray_view() # reordered

    print("Domain center: x={0}cm, y={1}cm, z={2}cm".format(x0, y0, z0))
    print("Domain width: dx={0}cm, dy={1}cm, dz={2}cm".format(dx, dy, dz))

    # Get levels and limits of all the grids
    n_levels = ds.index.max_level + 1
    levels = ds.index.grid_levels
    # reodered
    xmin, ymin, zmin = ds.index.grid_left_edge.in_units('cm').ndarray_view().transpose()
    xmax, ymax, zmax = ds.index.grid_right_edge.in_units('cm').ndarray_view().transpose()

    print("Re-centering simulation so that domain center is at (0, 0, 0)")
    xmin -= x0
    xmax -= x0
    ymin -= y0
    ymax -= y0
    zmin -= z0
    zmax -= z0

    # Loop over levels and add grids
    amr = AMRGrid()
    for ilevel in range(n_levels):

        # Add a new level
        level = amr.add_level()

        # Loop over yt grids that are at this level
        for igrid in np.nonzero(levels == ilevel)[0]:

            # Get yt grid
            yt_grid = ds.index.grids[igrid]

            # Add a new Hyperion grid
            grid = level.add_grid()

            grid.xmin, grid.xmax = xmin[igrid], xmax[igrid]
            grid.ymin, grid.ymax = ymin[igrid], ymax[igrid]
            grid.zmin, grid.zmax = zmin[igrid], zmax[igrid]

            grid.nx, grid.ny, grid.nz = yt_grid.shape #reordered
            
            # Transposes data as it's being read into hyperion AMR grid
            # Unclear if np.asfortranarray is necessary
            for output_quantity, input_field in six.iteritems(quantity_mapping):
                grid.quantities[output_quantity] = np.asfortranarray(np.transpose(yt_grid[input_field].in_units('g/cm**3').ndarray_view()))

    return amr

def enzo_m_gen(fname,field_add):
    
    reg,ds1 = enzo_grid_generate(fname,field_add)

    amr = yt_dataset_to_amr_grid_xyz(ds1, quantity_mapping={'density':('gas','dust_density')})

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

    dx = ds1.domain_width[0].in_units('cm')
    dy = ds1.domain_width[1].in_units('cm')
    dz = ds1.domain_width[2].in_units('cm')
    
    return m,xcent,ycent,zcent,dx,dy,dz,reg,ds1,boost

