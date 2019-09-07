Detailed Description of Parameters
**********

There are two parameters files - parameters_master and
parameters_model.  The distinction between where a parameter goes is
somewhat arbitrary, but is based on the following.  Parameters_master
tends to contain parmaeters that will likely be shared amongst all
snapshots in a given galaxy run, while parameters_model are parameters
that might change from run to run (like the snapshot name or the
galaxy center).

The parameters files are (somewhat non-traditionally) written as
python code files because this enables you to embed little snippets in
them that may be particular to your galaxy run.  For example, for
gadget snapshot naming conventions, it can be useful to have a snippet along the lines of::

  snapshot_num = 20
  snapnum_str = '{:03d}'.format(snapshot_num)

  snapshot_name = 'snapshot_'+snapnum_str+'.hdf5'


That ensures a 3 digit snapshot number, common to many gadget-style
simulations.


parameters_master
============

Resolution Keywords
------------

:oref:

   Over Refinement of the Octree.  For particle-based codes, 1 means
   each data holding cell (a False) gets refined one additional time,
   even after the octree refinement criteria has stopped.  Very heavy
   on the memory but can enable higher pixel resolution for images.
   Default is 0.

:n_ref:
   
   Refinement criteria for octree refinement for particle-based codes.
   This is the threshold number of particles to refine over.  When
   nparticles > n_ref the octree refines further.  Default is 64.

:zoom_box_len:

   Side length to zoom in on.  Is +/- zoom_box_len from the center.
   Units are proper kpc.  So, a grid centered on [0,0,0] with
   zoom_box_len = 200 would extend from [-200,200] kpc in physical
   units at the redshift of the simulation.

:bbox_lim:

   Initial bounding box of grid for particle simulations
   (+/- bbox_lim).  Units are kpc.  This must encompass all of the
   particles in a simulation currently.  This just has to be a big
   number, but you want to be careful of making *too* large as
   precision limitations only allow for up to 20 levels of refinement.



Parallelization
------------

:n_processes:

   Number of pool processes to run the radiative transfer on.  Note,
   the stellar population synthesis will only run on as many
   processors as are on a core since its parallelization is pool.map
   (not MPI)

:n_MPI_processes:

    Number of MPI tasks to run. For TORQUE this is best set as the same 
    as n_processes, while for SLURM this may not be the case.


RT Information
------------

For all photon counts, a decent rule of thumb is 10-100x the number of
grid cells that you have, though of course you should check the
convergence properties of your simulation.

:n_photons_initial:

   Number of photons to use in main iterations (for the whole grid)
   for specific energy and temperature calculations.

:n_photons_imaging:

   Number of photons to use for the SED/image calculation

:n_photons_raytracing_sources:

   If raytracing is set (which is the default hard-coded into the
   code), number of raytracing photons to use for source emission.

:n_photons_raytracing_dust:

   Similar to n_photons_raytracing_sources but for dust emission.

:FORCE_RANDOM_SEED:

    Boolean. True means the seed specified below will be used for random number
    generation in the Hyperion model. False means no seed will be set.

:seed:

    The seed with which to initialize random number generation in Hyperion. 
    Must be a negative integer.

Dust Information
------------

:dustdir:

   String. Path to where your dust files are.  String format -
   (e.g. '/home/desika/hyperion-dust-0.1.0/dust_files/')

:dustfile:
   
   String. Name of your main dust file.  String format -
   (e.g. 'd03_3.1_6.0_A.hdf5')

:PAH:

   Boolean - True means use model for PAHs, False means don't.

:dusttometals_ratio:

   Dust mass to metals mass ratio

:enforce_energy_range:

   Boolean. False ensures energy conservation.  But the emisivities
   may not be strictly correct if the energy in a cell is out of range
   of the emissivities.  True modifies the energy in the simulation,
   but ensures that the emissivities are consistent with the energy.
   See:
   <http://docs.hyperion-rt.org/en/latest/api/hyperion.model.Model.html?highlight=enforce#hyperion.model.Model.set_enforce_energy_range>

:SUBLIMATION:

    Boolean. If true, dust grains above the sublimation 
    temperature are automatically killed. Fast mode is to set this to False.

:SUBLIMATION_TEMPERATURE:

    The temperature in Kelvin above which dust grains are automatically killed. 
    Meaningless if SUBLIMATION == False.
   

Hydro Code Units
------------

Currently these are actually not used in `powderday
<https://github.com/dnarayanan/powderday.git>`_).  They remain in the
parameters file as a placeholder though as we may need them as an
over-ride if we find some HDF5 files don't contain this information.

:unit_mass:

   Mass code units for galaxy simulation.  Units: Msun/h

:unit_length:

   Length code unit for galaxy simulation.  Units: kpc/h

:unit_age:

   Stellar age units.  Units: Gyr/h

:unit_velocity:

   Velocity code unit for galaxy simulation.  Units: cm/s


Stellar SEDs Info
------------

:Force_Binning:

   Boolean.  True means force binning of the stellar SEDs (in bins of
   age and metallicity).  False means don't.  False results in an
   exact solution since the stellar SEDs are individually represented
   (as opposed to broken up into bins).  This said, this can be very
   slow to run, and extremely hard on the memory.

:imf_type:

   IMF parameter for stellar pops calculations.

   0. Salpeter
   1. Chabrier
   2. Kroupa
   3. Van Dokkum
   4. Dave

   Though note options 3 and 4 are currently not supported.

:pagb:

   Weight given to post AGB stars.  1 is the default.


Nebular Emission Info
------------

:add_neb_emission:

    Boolean. If set to True, nebular line emission from Cloudy lookup tables 
    (dev. by Nell Byler) will be added.

:add_agb_dust_model:




:CF_on:

   Boolean.  If set to True, then enables the Charlot & Fall
   birthcloud models for all stars with age younger than
   birth_cloud_clearing_age.

:birth_cloud_clearing_age:

   Stars with age < birth_cloud_clearing_age have Charlot & Fall
   birthclouds (if CF_on == True).  Meaningless if CF_on == False.
   Units: Gyr.

:Z_init:

   Forced metallicity increase in the newstar particles.  Useful for
   idealized galaxy simulations where the stars can form out of
   pristine gas.  Units are absolute (so 0.02 = Solar). Setting to 0
   (default) means that you use the stellar metallicities as they come
   in the simulation (i.e. for Cosmological simulations).

:disk_stars_age:

   Age in Gyr of disk stars for idealized simulations. Meaningless for
   cosmological simulations.  Note, if this is <=7, then these will
   live in Charlot & Fall birthclouds (if CF_on = True).

   Note, for Gadget simulations, stars are divided into newstars, disk
   stars and bulge stars.  For Tipsy outputs, the stars initalized
   with the simulation are auto-detected by their nonsensical ages,
   and assigned as disk stars.  So, if there are stars initalized with
   your Tipsy simulation, assign their ages (and metallicities below)
   as disk stars.

:bulge_stars_age:

   As disk_stars_age but for bulge stars.

:disk_stars_metals:

   Metallicity of disk stars in FSPS metallicity units.  See last page
   of FSPS manual for numbers.  (e.g. 20 = Solar for Padova + BaSeL
   tracks).  Meaningless for cosmological simulations.

:bulge_stars_metals:

   As disk_stars_metals but for bulge stars.

:N_STELLAR_AGE_BINS:
   
   Number of bins to bin the stellar ages in (boundaries are the
   oldest and youngest star particles; linear bins in log(age)).

:N_MASS_BINS:

   Meaningless parameter; place holder for future code additions.

:metallicity_legend:

   String.  Location of the metallicity maps in FSPS for the stellar
   libraries you use.  Currently Padova2007 is the default (hard coded
   into `powderday <https://github.com/dnarayanan/powderday.git>`_), so
   this should point to something like:
   "/Users/desika/fsps/ISOCHRONES/Padova/Padova2007/zlegend_basel.dat"
   

Black Holes
------------

:BH_SED:

    If true, `powderday <https://github.com/dnarayanan/powderday.git>`_ will 
    attempt to load black hole information from the snapshot.

:BH_eta:

    Used in calculating the black hole luminosity (bhluminosity = 
    BH_eta * mdot * c**2.)

:BH_model:

    BH model type, either Nenkova or other.

:BH_modelfile:

    The path to the Nenkova model file if BH_model is set to Nenkova. This file
    can be downloaded here and placed anywhere in the repository, as long as 
    the correct path is set in ``parameters_master``: 
    <https://www.clumpy.org/downloads/clumpy_models_201410_tvavg.hdf5>

:nenkova_params:

    Nenkova+ (2008) model parameters.


Images and SED Parameters
------------

:NTHETA:

   Number of polar angles to view galaxy at

:IMAGING:

   Must be set to ``True`` for 
   `powderday <https://github.com/dnarayanan/powderday.git>`_ to produce an image
   output file.

:filterdir:

   Directory where filter files are stored. They should be located in
   "/home/desika/powderday/filters/".

:filterfiles:

   A list of the names of all filters to be used. 
   `powderday <https://github.com/dnarayanan/powderday.git>`_ will run at each 
   wavelength in all specified filter files, and will produce a ``.hdf5`` file
   containing images convolved with each filter transmission function. Note 
   that this can be quite computationally intensive and scales with the number
   of wavelengths. Following the example in ``parameters_master``, additional 
   filters can be added to this list. In bash, ``cd`` into your ``filterdir`` 
   and use the following command to format the filenames for easy copying and 
   pasting into this list.
   
    .. code-block:: bash

       >>> shopt -s globstar; printf "#    '%s'\n" *.filter

:IMAGING_TRANSMISSION_FILTER:

   If enabled, filter convolution will be performed through 
   `Hyperion <http://www.hyperion-rt.org>`_ instead of through `powderday 
   <https://github.com/dnarayanan/powderday.git>`_. This is much faster, but is 
   still an experimental feature and does not seem to produce accurate 
   convolved images.


DEBUGGING
------------

You should probably never touch any of these.


parameters_model
============

:snapshot_name:

   String - currently the snapshot name of your galaxy run. (Naming
   will change as other front ends built).

:hydro_dir:

   Location of snapshots

:PD_output_dir:

   String - location of where `powderday
   <https://github.com/dnarayanan/powderday.git>`_ output files should go.

:Auto_TF_file:

   String - name of the TF logical file to be written (doesn't need a
   path - will go into PD_output_dir)

:Auto_dustdens_file:

   String - name of the dust density ascii file to be written (doesn't
   need a path - will go into PD_output_dir)

:inputfile:

   String - name of the input HDF5 (rtin) file for `powderday
   <https://github.com/dnarayanan/powderday.git>`_ to write before
   radiative transfer begins.

:outputfile:

   String - name of the output HDF5 (rtout) file after radiative transfer

:x_cent:

   Location in grid coordinates of the x-coordinate of the center of
   your galaxy.  Only pertinenet if MANUAL_CENTERING==True.  Otherwise
   ignored by `powderday <https://github.com/dnarayanan/powderday.git>`_.

:y_cent:

   As x_cent but for the y-coordinate

:z_cent:

   As x_cent but for the z-coordinate.
