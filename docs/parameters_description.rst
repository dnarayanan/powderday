Detailed Description of Parameters
**********

parameters_master
============

Resolution Keywords
------------

:oref:

   Over Refinement of the Octree - 1 means each data holding cell (a
   False) gets split into 8 one more time.  Very heavy on the memory.
   Default is 0.

:n_ref:
   
   Threshold number of particles to refine over.  When nparticles >
   n_ref the octree refines further.  Default is 64.

:zoom:
   
   Boolean. True means use a specified region of the grid; False means
   use the whole grid.  If set to True, then the grid will have side
   length +/- zoom_box_len

:zoom_box_len:

   Side length (if zoom == True) is +/- zoom_box_len from the center.
   Units are kpc.  So, a grid centered on [0,0,0] with zoom_box_len =
   200 would extend from [-200,200].

:bbox_lim:

   Initial bounding box of grid for SPH simulations (+/- bbox_lim).
   Units are kpc.  This must encompass all of the particles in a
   simulation currently.



Parallelization
------------

:n_processes:

   Number of MPI processes to run the radiative transfer on.  Note,
   the stellar population synthesis will only run on as many
   processors as are on a core since its parallelization is pool.map
   (not MPI)


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

Dust Information
------------

:dustdir:

   path to where your dust files are.  String format -
   (e.g. '/home/desika/hyperion-dust-0.1.0/dust_files/')

:dustfile:
   
   Name of your main dust file.  String format - (e.g. 'd03_3.1_6.0_A.hdf5')

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

   
Hydro Code Units
------------
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

:COSMOFLAG:

   Boolean.  True means this is a cosmological simulation, False means
   idealized galaxy simulation.

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
   live in Charlot & Fall birthclouds (if CF_on = True)

:bulge_stars_age:

   As disk_stars_age but for bulge stars.

:disk_stars_metals:

   Metallicity of disk stars in FSPS metallicity units.  See last page
   of FSPS manual for numbers.  (e.g. 20 = Solar for Padova + BaSeL
   tracks).  Meaningless for cosmological simulations.

:bulge_stars_metals:

   As disk_stars_metals but for bulge stars.
   
   
parameters_model
============
