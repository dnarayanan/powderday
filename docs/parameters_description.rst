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
   
   True means use a specified region of the grid; False means use the
   whole grid.  If set to True, then the grid will have side length
   +/- zoom_box_len

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
   
parameters_model
============
