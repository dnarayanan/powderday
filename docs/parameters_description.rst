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

   

parameters_model
============
