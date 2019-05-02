Model Tests
**********

Here, we describe some tests of the model.  This growing suite of
tests are packaged with the repository so that developers can examine
the impact of any code updates.



SKIRT tests
============

 `SKIRT <http://www.skirt.ugent.be/root/index.html>`_ is an
 outstanding dust radiative transfer code developed by the Ghent group
 with a number of features that distinguish it from `powderday
 <https://bitbucket.org/desika/powderday>`_.  Here, we describe how to
 run model tests between the two codes.


A few notes:

:Relevant Hash:
   
   The `powderday <https://bitbucket.org/desika/powderday>`_ hash for
   the most recent set of tests is `21204f2
   <https://bitbucket.org/desika/powderday/commits/21204f296c2765322751c3631bb825d46d40ad14>`_
   and the relevant `SKIRT
   <http://www.skirt.ugent.be/root/index.html>`_ hash is
   96e8be9761c3021498b847854bf84c3da3129555.
    

:Simulation Codes:

  The tests done here have only been done for Gadget/Gizmo simulations.

:Downloads:

 Currently, the `SKIRT <http://www.skirt.ugent.be/root/index.html>`_
 and `powderday <https://bitbucket.org/desika/powderday>`_ model
 outputs from our tests are available on the `BitBucket downloads page
 <https://bitbucket.org/desika/powderday/downloads/>`_.  The relevant
 parameter files (the ``.ski`` input file for `SKIRT
 <http://www.skirt.ugent.be/root/index.html>`_ and the
 ``parameters*.py`` files for `powderday
 <https://bitbucket.org/desika/powderday>`_ are packaged with the
 repository.

 
Actually Running the Tests: Powderday:
--------------

To run the `powderday <https://bitbucket.org/desika/powderday>`_ test
code, run::
  python pd_front_end.py tests/SKIRT/mw_zoom parameters_master_401 parameters_model_401

After having edited the ``parameters_master`` and ``parameters_model``
file for your specific paths.  The model snapshot is the same as in
the `Gadget/Gizmo quickstart
<https://powderday.readthedocs.io/en/latest/quickstart.html#gadget-gizmo>`_ section.

The tests run here were done with hash `16f281e <https://bitbucket.org/desika/powderday/commits/16f281e9fa156d7ef0d412a8acbc253bd1aa1389>`_

Actually Running the Tests: SKIRT:
--------------

`SKIRT <http://www.skirt.ugent.be/root/index.html>`_ needs to know
the stellar and gas particle information from the hydrodynamic
simulation.  In ```analytics.py``, there's a function
``SKIRT_data_dump()`` that dumps these files to your  `powderday <https://bitbucket.org/desika/powderday>`_ output directory.  Here, it should have dumped files with a path like::
  tests/SKIRT/mw_zoom/SKIRT.134.gas.particles.txt

`SKIRT <http://www.skirt.ugent.be/root/index.html>`_ also requires an
input parameter file which can be created via a terminal input
following their `tutorials
<http://www.skirt.ugent.be/tutorials/index.html>`_.  We have created the ``.ski`` input file::
  tests/SKIRT/pd_test.dust.ski
  
`SKIRT <http://www.skirt.ugent.be/root/index.html>`_ requires paths
relative to the run directory.  After compiling `SKIRT
<http://www.skirt.ugent.be/root/index.html>`_, one possible setup might be::
  cd SKIRT/run
  mkdir pd_test
  cp <path_to_powderday>/tests/SKIRT/pd_test.dust.ski pd_test
  ../release/SKIRT/main/skirt pd_test/pd_test.dust.ski

In the ``run`` directory, this should create a file like:
``/home/desika.narayanan/SKIRT/run/test_i90_sed.dat``.  This has
chosen one particular viewing angle (which, in fact may not be the
same viewing angle as the `powderday
<https://bitbucket.org/desika/powderday>`_ simulation - in fact it was
arbitrarily chosen).  Then you can run::
  cd <path_to_powderday>
  python tests/SKIRT/skirt_sed_plot.py

and it should produce an image like
..image::images/powderday_skirt_comparison.png (in case the image
doesn't actually embed, this is also found in the repo in the
tests/SKIRT directory).
	   

  


