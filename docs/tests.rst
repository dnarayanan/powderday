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
 <https://bitbucket.org/desika/powderday>`_.  ``analytics.py`` dumps out
 the stellar and gas information via a function ``SKIRT_data_dump``.
 This will dump files similar to, for example::
   tests/SKIRT/mw_zoom/SKIRT.134.gas.particles.txt

 In::
   tests/SKIRT/mw_zoom

   you'll find the ``parameters_master`` and
 ``parameters_model`` parameter files to successfully run the a
 `powderday <https://bitbucket.org/desika/powderday>`_ run performed
 here.  The tests run here were done with hash `16f281e
 <https://bitbucket.org/desika/powderday/commits/16f281e9fa156d7ef0d412a8acbc253bd1aa1389>`_
 using the ``.ski`` file::
   tests/SKIRT/pd_test.dust.ski

Currently, the `SKIRT <http://www.skirt.ugent.be/root/index.html>`_
tests have only been developed for ``Gadget/Gizmo`` simulations.
   	    
Actually Running the Tests
--------------


Powderday:
--------------

To run the `powderday <https://bitbucket.org/desika/powderday>`_ test
code, run::
  python pd_front_end.py tests/SKIRT/mw_zoom parameters_master_401 parameters_model_401

After having edited the ``parameters_master`` and ``parameters_model``
file for your specific paths.  The model snapshot is the same as in
the `Gadget/Gizmo quickstart
<https://powderday.readthedocs.io/en/latest/quickstart.html#gadget-gizmo>`_ section.



SKIRT:
--------------
