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
   tests/SKIRT/pd_test.dust.ski`.   

 
