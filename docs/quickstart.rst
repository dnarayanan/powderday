Quickstart
**********

In the examples subdirectory of the `powderday
<https://bitbucket.org/desika/powderday>`_ root directory are some
example snapshots for different hydro codes suppported thusfar.  This
will likely change over time as the code evolves and parameter files
change.  Also, eventually the examples will migrate to the `agora
<https://sites.google.com/site/santacruzcomparisonproject/>`_ project
snapshots.

For each example file there should be two parameter files that will be
reasonable for the associated snapshot, though you'll need to edit the
hard linked directories that specify where (e.g.) dust files are and
output should go.  To run, type (in the `powderday
<https://bitbucket.org/desika/powderday>`_) source directory::

  >python pd_front_end.py <example directory> <parameters_master_file>
  <parameters_model_file>

Note - the .py extensions on the parameter files need to be left off.



Gadget/Gizmo
============
The example simulation is a cosmological zoom simulation of a Milky Way mass galaxy that can be downloaded here (6 GB download):

 <https://www.dropbox.com/s/g6d47z3pm8l18p7/snapshot_134.hdf5?dl=0>

To run the code, you would type::

  >python pd_front_end.py examples/gadget/mw_zoom parameters_master_401 parameters_model_401

The SED (placed at z = 3 with a Planck13 cosmology) looks like:

.. image:: images/gadget_sed.png
   :align: center

and an example plotting code can be found in the convenience
subdirectory of the `powderday
<https://bitbucket.org/desika/powderday>`_ root directory.

Gasoline/Changa
============
