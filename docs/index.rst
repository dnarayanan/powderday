.. Powderday documentation master file, created by
   sphinx-quickstart on Sun Mar  8 13:46:18 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: images/pd_front_image.png
   :align: center


Welcome to Powderday's documentation!
=====================================

`powderday <https://bitbucket.org/desika/powderday>`_) is a dust
radiative transfer package designed to interface with galaxy formation
simulations in order to produce realistic spectral energy
distributions and images. The code utilizes`fsps
<https://code.google.com/p/fsps/source/checkout>`_ (and its python
hooks, `python-fsps
<http://dan.iel.fm/python-fsps/current/installation/>`_) for stellar
SEDs, and `Hyperion <http://www.hyperion-rt.org>`_ for dust radiative
transfer.  Threaded throughout is `yt <http://yt-project.org>`_.

The code has two principle aims - to be flexible (and extremely
modular), and to be easy to use.  The package is written in python,
though makes use of fortran, C and cython.  

As of now, the code works principally for Gadget, Gizmo and other
variants.  Support for Gasoline, RAMSES, ART and Enzo are planned
however.  The current timeline is for TIPSY front ends to be built
over Spring 2015, with AMR support to follow

Finally, this code benefits from the contributions, either directly to
this software, or indirectly by contributions to dependency software,
by many astrophysicists, including Matt Turk, Tom Robitaille, Robert
Thompson, Charlie Conroy, Dan Forman-Mackey and Phil Hopkins.

.. important:: **Disclaimer**:

   While the developers have made every effort to ensure
   that the code is bug-free, bugs invariably find their
   way in.  The developers are not responsible for
   incorrect results, arising either from inherent bugs
   in the code, or mistakes in usage.  Any questions
   regarding the code or usage should be sent to the
   mailing list.



Contents:

.. toctree::
   :maxdepth: 2

   installation.rst
   quickstart.rst
   parameters_description.rst
   algorithm_description.rst
   convenience.rst
   referencing.rst 
   troubleshooting.rst
  

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

