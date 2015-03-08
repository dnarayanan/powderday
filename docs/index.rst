.. Powderday documentation master file, created by
   sphinx-quickstart on Sun Mar  8 13:46:18 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: images/pd_front_image.png
   :align: center


Welcome to Powderday's documentation!
=====================================

Powderday is a dust radiative transfer package designed to interface
with galaxy formation simulations in order to produce realistic
spectral energy distributions and images. The code utilizes ''FSPS''
(and its python hooks, ''python-fsps'') for stellar SEDs, and
''Hyperion'' for dust radiative transfer.  Threaded throughout is
''yt''.

The code has two principle aims - to be flexible (and extremely
modular), and to be easy to use.  The package is written in python,
though makes use of fortran, C and cython.  

As of now, the code works principally for Gadget, Gizmo and other
variants.  Support for Gasoline, RAMSES, ART and Enzo are planned
however.  The current timeline is for TIPSY front ends to be built
over Spring 2015, with AMR support to follow

Finally, this code benefits from the contributions, either directly to
this software, or indirectly by contributions to dependency software,
by many astrophysicists:

Matthew Turk (NCSA)

Tom Robitaille (MPIA)

Robert Thompson (University of Western Cape)

Philip Hopkins (Caltech)

Charlie Conroy (Harvard)

Daniel Foreman-Mackey (NYU)

Contents:

.. toctree::
   :maxdepth: 2



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

