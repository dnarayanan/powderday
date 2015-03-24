Getting Started
**********

Overview of Requirements
============

	* **python2.x**

	  * numpy
	  * scipy
	  * astropy
	  * Atpy
	  * h5py
	  * ipdb

	* **compilers**

	  * gcc
	  * gfortran


	* **Additional Packages (with Instructions Below)**
	  * mercurial  <http://mercurial.selenic.com/>
	  * powderday <http://bitbucket.org/desika/powderday>
	  * yt <http://yt-project.org>
	  * FSPS <https://code.google.com/p/fsps/source/checkout>
	  * python-fsps <http://dan.iel.fm/python-fsps/current/>
	  * Hyperion <http://www.hyperion-rt.org/>
	  * Hyperion Dust Files <http://docs.hyperion-rt.org/en/stable/dust/dust.html>


Installation
============

We're working on an all-in-one installer for `powderday
<https://bitbucket.org/desika/powderday>`_.  For the time being, for
better or worse, the installation of individual packages that
`powderday <https://bitbucket.org/desika/powderday>`_ depends on is
manual.  What follows is a self-contained installation manual, though
for problematic installs of any of the sub packages, it's definitely
best to visit the main docs on the main software site (which are
always linked below in each subsection).

.. _python:

python
--------------

`powderday <https://bitbucket.org/desika/powderday>`_ currently only
works with python 2.x.  The code was developed on, and principally
tested with python 2.7. 

.. _mercurial:


mercurial
--------------

You'll need this to clone `powderday
<https://bitbucket.org/desika/powderday>`_ using mercurial (hg).  If
you don't want to install mercurial, then first install `yt
<http://yt-project.org>`_ (before `powderday
<https://bitbucket.org/desika/powderday>`_ ) - `yt
<http://yt-project.org>`_ ships with its own hg which you can
optionally use.


.. _powderday:


powderday
--------------

Simply clone the latest and greatest from the repo::

  >hg clone https://desika@bitbucket.org/desika/powderday

And that's it!  Once it's cloned, there's no subsidiary installation commands.

.. _yt:

yt
--------------

`powderday <https://bitbucket.org/desika/powderday>`_ has `yt
<http://yt-project.org>`_ threaded throughout, and thus needs the
software to be installed to function.  If you already have `yt
<http://yt-project.org>`_ installed, be sure that it's `yt
<http://yt-project.org>`_ 3.x (i.e at least the 'stable' branch).

If you don't have `yt <http://yt-project.org>`_ already installed,
it's super easy!  There are many ways to do this (as directed on the
`yt <http://yt-project.org>`_ project website.  We require the
standalone installer script for the Development version.  This will
ensure that bugs caught in or bleeding-edge features of `yt
<http://yt-project.org>`_ are included.  

.. _fsps:

fsps
--------------

`fsps <https://code.google.com/p/fsps/source/checkout>`_ can be checked out with::
  
  >svn checkout http://fsps.googlecode.com/svn/trunk/ fsps

and directions to the installation are in the `Manual <https://www.cfa.harvard.edu/~cconroy/ FSPS_files/MANUAL.pdf>`_

Currently, `python-fsps
<http://dan.iel.fm/python-fsps/current/installation/>`_ is known to
work with revision 145 of `fsps
<https://code.google.com/p/fsps/source/checkout>`_ , so it is safest
to revert to this version via::
  
  >svn update -r 145

in the src directory of fsps.  It is almost certainly going to be necessary downstream when installing  `python-fsps
<http://dan.iel.fm/python-fsps/current/installation/>`_ to have the -fPIC flags set in `fsps <https://code.google.com/p/fsps/source/checkout>`_ when making.  So, in the Makefile of `fsps <https://code.google.com/p/fsps/source/checkout>`_ , set::
  
  >F90FLAGS = -O -cpp -fPIC

Finally, the SPS_HOME variable must be set in your environment to point to the FSPS/src directory.  For example, if your environment is bash, in your .bashrc set something along the lines of::
   
  >export SPS_HOME=/Users/desika/fsps/



.. _python-fsps:

python-fsps
--------------

`powderday <https://bitbucket.org/desika/powderday>`_ depends on
python hooks for `fsps
<https://code.google.com/p/fsps/source/checkout>`_ written by Daniel
Foreman-Mackey and others called `python-fsps
<http://dan.iel.fm/python-fsps/current/installation/>`_.  There are a
few ways to install it.  Perhaps the easiest is via a pip installer::
  >pip install fsps
  
Though you could also install the development version::
  
  >git clone https://github.com/dfm/python-fsps.git
  >cd python-fsps
  >python setup.py install

You can test the installation by opening python and typing::

>import fsps

.. _Hyperion:

Hyperion
--------------

`Hyperion <http://www.hyperion-rt.org>`_ is the main work horse of
`powderday <https://bitbucket.org/desika/powderday>`_.  The full
directions for installation are well-described on the main `Hyperion
<http://www.hyperion-rt.org>`_ website.  Here, we summarize the
installation which should get most users through without any real
difficulty.


1. First download the tarball and unpack it.::

     >tar -xzvf hyperion.xxx
     >cd hyperion.xxx
     
2. Install the fortran dependencies::

   >cd deps/fortran
   >python install.py <prefix>

where <prefix> is where you want the libraries to be installed.  To
avoid conflicts with other packages, I usually install somewhere
like::

  >python install.py /usr/local/hyperion

as suggested by the `Hyperion <http://www.hyperion-rt.org>`_ docs.  Ensure that the
following commands return something sensible::

  >which mpif90
  >which h5fc

if not, your path probably needs to include wherever the <prefix> directory pointed to.
  

 
3. Install any remaining python dependencies. These are listed `here <http://docs.hyperion-rt.org/en/stable/installation/python_dependencies.html>`_  
   
4. Install `Hyperion <http://www.hyperion-rt.org>`_  itself.  To do this::
     
     >cd hyperion.xxx
     >python setup.py install

or::

  >python setup.py install --user

if you don't have root access.  At this point::

  >import hyperion

from within python should work, and typing::

  >hyperion

at the command line should return something along the lines of::

  >usage: hyperion [-h] [-f] [-m n_cores] input output
  >hyperion: error: too few arguments

if not, check the the path that is near one of the last lines of the
setup.py installation (that is something associated with the
number 755) and make sure it's in your path.  Ir's most likely to be a
python binaries directory.

You then have to install the Fortran Binaries::

  >./configure  --prefix=prefix
  >make
  >make install

where the prefix is wherever you installed the Fortran libraries
before.  Make sure this works by typing at the command line::

  >hyperion_sph

which should return something like::

  >Usage: hyperion_sph [-f] input_file output_file


  .. _Hyperion_dust:

Hyperion Dust Files
--------------

Unless you've written your own dust files, you will likely want to use
the pre-compiled dust files developed by Tom Robitaille (though don't
ship with `Hyperion <http://www.hyperion-rt.org>`_ due to their size).
To install these download them here:
http://docs.hyperion-rt.org/en/stable/dust/dust.html.  Then to
install::

  >tar xvzf hyperion-dust-xxx.tar.gz
  >cd hyperion-dust-0.1.0
  >python setup.py build_dust

If you want to use the PAH model in `powderday
<https://bitbucket.org/desika/powderday>`_, you'll additionally need
these files in the same dust directory.  To download, click on the link,
then click 'raw' on the right side of each page.

1. https://github.com/hyperion-rt/paper-galaxy-rt-model/blob/master/dust/big.hdf5
2. https://github.com/hyperion-rt/paper-galaxy-rt-model/blob/master/dust/vsg.hdf5
3. https://github.com/hyperion-rt/paper-galaxy-rt-model/blob/master/dust/usg.hdf5
