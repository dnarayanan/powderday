Getting Started
**********

Overview of Requirements
============

	* **python**

	  * numpy
	  * scipy
	  * astropy
	  * Atpy
	  
	* **yt** <http:yt-project.org>
	* **FSPS** <https://code.google.com/p/fsps/source/checkout>
	* **python-fsps** <http://dan.iel.fm/python-fsps/current/>
	* **Hyperion** <http://www.hyperion-rt.org/>

	  
Installation
============

We're working on an all-in-one installer for `powderday
<https://bitbucket.org/desika/powderday>`_.  For the time being, for
better or worse, the installation of individual packagesa that
`powderday <https://bitbucket.org/desika/powderday>`_ depends on is
manual. Here, we'll detail the order of operations, though it's always
wise to check the parent site for installation.

.. _yt:

yt
--------------

First and foremost, you need to have `yt <http://yt-project.org>`_
installed.  There are many ways to do this (as directed on the `yt
<http://yt-project.org>`_ project website.  Utilizing `anaconda
<https://store.continuum.io/>`_) is certainly a good option
because (a) it comes with all of the required python subpackages,
and (b) the `yt <http://yt-project.org>`_ installation is
exceptionally easy then.

Note, that you will need yt 3.x, meaning at least the stable or
development versions.

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
   
  >export SPS\_HOME=/Users/desika/fsps/



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
