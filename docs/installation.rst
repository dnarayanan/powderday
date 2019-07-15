Getting Started
**********

Overview of Requirements
============

	* **python2.x**

	  * numpy (any version except 1.10.*)
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
	  * git  <http://git-scm.com/>
	  * powderday <http://bitbucket.org/desika/powderday>
	  * yt <http://yt-project.org>
	  * FSPS <https://code.google.com/p/fsps/source/checkout>
	  * python-fsps <http://dan.iel.fm/python-fsps/current/>
	  * Hyperion <http://www.hyperion-rt.org/>
	  * Hyperion Dust Files <http://docs.hyperion-rt.org/en/stable/dust/dust.html>

Installation
============
	    


Manual Installation
--------------

What follows is a self-contained installation manual, though for
problematic installs of any of the sub packages, it's definitely best
to visit the main docs on the main software site (which are always
linked below in each subsection).

.. _python:

python
--------------

`powderday <https://bitbucket.org/desika/powderday>`_ currently only
works with python 2.x with provisional support for python 3.x (i.e. it
*should* work though there may be issues yet outstanding -- please do
file an issue in on the BitBucket site if you find a 3.x issue)).  The
code was developed on, and principally tested with python 2.7.

This said, the code will be robust with python 3.x starting in Summer
of 2019.

As you will see, `powderday <https://bitbucket.org/desika/powderday>`_
currently requires a particular branch of `yt
<http://yt-project.org>`_. As a result, one path that we have seen work
well for users is to set up a different python environment for the
`powderday <https://bitbucket.org/desika/powderday>`_ installation.   This could look something like (assuming a ``conda`` installation of python)::

  >conda create --name pd_environment python=2.7
  >source activate pd_environment

(And then when you want to exit the environment, you can type)::

  >source deactivate pd_environment

Then, whenever you're in the ``pd_environment``, everything you
install will remain contained to that particular installation of
python, and not conflict with any of your other installed packages.

.. _mercurial:


mercurial
--------------

You'll need this to clone `powderday
<https://bitbucket.org/desika/powderday>`_ using mercurial (hg).  If
you don't have mercurial, you can install it in your
``pd_environment`` via::

  >source activate pd_environment
  >conda install mercurial


.. _powderday:


powderday
--------------

Simply clone the latest and greatest from the repo::

  >hg clone https://desika@bitbucket.org/desika/powderday

And that's it!  Once it's cloned, there's no subsidiary installation commands.

.. _yt:

yt
--------------

`yt <http://yt-project.org>`_ 3.x comes bundled with 
`Hyperion <http://www.hyperion-rt.org>`_, so it is not necessary to install it 
on its own. However, as of Summer 2019, development has begun on making
`powderday <https://bitbucket.org/desika/powderday>`_ use
`yt <http://yt-project.org>`_ 4.x, the 'demeshed' 
version of `yt <http://yt-project.org>`_. The 
`powderday <https://bitbucket.org/desika/powderday>`_ -
`yt <http://yt-project.org>`_ 4.x update is in its early stages and is not 
expected to function properly just yet. That being said, development is taking
place on the ``pd-4.x`` branch of 
`powderday <https://bitbucket.org/desika/powderday>`_, and instructions for
setting up  to run with
`yt <http://yt-project.org>`_ 4.x are at the bottom of this section.


.. _fsps:

fsps
--------------

`fsps <https://code.google.com/p/fsps/source/checkout>`_ can be checked out with::
  
  > git clone https://github.com/cconroy20/fsps

and directions to the installation are in the `Manual <https://www.cfa.harvard.edu/~cconroy/ FSPS_files/MANUAL.pdf>`_.

It's likely going to be necessary downstream when installing  `python-fsps
<http://dan.iel.fm/python-fsps/current/installation/>`_ to have the -fPIC flags set in `fsps <https://code.google.com/p/fsps/source/checkout>`_ when making.  So, in the Makefile of `fsps <https://code.google.com/p/fsps/source/checkout>`_, set::
  
  >F90FLAGS = -O -cpp -fPIC

if your ``gcc`` version is lower than 4.3.0, or::

  >F90FLAGS = -03 -mtune=native -cpp -fPIC

if ``gcc`` is version 4.3.0 or higher. This can be checked with 
``gcc --version``. Additionally, at this time 
`powderday <https://bitbucket.org/desika/powderday>`_  doesn't work with the 
default MIST Isochrones.  To fix this, you'll need to edit sps_vars.f90 in 
`fsps <https://code.google.com/p/fsps/source/checkout>`_  to look like::
  
  !------set the isochrone library------!
  #define MIST 0
  !Padova models circa 2008
  #define PADOVA 1
  #define PARSEC 0
  #define BASTI 0
  #define GENEVA 0

To explicitly compile::

  make clean
  make
  
Finally, the SPS_HOME variable must be set in your environment to point to the FSPS/src directory.  For example, if your environment is bash, in your .bashrc set something along the lines of::
   
  >export SPS_HOME=/Users/desika/fsps/



.. _python-fsps:

python-fsps
--------------

`powderday <https://bitbucket.org/desika/powderday>`_ depends on
python hooks for `fsps
<https://code.google.com/p/fsps/source/checkout>`_ written by Daniel
Foreman-Mackey and others called `python-fsps
<http://dan.iel.fm/python-fsps/current/installation/>`_.  You can install from the GitHub page::
  
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
directions for installation are well-described on the main
`Installation page for Hyperion
<http://docs.hyperion-rt.org/en/stable/installation/installation.html>`_.
Here, we summarize the installation which should get most users
through without any real difficulty.

There are two ways to install `Hyperion <http://www.hyperion-rt.org>`_.  The first is via ``conda``::

  >conda install -c conda-forge hyperion

Note, this will eventually become deprecated for `powderday
<https://bitbucket.org/desika/powderday>`_ (or at least modified as
the `Hyperion <http://www.hyperion-rt.org>`_ ``conda`` install ships
with `yt 3.x <http://yt-project.org>`_, and eventual upgrade to `yt
4.x <http://yt-project.org>`_ is coming in Summer 2019.

The second and manual way to install `Hyperion
<http://www.hyperion-rt.org>`_ follows:
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

  >tar -xvzf hyperion-dust-xxx.tar.gz
  >cd hyperion-dust-0.1.0
  >python setup.py build_dust

If you want to use the PAH model in `powderday
<https://bitbucket.org/desika/powderday>`_, you'll additionally need
these files in the same dust directory.  To download, click on the link,
then click 'raw' on the right side of each page.

1. https://github.com/hyperion-rt/paper-galaxy-rt-model/blob/master/dust/big.hdf5
2. https://github.com/hyperion-rt/paper-galaxy-rt-model/blob/master/dust/vsg.hdf5
3. https://github.com/hyperion-rt/paper-galaxy-rt-model/blob/master/dust/usg.hdf5

Please note the caveat that the PAH files are generated using some
approxmations described in `Robitaille et
al. <http://www.aanda.org/articles/aa/abs/2012/09/aa19073-12/aa19073-12.html>`_,
and we encourage the user of these PAH files to read this paper,
especially section 3.4.2.


yt-4.x configuration [WIP]
--------------------

First, it is recommended to make a new python environment in which to run the 
4.x development branch::

    > conda create -n pd4env python=2.7
    > conda activate pd4env

Since `Hyperion <http://www.hyperion-rt.org>`_ comes with 
`yt <http://yt-project.org>`_ 3.x, it must be installed using the ``--no-deps``
flag, since you will install the dependencies manually in the next step::

    > conda install --no-deps -c conda-forge hyperion

Now, install all of the dependencies `Hyperion <http://www.hyperion-rt.org>`_
needs, *except* `yt <http://yt-project.org>`_::

    > conda install -c conda-forge astropy atomicwrites attrs backports backports.functools_lru_cache backports.shutil_get_terminal_size backports_abc configparser contextlib2 cycler cython dbus decorator enum34 expat fastcache fontconfig freetype funcsigs functools32 futures gettext glib gmp gmpy2 gst-plugins-base gstreamer h5py hdf5 hyperion-fortran icu importlib_metadata ipython ipython_genutils jpeg kiwisolver libblas libcblas libgfortran-ng libiconv liblapack libpng libuuid libxcb libxml2 linecache2 matplotlib more-itertools mpc mpfr mpi mpich mpmath numpy openblas packaging pathlib2 pcre pexpect pickleshare pluggy prompt_toolkit pthread-stubs ptyprocess py pygments pyparsing pyqt pytest python-dateutil pytz qt scandir simplegeneric singledispatch sip six subprocess32 sympy tornado traceback2 traitlets unittest2 wcwidth xorg-libxau xorg-libxdmcp xz zipp

Now, clone the 4.x development branch from the `yt <http://yt-project.org>`_
source repository and build it::

    > git clone https://github.com/yt-project/yt.git
    > cd yt/
    > git checkout yt-4.0
    > python setup.py develop

Now, to check that everything worked, make sure the output of the following 
commands look something like this::

    > which yt
    ~/miniconda3/envs/myenv/bin/yt

    > ipython
    In [1]: import yt
    In [2]: yt.__version__
    Out[2]: '4.0.dev0'

As long as the rest of `powderday <https://bitbucket.org/desika/powderday>`_'s
dependencies have been installed, at this point you should be good to go.


Troubleshooting your Installation
============

  .. _python-fsps installation issues:

python-fsps installation issues
--------------

1.  `python-fsps
<http://dan.iel.fm/python-fsps/current/installation/>`_ can't find f2py
   
   f2py is a numpy package that is sometimes named f2py2.7 by numpy.
   At the same time, `python-fsps
   <http://dan.iel.fm/python-fsps/current/installation/>`_ expects it
   to be called f2py (as it sometimes is; for example in Anaconda).
   So, you might need to locate f2py (it ships with `yt
   <http://yt-project.org>`_, so if you for example use the `yt
   <http://yt-project.org>`_ python) you need to link the following
   files::

   >cd /Users/desika/yt-x86_64/bin
   >ln -s f2py2.7 f2py

   and::

   >cd /Users/desika/yt-x86_64/lib/python2.7/site-packages
   >ln -s numpy/f2py/ f2py

   This should hopefully fix it.


2. Issues with 'f2py' in the  `python-fsps
   <http://dan.iel.fm/python-fsps/current/installation/>`_ installation:

   Numpy has made some changes to f2py in the 1.10.x version of numpy.
   The easiest fix is to use a non 1.10.* version of numpy (thanks to
   Ben Johnson for finding this).

3.  `python-fsps
<http://dan.iel.fm/python-fsps/current/installation/>`_ has mysterious
installation failures.  Often this has to do with a bad `FSPS
<https://github.com/cconroy20/fsps>`_ compilation. Even if it seems
like `FSPS <https://github.com/cconroy20/fsps>`_ has compiled, it may
not actually execute properly if the correct compilers aren't set in
the MakeFile.  Thanks to Ena Choi for pointing this one out.

Hyperion Installation Issues
---------------

1. Manual installations seem to not be fully updated from the
   `Hyperion <http://www.hyperion-rt.org>`_ website.  The following
   issues are known (uncovered by Katarina Kraljic)

   a. Hyperion-0.9.10 does not contain the /deps/fortran directory. It
      will be necesary to take this from version 0.9.9

   b. /deps/fortran/install.py hardcodes some links that do not exist
      anymore.  The URLs should be updated as:

      ZLIB_URL = "http://zlib.net/zlib-1.2.11.tar.gz"  
      HDF5_URL = 'http://www.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.gz'

   c. Hyperions configure file doesn't have an option for an MPI
    compiler that is mpif90.openmpi.  One option is to add this to the configure file around line 1940::

      if test "$mpi_compiler" == mpif90.openmpi
         then
            mpi_compiler=`basename $(mpif90 -show | awk {'print $1'})`
      fi
   
