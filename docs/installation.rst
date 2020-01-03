Getting Started
**********

Overview of Requirements
============

	* **python>=2.7,>=3.5**

	  * numpy (any version except 1.10.*)
	  * scipy
	  * astropy
	  * h5py

	* **compilers**

	  * gcc
	  * gfortran


	* **Additional Packages (with Instructions Below)**

	  * git  <http://git-scm.com/>
	  * powderday <https://github.com/dnarayanan/powderday.git>
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

`powderday <https://github.com/dnarayanan/powderday.git>`_ should work with python >=2.7 or >=3.5, but has not been heavily used under >=3.5.
Please file an issue if you encounter one.

As you will see, `powderday <https://github.com/dnarayanan/powderday.git>`_
currently requires a particular branch of `yt
<http://yt-project.org>`_. As a result, one path that we have seen work
well for users is to set up a different python environment for the
`powderday <https://github.com/dnarayanan/powderday.git>`_ installation.   This could look something like (assuming a ``conda`` installation of python)::

  >conda create --name pd_environment
  >source activate pd_environment

(And then when you want to exit the environment, you can type)::

  >source deactivate pd_environment

Then, whenever you're in the ``pd_environment``, everything you
install will remain contained to that particular installation of
python, and not conflict with any of your other installed packages.

.. _powderday:


powderday
--------------

Simply clone the latest and greatest from the repo::

  >git clone https://github.com/dnarayanan/powderday.git

To install, `cd` into the cloned repository and run the usual `python setup.py install`.

.. _yt:

yt
--------------

`yt <http://yt-project.org>`_ 3.x comes bundled with 
`Hyperion <http://www.hyperion-rt.org>`_, so it is not necessary to install it 
on its own. However, as of Summer 2019, development has begun on making
`powderday <https://github.com/dnarayanan/powderday.git>`_ use
`yt <http://yt-project.org>`_ 4.x, the 'demeshed' 
version of `yt <http://yt-project.org>`_. The 
`powderday <https://github.com/dnarayanan/powderday.git>`_ -
`yt <http://yt-project.org>`_ 4.x update is in its early stages and is not 
expected to function properly just yet. That being said, development is taking
place on the ``pd-4.x`` branch of 
`powderday <https://github.com/dnarayanan/powderday.git>`_, and instructions for
setting it up to run with
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
`powderday <https://github.com/dnarayanan/powderday.git>`_  doesn't work with the 
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

`python-fsps <https://github.com/dfm/python-fsps>`_  will be installed automatically by the `powderday` setup.py script.
  
You can test the installation by opening python and typing::

>import fsps

.. _Hyperion:

Hyperion
--------------

`Hyperion <http://www.hyperion-rt.org>`_ is the main work horse of
`powderday <https://github.com/dnarayanan/powderday.git>`_.  The full
directions for installation are well-described on the main
`Installation page for Hyperion
<http://docs.hyperion-rt.org/en/stable/installation/installation.html>`_.
Here, we summarize the installation which should get most users
through without any real difficulty.

There are two ways to install `Hyperion <http://www.hyperion-rt.org>`_.  The first is via ``conda``::

  >conda install -c conda-forge hyperion

Note, this will eventually become deprecated for `powderday
<https://github.com/dnarayanan/powderday.git>`_ (or at least modified as
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
<https://github.com/dnarayanan/powderday.git>`_, you'll additionally need
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

In the future, `yt <http://yt-project.org>`_ will eventually
transition from 3.x to 4.x. The latter offers a number of advantages
including a demeshed handling of particle datasets, as well as an
`arepo <https://www.h-its.org/2014/10/28/arepo/>`_ front end.  We are
happy to announce that as of December 31st, 2019 via hash
`59315f311535b5f2309c705f5a71519148aa4f29
<https://github.com/dnarayanan/powderday/commit/59315f311535b5f2309c705f5a71519148aa4f29>`_,
`powderday <https://github.com/dnarayanan/powderday.git>`_ is now `yt
<http://yt-project.org>`_ 4.x compliant.

The following offer documentation for installing `yt
<http://yt-project.org>`_ 4.x, as well as a slightly different order
of operations for `powderday
<https://github.com/dnarayanan/powderday.git>`_ code dependency
installation.  Please note, `yt <http://yt-project.org>`_ 4.x is still
in development mode, and should be treated as such.  The following
documentation utilizes not-yet-merged in branches of `yt
<http://yt-project.organological>`_ for example, and all 4.x features
should be considered experimental currently.

First, it is recommended to make a new python environment in which to run the 
4.x development branch::

    > conda create -n pd4env
    > conda activate pd4env

Our entire goal here is to install `powderday <https://github.com/dnarayanan/powderday.git>`_  that is dependent on `yt <http://yt-project.org>`_ 4.x.  Since `Hyperion <http://www.hyperion-rt.org>`_ comes with 
`yt <http://yt-project.org>`_ 3.x,  we a few options, neither of which are particularly elegant.

The easiest way to handle this is to install `Hyperion <http://www.hyperion-rt.org>`_  as normal, via::

    >conda install -c conda-forge hyperion

And then install `yt <http://yt-project.org>`_ 4.x.  The latter is installed via first ensuring that you have all of the dependencies installed::

  >pip install numpy jupyter sphinx gitpython h5py matplotlib cython nose scipy astropy sympy mpi4py
followed by actually installing  `yt <http://yt-project.org>`_::
  
  >git clone https://github.com/AshKelly/yt.git
  >cd yt
  >git checkout yt-4.0-new-octree
  >git pull
  >pip install -e .

If we do this, this will overwrite the `yt <http://yt-project.org>`_ 3.x installation that ships with `Hyperion <http://www.hyperion-rt.org>`_, and you should be good to go.  To check that everything worked, make sure the output of the following 
commands look something like this::

    > ipython
    In [1]: import yt
    In [2]: yt.__version__
    Out[2]: '4.0.dev0'
  
Note, in the above we are actually not yet installing the master
branch of `yt <http://yt-project.org>`_ 4.x, but rather Ashley Kelly's
branch which is still experimental.  As Ash's branch gets merged into
the master `yt <http://yt-project.org>`_ 4.x branch, we will update
these docs.

The second way of handling this, which is no longer recommended as it
can be a bit painful and with varying success rates, is to install `Hyperion
<http://www.hyperion-rt.org>`_ with the ``--no-deps`` flag, since you
will install the dependencies manually in the next step::

    > conda install --no-deps -c conda-forge hyperion

Now, install all of the dependencies `Hyperion <http://www.hyperion-rt.org>`_
needs, *except* `yt <http://yt-project.org>`_::

    > conda install -c conda-forge astropy atomicwrites attrs backports backports.functools_lru_cache backports.shutil_get_terminal_size backports_abc configparser contextlib2 cycler cython dbus decorator enum34 expat fastcache fontconfig freetype funcsigs functools32 futures gettext glib gmp gmpy2 gst-plugins-base gstreamer h5py hdf5 hyperion-fortran icu importlib_metadata ipython ipython_genutils jpeg kiwisolver libblas libcblas libgfortran-ng libiconv liblapack libpng libuuid libxcb libxml2 linecache2 matplotlib more-itertools mpc mpfr mpi mpich mpmath numpy openblas packaging pathlib2 pcre pexpect pickleshare pluggy prompt_toolkit pthread-stubs ptyprocess py pygments pyparsing pyqt pytest python-dateutil pytz qt scandir simplegeneric singledispatch sip six subprocess32 sympy tornado traceback2 traitlets unittest2 wcwidth xorg-libxau xorg-libxdmcp xz zipp

At this point, test your `Hyperion <http://www.hyperion-rt.org>`_ installation by trying_::

  import hyperion

If this doesn't work, repeat:_::

  > conda install --no-deps -c conda-forge hyperion

Now, install clone the 4.x development branch from the `yt <http://yt-project.org>`_ 4.x as above.

As long as the rest of `powderday <https://github.com/dnarayanan/powderday.git>`_ 's
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
   
