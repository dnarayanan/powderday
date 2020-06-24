Getting Started
**********

Overview of Requirements
============

	* **python>=2.7,>=3.5**

	  * numpy (any version except 1.10.*)
	  * scipy
	  * astropy (3.2.3)
	  * h5py
	  * scikit-learn

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
1. First clone the main repository.::

     >git clone https://github.com/astrofrog/hyperion.git
     
2. Install the python module::

   >cd hyperion
   >pip install .


3. Ensure that if you type::
     >hyperion

it returns a sensible output.  It should return something along the lines of::

  >usage: hyperion [-h] [-f] [-m n_cores] input output
  >hyperion: error: too few arguments

If it can't find `Hyperion <http://www.hyperion-rt.org>`_, check the
the path that is near one of the last lines of the setup.py
installation (that is something associated with the number 755) and
make sure it's in your path.  Ir's most likely to be a python binaries
directory.

4. Install the submodules manually::

   >git submodule init
   >git submodule update

5. Install the Fortran binaries::

     > ./configure
or::

  >./configure --prefix=$HOME/local

or some such path if you aren't administrator on your computer.  Note
for this step you'll need your compilers, MPI and HDF5 installations
active (so, on a supercomputer you might need to load these modules
such as [for example, on the University of Florida HiPerGator
supercomputer])::

  >module load git/2.14.1  intel/2018.1.163  openmpi/3.1.0  libz/1.2.8 hdf5/1.10.1

6. Compile the code::

   > make
   > make install
   


Make sure this works by typing at the command line::

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


yt-4.x configuration 
--------------------

`yt <http://yt-project.org>`_ has recently (June 22nd, 2010)
transitioned from 3.x to 4.x. The latter offers a number of advantages
including a demeshed handling of particle datasets, as well as an
`arepo <https://www.h-its.org/2014/10/28/arepo/>`_ front end.  We are
happy to announce that as of December 31st, 2019 via hash
`59315f311535b5f2309c705f5a71519148aa4f29
<https://github.com/dnarayanan/powderday/commit/59315f311535b5f2309c705f5a71519148aa4f29>`_,
`powderday <https://github.com/dnarayanan/powderday.git>`_ is now `yt
<http://yt-project.org>`_ 4.x compliant.

Unfortunately, there are currently there are some bugs in the mainline
`yt <http://yt-project.org>`_ 4.x branch that are impacting how the
octree is built for particle-like simulations.  These, however, are
fixed in Ashley Kelly's fork (and will be merged into the main `yt
<http://yt-project.org>`_ branch soon enough).  In the mean time, if
you are eager to use `yt <http://yt-project.org>`_ 4.x for particle
simulations, please see below.  


First, it is recommended to make a new python environment in which to run the 
4.x development branch::

    > conda create -n pd4env
    > conda activate pd4env

Next, install `Hyperion <http://www.hyperion-rt.org>`_  as normal, via::

    >conda install -c conda-forge hyperion

And then install `yt <http://yt-project.org>`_ 4.x.  The latter is installed via first ensuring that you have all of the dependencies installed::

  >pip install numpy jupyter sphinx gitpython h5py matplotlib cython nose scipy astropy sympy mpi4py
followed by actually installing  `yt <http://yt-project.org>`_::
  
  >git clone https://github.com/AshKelly/yt.git
  >cd yt
  >git checkout yt-4.0-new-octree
  >git pull
  >pip install -e .

If we do this, this will overwrite the `yt <http://yt-project.org>`_ installation that ships with `Hyperion <http://www.hyperion-rt.org>`_, and you should be good to go.  To check that everything worked, make sure the output of the following 
commands look something like this::

    > ipython
    In [1]: import yt
    In [2]: yt.__version__
    Out[2]: '4.0.dev0'
  
Note, in the above we are actually not yet installing the master
branch of `yt <http://yt-project.org>`_ 4.x, but rather Ashley Kelly's
fork, and the yt-4.0-new-octree branch of Ash's fork.  As this branch gets merged into
the master `yt <http://yt-project.org>`_ 4.x branch, we will update
these docs.


.. _fsps:

fsps
--------------

`fsps <https://code.google.com/p/fsps/source/checkout>`_ can be checked out with::
  
  > git clone https://github.com/cconroy20/fsps

and directions to the installation are in the `Manual <https://www.cfa.harvard.edu/~cconroy/ FSPS_files/MANUAL.pdf>`_.

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

4. `Hyperion <http://www.hyperion-rt.org>`_ does not run with MPI
   support.  Some users have found that when manually installing
   `Hyperion <http://www.hyperion-rt.org>`_ (as opposed to using the
   conda installation) that MPI will then work.

5. On the HiPerGator supercomputing cluster, `Hyperion
   <http://www.hyperion-rt.org>`_ will not install manually, complaining about mismatches with gcc and HDF5.  The solution to this is to use the intel compilers.  The following is a known pacakge combination that works::

  >module load git/2.14.1  intel/2018.1.163  openmpi/3.1.0  libz/1.2.8 hdf5/1.10.1


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
   
