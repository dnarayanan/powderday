Troubleshooting and Known Issues
**********

1. yt compilation complains with something like: "Error compiling Cython file"
--------------

Answer - try running::

  pip install -U Cython

2. Python-fsps compilation returns errors along the lines of
--------------
Error::

     compiling Fortran sources
     Fortran f77 compiler: /usr/bin/gfortran -Wall -ffixed-form -fno-second-underscore -fPIC -O3 -funroll-loops
     Fortran f90 compiler: /usr/bin/gfortran -fPIC -fPIC -O3 -funroll-loops

Try re-compiling with the flag::

  -fPIC

3. When running pd via a SLURM scheduler, you get the error when importing fsps
--------------
Error::

   build/bdist.linux-x86_64/egg/fsps/__init__.py in <module>()
   ImportError: Your FSPS version does not seem to be under git version control. FSPS is available on github at https://github.com/cconroy20/fsps and should be cloned from there

Comment out the lines in python-fsps/fsps/__init__.py surrounding the
checking of githashes.  (h/t to 'Ena Choi <http://www.physics.rutgers.edu/~enachoi/EC/Ena_Choi.html>'_ for uncovering this one)
