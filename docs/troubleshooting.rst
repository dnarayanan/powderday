Troubleshooting and Known Issues
**********

1. yt compilation complains with something like: "Error compiling Cython file"

Answer - try running::

pip install -U Cython

2. Errors along the lines of::


     compiling Fortran sources
     Fortran f77 compiler: /usr/bin/gfortran -Wall -ffixed-form -fno-second-underscore -fPIC -O3 -funroll-loops
     Fortran f90 compiler: /usr/bin/gfortran -fPIC -fPIC -O3 -funroll-loops

Try re-compiling with the flag::

  -fPIC

3. 
