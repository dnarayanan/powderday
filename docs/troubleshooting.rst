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
checking of githashes.  (h/t to `Ena Choi <http://www.physics.rutgers.edu/~enachoi/EC/Ena_Choi.html>`_ for uncovering this one).

4. If in powderday you start getting memory errors on the Pool process call that look like this:
--------------
Error::

  calculating the SEDs for  42  bins
  Traceback (most recent call last):
  File "/ufrc/narayanan/desika.narayanan/pd/pd_front_end.py", line 139, in <module>
    m=add_binned_seds(df_nu,stars_list,diskstars_list,bulgestars_list,m)
  File "/ufrc/narayanan/desika.narayanan/pd/source_creation.py", line 244, in add_binned_seds
    #of star particles that go in every [wz,wa,wm]
  File "/ufrc/narayanan/desika.narayanan/pd/SED_gen.py", line 225, in allstars_sed_gen
    p = Pool(processes = cfg.par.n_processes)
  File "/ufrc/narayanan/desika.narayanan/miniconda2/lib/python2.7/multiprocessing/__init__.py", line 232, in Pool
    return Pool(processes, initializer, initargs, maxtasksperchild)
  File "/ufrc/narayanan/desika.narayanan/miniconda2/lib/python2.7/multiprocessing/pool.py", line 159, in __init__
    self._repopulate_pool()
  File "/ufrc/narayanan/desika.narayanan/miniconda2/lib/python2.7/multiprocessing/pool.py", line 223, in _repopulate_pool
    w.start()
  File "/ufrc/narayanan/desika.narayanan/miniconda2/lib/python2.7/multiprocessing/process.py", line 130, in start
    self._popen = Popen(self)
  File "/ufrc/narayanan/desika.narayanan/miniconda2/lib/python2.7/multiprocessing/forking.py", line 121, in __init__
    self.pid = os.fork()
  OSError: [Errno 12] Cannot allocate memory

Try reducing the number of n_processes in the parameters_master.py file
