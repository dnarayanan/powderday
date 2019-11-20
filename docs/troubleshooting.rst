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
  
2. Python-fsps compilation can't find 'f95'
--------------
Error (see `here <https://github.com/dfm/python-fsps/issues/56>`_)::

     Could not locate executable f95
     ...
     error: Command "f95 -fPIC -fPIC -O3 ..." failed with exit status 127

Try setting an alias to gfortran in your :code:`.bashrc`::

    alias f95='gfortran'
 
If this doesn't work, create a soft link to the gfortran executable in your path, e.g.::
    
    ln -s /apps/compilers/gcc/6.3.0/bin/gfortran ~/bin/f95


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

5. Assertion Error in SED generation related to zmet being out of range
--------------
Error::

  assigning stars to SED bins
  Running SPS for Binned SEDs
  calculating the SEDs for  105  bins
  Traceback (most recent call last):
  File "pd_front_end.py", line 134, in <module>
  m=add_binned_seds(df_nu,stars_list,diskstars_list,bulgestars_list,m)
  File "source_creation.py", line 300, in add_binned_seds
  binned_stellar_nu,binned_stellar_fnu_has_stellar_mass,disk_fnu,bulge_fnu = sg.allstars_sed_gen(sed_bins_list_has_stellar_mass,diskstars_list,bulgestars_list)
  File "SED_gen.py", line 216, in allstars_sed_gen
    add_neb_emission = cfg.par.add_neb_emission, add_agb_dust_model=cfg.par.add_agb_dust_model)
  File "build/bdist.linux-x86_64/egg/fsps/fsps.py", line 468, in __init__
    self.params[k] = kwargs.pop(k, v)
  File "build/bdist.linux-x86_64/egg/fsps/fsps.py", line 1093, in __setitem__
    self.check_params()
  File "build/bdist.linux-x86_64/egg/fsps/fsps.py", line 1067, in check_params
    "zmet={0} out of range [1, {1}]".format(self._params["zmet"], NZ)
    AssertionError: zmet=20 out of range [1, 12]

Recompile fsps with the libraries set to Padova (and not MIST) in sps_vars.f90. Also recompile python-fsps


6. Pool.map errors in powderday
--------------

Freezing during ``Pool.map`` and `'metallicity outside of range'` errors::

    Entering Pool.map multiprocessing for Stellar SED generation
    SSP_GEN ERROR: metallicity outside of range          14
    SSP_GEN ERROR: metallicity outside of range          15
    ...

Some installations have encountered this issue, but its cause has not yet been 
determined. One potential fix could be using 
`Miniconda <https://repo.continuum.io/miniconda/>`_ instead of 
`Anaconda <https://www.anaconda.com/distribution/>`_ Python, although this has 
not been confirmed to be the source of the problem. If something went wrong at 
any point in the installation process, starting from scratch and doing a fresh 
installation may also fix the issue.


6. Can't find "builtins"
--------------

During run, we get an error that looks like::

    
  (pd4env) [desika.narayanan@login2 pd_git]$ pd_front_end.py  tests/SKIRT/gizmo_mw_zoom/ parameters_master_401 parameters_model_401
  Traceback (most recent call last):
  File "/home/desika.narayanan/miniconda3/envs/pd4env/bin/pd_front_end.py", line 4, in <module>
    __import__('pkg_resources').run_script('powderday==0.1.0', 'pd_front_end.py')
  File "/home/desika.narayanan/miniconda3/envs/pd4env/lib/python2.7/site-packages/pkg_resources/__init__.py", line 666, in run_script
    self.require(requires)[0].run_script(script_name, ns)
  File "/home/desika.narayanan/miniconda3/envs/pd4env/lib/python2.7/site-packages/pkg_resources/__init__.py", line 1469, in run_script
    exec(script_code, namespace, namespace)
  File "/home/desika.narayanan/miniconda3/envs/pd4env/lib/python2.7/site-packages/powderday-0.1.0-py2.7.egg/EGG-INFO/scripts/pd_front_end.py", line 7, in <module>

  File "build/bdist.linux-x86_64/egg/powderday/__init__.py", line 4, in <module>
  File "build/bdist.linux-x86_64/egg/powderday/SED_gen.py", line 19, in <module>
  File "build/bdist.linux-x86_64/egg/powderday/nebular_emission/cloudy_model.py", line 2, in <module>
  File "build/bdist.linux-x86_64/egg/powderd

    ...

Try the following::
  pip install future

