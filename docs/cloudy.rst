Nebular Emission 
**********
This section explains how we include nebular line emission within powderday. We use `CLOUDY <https://www.nublado.org/>`_ 
photoionization code for all our calculations. There are two different methods that are available for calculating nebular emission. 
These can be selected by setting *use_cloudy_tables* flag in parameters_master.py .

Look up tables
--------------
If *use_cloudy_tables* is set to True, CLOUDY look up tables (dev. by Nell Byler) that come pre-packaged with `FSPS <https://github.com/cconroy20/fsps>`_  
are used to calculate nebular emission.


Running CLOUDY (Under active development)
--------------
If *use_cloudy_tables* is set to False, we calculate nebular emission by running CLOUDY for each young star particle. 
This makes the code slower but allows much more freedom over the model parameters. To use this option you 
need to have CLOUDY (version 17.00) installed on your machine . In order to install CLOUDY go this 
`link <https://nublado.org/wiki/StepByStep/>`_ and follow the step by step guide. 

Please make sure that you install *version 17.00*, using any other verison may lead to errors or inconsistent results. 
Also, make sure that you have CLOUDY_EXE and CLOUDY_DATA_PATH environmant variables set. If for any reason the CLOUDY run crashes 
the code reverts back to using look up tables for that star particle.

There are a whole host of parameters that you can change as per your need. A complete list of all the tunable parameters 
is described in parameters_description section. 
