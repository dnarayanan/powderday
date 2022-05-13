from __future__ import print_function
import os
import sys
import powderday.config as cfg

def check_parameter_compatibility():
    if cfg.par.otf_extinction==True:
        try:
            assert(cfg.par.dust_grid_type == "manual")
        except AssertionError:
            raise AssertionError("otf_extinction is set in parameters_master: this means the dust_grid_type must be manual.  it is currently set as something else.")


def file_exist(fname):
    if os.path.isfile(fname) == True: pass
    
    else: 
        print ("\n\n\n")
        print ("======================================================")
        print ("File: %s doesn't exist; Powderday Crash!!!"%fname)
        print ("======================================================")
        print ("\n\n\n")
        sys.exit()
        
    
