from __future__ import print_function
import os
import sys

def file_exist(fname):
    if os.path.isfile(fname) == True: pass
    
    else: 
        print ("\n\n\n")
        print ("======================================================")
        print ("File: %s doesn't exist; Powderday Crash!!!"%fname)
        print ("======================================================")
        print ("\n\n\n")
        sys.exit()
        
    
