import os

def file_exist(fname):
    if os.path.isfile(fname) == True: pass
    
    else: 
        print "\n\n\n"
        print "======================================================"
        print "File: %s doesn't exist; Powderday Crash!!!"%fname
        print "======================================================"
        print "\n\n\n"

        
    
