import pdb

def sanity_check(refined_string):
    print 'Entering the Octree Sanity Check'
    content = refined_string
    prev = content
 
    for iter in range(16):
 
        content = content.replace('TFFFFFFFF', 'F')
 
        print "Length: {0}".format(len(content))
        if content == prev:
            break
 
            prev = content
 

    print 'Converged: {0} (should be F)'.format(content)

        
