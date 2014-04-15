import pdb

def sanity_check(refined_string,max_level):
    print 'Entering the Octree Sanity Check'
    content = refined_string
    prev = content
 
    print 'inside sanity_check: max_level = ',max_level

    for iter in range(max_level):
 
        content = content.replace('TFFFFFFFF', 'F')
 
        print "Length: {0}".format(len(content))
        if content == prev:
            break
 
            prev = content
 

    print 'Converged: {0} (should be F)'.format(content)

        
