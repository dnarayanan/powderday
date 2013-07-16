#When I restart in Aspen:

#1. Run the grid as follows in ipython
#    a. from grid_construction import *
#    b. positions = {}
#    c. refined = [True (8 falses)]
#    d. refineddict={0:0}
#    dum = recursive_refine(positions,1,0,1,0,1,0,1,refined,0,refineddict)

#this should give us the positions of the 8 subcells (though not the base grid).

#To do:

#the refined index returns back to the master value when we're done
#recursively refining, which is not right and will cause cellvar
#assignment 



#turn the mastercell into [0]
#turn refined into combining with refined_nomaster by refined = [True] + refined_nomaster



import random
import numpy as np
import pfh_readsnap
import parameters as par
from datetime import datetime
from astropy.table import Table
from astropy.io import ascii


random.seed('octree-demo')
 
import pdb


def gadget_logical_generate(sdir,snum):
    
   

    #do gas
    ptype = 0
    
    print 'reading in the snapshot with pfh_readsnap'
    gas_dict = pfh_readsnap.readsnap(sdir,snum,ptype)


    #get the dust mass
    #dust mass is the gas mass * gadgetunits * metallicity * 0.4 (= dust to metals ratio)

    z = gas_dict['z']
    z = z[:,0]
    m = gas_dict['m']
    pos = gas_dict['p']
    hsml = gas_dict['h']

    dustmass = m * z * 1.e10 * 0.4

    
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]

    #construct the octree based on the current gadget grid
    mastercell=[0] #the master cell we'll loop through in the recursive construct_octree

    refined = [True, False, False, False, False, False, False, False, False]
    coordinates = position_calculate(-1.*par.dx,par.dx,
                                     -1.*par.dy,par.dy,
                                     -1.*par.dz,par.dz,
                                     refined)

    print 'constructing the octree: starting at ',str(datetime.now())
    t1 = datetime.now()

    refined = construct_octree(x,y,z,hsml,coordinates,mastercell,refined)

    t2 = datetime.now()
    print 'total time taken for octree construction: '+str(t2-t1)
    
    pdb.set_trace()
    return refined


def construct_octree(x,y,z,hsml,coordinates_in,mastercell=[0],
                     refined=[True,False,False,False,False,False,False,False,False]):


    global coordinates
    coordinates = coordinates_in

    for subcell in range(8):
       
        mastercell[0] += 1

        #criteria for subdividing
        
        #this is a speedup block.  we only want to calculate the
        #coordinates if the grid is refined since this is the rate
        #limiting step (the second coordinate calculation that comes
        #in the conditional) since this is when the 'refined' list is
        #updated.  this said, if the code goes up in one of the
        #'construct_octree' nests, then it will lose information about
        #what the coordinates array is.  the conditional below figures
        #out when this has happened and recalculates the coordinates
        #only these times (versus every time the subcell loop is moved
        #forward), saving over a factor of 3 in computation time.

        '''       
        if len(coordinates) != len(refined)-1:
        print 'recalculating coordinates'
        coordinates = position_calculate(-1.*par.dx,par.dx,
        -1.*par.dy,par.dy,
        -1.*par.dz,par.dz,
        refined)
        '''            
        assert len(coordinates) == len(refined)-1

       
        
      
        

       #we do mastercell[0] -1 since on the first loop through,
       #mastercell = 1 (which is done since the first element of
       #'refined' is the master grid, which we ignore henceforth).
       #The coordinates array is len(refined)-1 rows long, since it
       #doesn't include the master grid

        piic = np.where((x > coordinates[mastercell[0]-1][0]) &
                        (x < coordinates[mastercell[0]-1][1]) & 
                        (y > coordinates[mastercell[0]-1][2]) & 
                        (y < coordinates[mastercell[0]-1][3]) & 
                        (z > coordinates[mastercell[0]-1][4]) & 
                        (z < coordinates[mastercell[0]-1][5]))[0]
        
        #if there are particles inside the cell
        if len(piic) != 0:
            #if none of those particles extend beyond the cell wall
            if np.where((x[piic]-hsml[piic] > coordinates[mastercell[0]-1][0]) &
                        (x[piic]+hsml[piic] < coordinates[mastercell[0]-1][1]) & 
                        (y[piic]-hsml[piic] > coordinates[mastercell[0]-1][2]) &
                        (y[piic]+hsml[piic] < coordinates[mastercell[0]-1][3]) &
                        (z[piic]-hsml[piic] > coordinates[mastercell[0]-1][4]) &
                        (z[piic]+hsml[piic] < coordinates[mastercell[0]-1][5])):

        
                print mastercell
                print 'len(refined) = ',len(refined)
        
               
                refined[mastercell[0]] = True
                for i in range(1,9): refined.insert(mastercell[0]+i,False)
            
               
                #this is a speedup - if we find a True, and are
                #refining, then we only calculate the coordinates for
                #the 8 subcells that correspond to the 8 False's added
                #after the True was inserted (and then insert them
                #into the original coordinates grid). This is *much*
                #quicker than re-calculating the coordinates for the
                #entire refined list, and allows the code to continue
                #at the same pace, rather than getting slowe r and
                #slower (and eventually grinding to a halt) when the
                #refined list gets big enough (>~1e4 elements).

                sub_coordinates = position_calculate(coordinates[mastercell[0]-1][0],
                                                     coordinates[mastercell[0]-1][1],
                                                     coordinates[mastercell[0]-1][2],
                                                     coordinates[mastercell[0]-1][3],
                                                     coordinates[mastercell[0]-1][4],
                                                     coordinates[mastercell[0]-1][5],
                                                     [True,False,False,False,False,False,False,False,False])
                
               
                coordinates = np.insert(coordinates,mastercell[0],sub_coordinates,0)
                



               

                construct_octree(x,y,z,hsml,coordinates,mastercell,refined)


    return refined



'''
def construct_octree2(pos,hsml,mastercell,
                     refined=[True,False,False,False,False,False,False,False,False]):

   

    #debug
    print len(refined)
    print 'mastercell = ',mastercell


    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]

    refined_nomaster = refined[1::] #this is the refined array without the
                               #first "True" representing the master
                               #grid in it.  We need this for looping
                               #through the oct, though need a copy
                               #with the master 'True' in it for the
                               #location calculation.

    #this is the index that we march through in the refined array (and
    #refined_nomaster) as we examine the refined criteria for each
    #cell.  we start with value -1 because we increment early, before
    #we actually start examining individual cells.

   
    if mastercell[0] > 220: 
        print '================'
        print 'mastercell 1 = ',mastercell[0]

    for subcell in range(8):
        
        mastercell[0] += 1

        #======================================================
        #criteria for subdivisions - swap this out if you like
        #======================================================


        #divide = random.random() < 0.12

       

        coordinates = position_calculate(-1.*par.dx,par.dx,
                                         -1.*par.dy,par.dy,
                                         -1.*par.dz,par.dz,
                                         refined)

        if mastercell[0] > 220: 
            print len(refined)
            print len(refined_nomaster)
            print 'mastercell 2 = ',mastercell[0]

 #if the sizes of any particles that are in the cell are
        #smaller than the cell, subdivide

        #find the list of particles that are in the current mastercell
        if mastercell[0] == 289: pdb.set_trace()
        #piic = particle indices inside cell 


        piic = np.where((x > coordinates[mastercell[0]][0]) &
                 (x < coordinates[mastercell[0]][1]) & 
                 (y > coordinates[mastercell[0]][2]) & 
                 (y < coordinates[mastercell[0]][3]) & 
                 (z > coordinates[mastercell[0]][4]) & 
                 (z < coordinates[mastercell[0]][5]))[0]
        

        #if there are particles inside the cell
        if len(piic) != 0:
            #if none of those particles extend beyond the cell wall
            if np.where((x[piic]-hsml[piic] > coordinates[mastercell[0]][0]) &
                        (x[piic]+hsml[piic] < coordinates[mastercell[0]][1]) & 
                        (y[piic]-hsml[piic] > coordinates[mastercell[0]][2]) &
                        (y[piic]+hsml[piic] < coordinates[mastercell[0]][3]) &
                        (z[piic]-hsml[piic] > coordinates[mastercell[0]][4]) &
                        (z[piic]+hsml[piic] < coordinates[mastercell[0]][5])):
                
                
                #======================================================
                #END criteria for subdivisions 
                #====================================================== 

                if mastercell[0] > 200: print 'mastercell 3 = ',mastercell[0]


                #we have passed the criteria for subdivisoins: add 8
                #subcells in the refined list, and change the existing
                #boolean that we're on to a True. When adding 8
                #subcell booleans, we add Falses by default, then
                #switch them to True if necessary later.


                refined_nomaster[mastercell[0]] = True
                
                for i in range(1,9):
                    refined_nomaster.insert(mastercell[0]+i,False)
                    if mastercell[0] > 220: print 'mastercell 4 = ',mastercell[0]

                if mastercell[0] >= 288: pdb.set_trace()
                    #refined[1::] = refined_nomaster
                refined = [True] + refined_nomaster
                if mastercell[0] >= 288: pdb.set_trace()
                
                if mastercell[0] > 220: 
                    print 'bottom pdb'
                    print len(refined)
                    print 'mastercell 5 = ',mastercell[0]

                

                #since we're subdividing, we need to recursively refine
                construct_octree(pos,hsml,mastercell,refined)

       
       
   
       
    
    return refined
        

'''








'''MASTER COPY OF CONSTRUCT_OCTREE
def construct_octree(refined = [True]):
    
    for subcell in range(8):
        
        #criteria for subdivisions
        divide = random.random() < 0.12
         
        #append the boolean to the refined list
        refined.append(divide)
         
         
        #if the cell is subdivided, recursively refine
         
        if divide:
            construct_octree(refined)
            
    return refined
 

'''






def position_calculate(minx,maxx,miny,maxy,minz,maxz, refined = [True]):

    #define the positions dictionary, and the base number of the root cell
    positions = {}

    base_number = 0.
    refinedindex = 0.

    refineddict={0:0}
    


    if refined[0] == True: 
        #refinedindex += 1 #the zeroeth cell is the True which is the base grid
        coordinates = recursive_refine(positions,base_number,minx,maxx,miny,maxy,minz,maxz,
                                       refined,refinedindex,refineddict)

    return coordinates


def recursive_refine(positions,base_number,minx,maxx,miny,maxy,minz,maxz,refined,refinedindex,refineddict):
    #the [minx...maxz] variables are the corners of the grid that is subdividing into an oct of cells



    oct_cell_counter = 1



    #define the 6 wall coordinates of each subcell in the oct
    for ix in range(2):
        for iy in range(2):
            for iz in range(2):

                refineddict[0] += 1
                refinedindex = refineddict[0]
                
                cellvar = 'ref'+str(base_number)+str(oct_cell_counter)
                currentbase = str(base_number)+str(oct_cell_counter)

                #debug print cellvar
                

                cell_minx = ((maxx-minx)/2. * ix) + minx
                cell_maxx = ((maxx-minx)/2. * (ix+1))+minx
                cell_miny = ((maxy-miny)/2. * iy) + miny
                cell_maxy = ((maxy-miny)/2. * (iy+1))+miny
                cell_minz = ((maxz-minz)/2. * iz) + minz
                cell_maxz = ((maxz-minz)/2. * (iz+1)) + minz


                

                positions[cellvar] = [cell_minx,cell_maxx,cell_miny,cell_maxy,cell_minz,cell_maxz]

                




                '''
                the difference between oct_cell_counter and
                refinedindex is that oct_cell_counter resets every
                time we move to a new oct.  refinedindex keeps
                incrementing forward, to march through the 'refined'
                array.
                '''
                

               

                #                refinedindex = refinedindex+1
                
                '''
                print 'current base = ' , currentbase
                print 'oct_cell_counter = ', oct_cell_counter
                print 'refinedindex = ',refinedindex
                '''                

                if refined[refinedindex] == True: 

                    '''
                    coordinates = recursive_refine(positions,currentbase+'.'+str(oct_cell_counter),
                                                   cell_minx,cell_maxx,cell_miny,cell_maxy,
                                                   cell_minz,cell_maxz,refined,refinedindex)
                    '''
                  
                    coordinates = recursive_refine(positions,currentbase+'.',
                                                   cell_minx,cell_maxx,cell_miny,cell_maxy,
                                                   cell_minz,cell_maxz,refined,refinedindex,refineddict)
                    

                oct_cell_counter += 1



    #once we're done recursively refining, 
    #remap the positions to be in order, since the keys get stuck in random order in the dictionary

    sorted_keys = sorted(positions) #these are the keys in numeric order




    #these are the vectors of positions that are ncells long
    x_min_out = []
    y_min_out = []
    z_min_out = []
    x_max_out = []
    y_max_out = []
    z_max_out = []





    for i in range(len(sorted_keys)):
        cellvar = str(sorted_keys[i])
        x_min_out.append(positions[cellvar][0])
        x_max_out.append(positions[cellvar][1])
        y_min_out.append(positions[cellvar][2])
        y_max_out.append(positions[cellvar][3])
        z_min_out.append(positions[cellvar][4])
        z_max_out.append(positions[cellvar][5])


        
    #this is an array of positions that has dimensions [ncells,6], where 6 is for 6 walls

    pos_out = np.zeros((len(positions),6))
    
    for i in range(len(positions)):
        pos_out[i,0] = x_min_out[i]
        pos_out[i,1] = x_max_out[i]
        pos_out[i,2] = y_min_out[i]
        pos_out[i,3] = y_max_out[i]
        pos_out[i,4] = z_min_out[i]
        pos_out[i,5] = z_max_out[i]
        
   
    return(pos_out)

#    return(x_min_out,x_max_out,y_min_out,y_max_out,z_min_out,z_max_out)

        
