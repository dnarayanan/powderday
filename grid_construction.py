import random
import numpy as np
import pfh_readsnap
#import parameters as par
import config as cfg
from datetime import datetime
from astropy.table import Table
from astropy.io import ascii
import hyperion_octree_stats as hos

import constants as const

random.seed('octree-demo')
 
import pdb
import os.path
import sys



def yt_octree_generate():
    
    from yt.mods import *
    from yt.geometry.oct_container import OctreeContainer
    from yt.geometry.selection_routines import AlwaysSelector


    fname = cfg.par.hydro_dir+cfg.par.Gadget_snap_name


    #first get the bounding box size
    ptype = 0 #for gas
    print 'in yt_octree_generate: reading in the snapshot with pfh_readsnap'


    

    sdir = cfg.par.hydro_dir
    snum = cfg.par.Gadget_snap_num
    gas_dict = pfh_readsnap.readsnap(sdir,snum,ptype)
    


    metals = gas_dict['z']
    metals = metals[:,0]
    m = gas_dict['m']
    pos = gas_dict['p']
    hsml = gas_dict['h']

    dustmass = m * metals * 1.e10 * 0.4

    
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]

    bbox_lim = max([np.absolute(min(x)),np.absolute(max(x)),
                    np.absolute(min(y)),np.absolute(max(y)),
                    np.absolute(min(z)),np.absolute(max(z))])
    
    #the bounding box needs to be a good bit bigger than the actual
    #min/maxes of the particle locations for constructing the octree
 #   bbox = [[min(x)*2,max(x)*2],
 #           [min(y)*2,max(y)*2],
 #           [min(z)*2,max(z)*2]]

    bbox = [[-2*bbox_lim,2*bbox_lim],
            [-2*bbox_lim,2*bbox_lim],
            [-2*bbox_lim,2*bbox_lim]]
             
    

    unit_base = {'UnitLength_in_cm'         : 3.08568e+21,
                 'UnitMass_in_g'            :   1.989e+43,
                 'UnitVelocity_in_cm_per_s' :      100000}

    print 'NOTE: this assumes the following Gaget parameters which are hard-coded into yt_octree_generate:'
    print unit_base

    #==================================
    #turk's code
    #==================================

    pf = load(fname,unit_base=unit_base,bounding_box=bbox,over_refine_factor=0)
    from yt.data_objects.particle_unions import ParticleUnion
    pu = ParticleUnion("all", list(pf.particle_types_raw))
    
    saved = pf.index.oct_handler.save_octree(always_descend=True)
    
    always = AlwaysSelector(None)
    ir1 = pf.index.oct_handler.ires(always)  #refinement levels
    fc1 = pf.index.oct_handler.fcoords(always)  #coordinates in kpc
    fw1 = pf.index.oct_handler.fwidth(always)  #width of cell in kpc




    #==================================



    refined = saved['octree']
    refined2 = []

 
    for i in range(len(refined)):
        if refined[i] == 1: refined2.append(True)
        if refined[i] == 0: refined2.append(False)

    refined = refined2


    #smooth the data on to the octree
    
    volume = np.zeros(len(refined))
    wTrue = np.where(np.array(refined) == True)[0]
    wFalse = np.where(np.array(refined) == False)[0]
    volume[wFalse] = (fw1 * const.pc * 1.e3)**3.
    

    max_level = pf.index.oct_handler.max_level
    
    


    #get the dust mass via smooth_operator linalg 
    #    dust_mass_grid = smooth_operator.particle_smooth_linalg(x,y,z,hsml,fc1,pos,dustmass,refined)

    
#    import particle_smooth_cython as psnc

    #mass_grid is the smoothed mass_grid; this is just a dummy array
    #that has to get fed into psnc.particle_smooth_new

    mass_grid = np.zeros(len(wFalse))


    if cfg.par.CONSTANT_DUST_GRID == False: #this is the default; if True is set, then we'll 


        '''PSNC.PARTICLE_SMOOTH_NEW STUFF
        print 'Entering psnc.particle_smooth_new'
        temp_dust_mass_grid = psnc.particle_smooth_new(x,y,z,hsml,fc1,dustmass,refined,mass_grid)
        #normalizing for mass conservation
        
        temp_dust_mass_grid /= sum(temp_dust_mass_grid)/sum(m)
    
        #copy over the temp_mass_grid to a grid that is as big as refined
        #(and not just as big as wFalse) for the hyperion calculation
        #(which requires the octree grid include the True's)
        
        dust_mass_grid = np.zeros(len(refined))
        dust_mass_grid[wFalse] = temp_dust_mass_grid
        
        dust_density_grid = dust_mass_grid*const.msun/volume #in gm/cm^-3
        #since volume = 0 where there's a True, dust_density_grid is nan
        #where there's trues, so we have to fix this
        dust_density_grid[wTrue] = 0
        '''


        '''USING YT SMOOTHING'''

        from particle_smooth_yt import yt_smooth
        metallicity_smoothed,mass_smoothed,density_smoothed = yt_smooth(pf)
        dust_smoothed = np.zeros(len(refined))
        

        '''
        dust_smoothed[wFalse] = mass_smoothed * metallicity_smoothed * cfg.par.dusttometals_ratio
        dust_density_grid = dust_smoothed/volume #in gm/cm^-3       
        '''
              
        dust_smoothed[wFalse] = metallicity_smoothed * density_smoothed
        dust_density_grid = dust_smoothed #in gm/cm^-3       
        #since volume = 0 where there's a True, dust_density_grid is nan        
        #where there's trues, so we have to fix this                            
        dust_density_grid[wTrue] = 0
         
        
    else:
        print 'cfg.par.CONSTANT_DUST_GRID=True'
        print 'setting constant dust grid to 4.e-22'
        dust_density_grid = np.zeros(len(refined))+4.e-27
        #since volume = 0 where there's a True, dust_density_grid is nan
        #where there's trues, so we have to fix this
        dust_density_grid[wTrue] = 0






     #file I/O
    print 'Writing Out the Coordinates and Logical Tables'

    xmin = fc1[:,0]-fw1[:,0]
    xmax = fc1[:,0]+fw1[:,0]
    ymin = fc1[:,1]-fw1[:,1]
    ymax = fc1[:,1]+fw1[:,1]
    zmin = fc1[:,2]-fw1[:,2]
    zmax = fc1[:,2]+fw1[:,2]

    coordinates_Table = Table([fc1[:,0]-fw1[:,0],fc1[:,0]+fw1[:,0],fc1[:,1]-fw1[:,1],
                               fc1[:,1]+fw1[:,1],fc1[:,2]-fw1[:,2],fc1[:,2]+fw1[:,2]],
                              names = ['xmin','xmax','ymin','ymax','zmin','zmax'])
    
    ascii.write(coordinates_Table,cfg.par.PD_output_dir+cfg.par.Auto_positions_file)

    logical_Table = Table([refined[:]],names=['logical'])
    ascii.write(logical_Table,cfg.par.PD_output_dir+cfg.par.Auto_TF_file)


    dust_dens_Table = Table([dust_density_grid[:]],names=['dust density'])
    ascii.write(dust_dens_Table,cfg.par.PD_output_dir+cfg.par.Auto_dustdens_file)
        

    return refined,dust_density_grid,xmin,xmax,ymin,ymax,zmin,zmax








def gadget_logical_generate(sdir,snum):
    
   

    #do gas
    ptype = 0
    
    print 'reading in the snapshot with pfh_readsnap'
    gas_dict = pfh_readsnap.readsnap(sdir,snum,ptype)


    #get the dust mass
    #dust mass is the gas mass * gadgetunits * metallicity * 0.4 (= dust to metals ratio)

    metals = gas_dict['z']
    metals = metals[:,0]
    m = gas_dict['m']
    pos = gas_dict['p']
    hsml = gas_dict['h']

    
    x = pos[:,0]
    y = pos[:,1]
    z = pos[:,2]




    n_particles = 5000
    x = x[0:n_particles]
    y = y[0:n_particles]
    z = z[0:n_particles]
    hsml = hsml[0:n_particles]
   



    #construct the octree based on the current gadget grid
    mastercell=[0] #the master cell we'll loop through in the recursive construct_octree

    refined = [True, False, False, False, False, False, False, False, False]
    coordinates = position_calculate(-1.*cfg.par.dx,cfg.par.dx,
                                     -1.*cfg.par.dy,cfg.par.dy,
                                     -1.*cfg.par.dz,cfg.par.dz,
                                     refined)

    print 'constructing the octree: starting at ',str(datetime.now())
    t1 = datetime.now()

    octree_return = construct_octree(x,y,z,hsml,coordinates,mastercell,refined)
    refined = octree_return[0]
    coordinates = octree_return[1]

    t2 = datetime.now()
    print 'total time taken for octree construction: '+str(t2-t1)

    

    #particle smoothing



    print 'smoothing particles on to grid'
    t1 = datetime.now()
    dust_mass_grid = particle_smooth(x,y,z,hsml,coordinates,pos,dustmass,refined)
    
    volume = ((coordinates[:,1]-coordinates[:,0])*const.pc*1.e3)*((coordinates[:,3]-coordinates[:,2])*const.pc*1.e3)*((coordinates[:,5]-coordinates[:,4])*const.pc*1.e3)
    

    base_grid_volume = ((float(cfg.par.dx) * const.pc * 1.e3) * (float(cfg.par.dy) * const.pc * 1.e3) * (float(cfg.par.dz) * const.pc * 1.e3))

    volume = np.insert(volume,0,base_grid_volume) #put the base grid in 
    dust_density_grid = dust_mass_grid*const.msun/volume #in gm/cm^-3

   
    t2 = datetime.now()
    print 'total time for particle smoothing: '+str(t2-t1)
    
    
    #file I/O
    print 'Writing Out the Coordinates and Logical Tables'

    coordinates_Table = Table([coordinates[:,0],coordinates[:,1],coordinates[:,2],
                              coordinates[:,3],coordinates[:,4],coordinates[:,5]],
                              names = ['xmin','xmax','ymin','ymax','zmin','zmax'])
    
    ascii.write(coordinates_Table,cfg.par.PD_output_dir+cfg.par.Auto_positions_file)

    logical_Table = Table([refined[:]],names=['logical'])
    ascii.write(logical_Table,cfg.par.PD_output_dir+cfg.par.Auto_TF_file)


    dust_dens_Table = Table([dust_density_grid[:]],names=['dust density'])
    ascii.write(dust_dens_Table,cfg.par.PD_output_dir+cfg.par.Auto_dustdens_file)








    return refined


def construct_octree(x,y,z,hsml,coordinates_in,mastercell=[0],
                     refined=[True,False,False,False,False,False,False,False,False],
                     refined_levels = 0):


    global coordinates
    coordinates = coordinates_in

  

    for subcell in range(8):
     
        if cfg.par.VERBOSE:
            if refined_levels > 10: print 'refined_levels = ',refined_levels
       

       

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
        coordinates = position_calculate(-1.*cfg.par.dx,cfg.par.dx,
        -1.*cfg.par.dy,cfg.par.dy,
        -1.*cfg.par.dz,cfg.par.dz,
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

            #what is the size of the particles within the current
            #cell?  do they extend past the wall?  if so, add them to
            #the particles_that_are_too_small array


            particles_that_are_too_small = np.where((x[piic]-hsml[piic] > coordinates[mastercell[0]-1][0]) |
                                                    (x[piic]+hsml[piic] < coordinates[mastercell[0]-1][1]) | 
                                                    (y[piic]-hsml[piic] > coordinates[mastercell[0]-1][2]) |
                                                    (y[piic]+hsml[piic] < coordinates[mastercell[0]-1][3]) |
                                                    (z[piic]-hsml[piic] > coordinates[mastercell[0]-1][4]) |
                                                    (z[piic]+hsml[piic] < coordinates[mastercell[0]-1][5]))[0]
            
           


            #if there are any particles that are too small (that fit wholly within the cell) then refine
            if len(particles_that_are_too_small) != 0:
                if cfg.par.VERBOSE: 
                    print 'len (particles_that_are_too_small) = ',len(particles_that_are_too_small)

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
                



               


                construct_octree(x,y,z,hsml,coordinates,mastercell,refined,refined_levels+1)


    return refined,coordinates
















def particle_smooth(x,y,z,hsml,coordinates,pos,m,refined):
 
    #define the grid in terms of refined_nomaster, to get the
    #dimensions right with respect to the coordinates array.  then,
    #later, we'll add a '0' back in the place of the base grid in
    #the final mass_grid


    refined_nomaster = refined[1::]

    mass_grid = np.zeros(len(refined_nomaster))

    wTrue = np.where(np.array(refined_nomaster) == True)[0]

    #loop through the particles and see which cells correspond to
    #being within 3 smoothing lengths of the center of the particle.
    #we ignore any True in the refined list, since this doesn't hold
    #data. but we hvae to include it, because Hyperion expects a value
    #(of 0) in the True spaces.



    for p in range(len(x)): 
        dist = np.sqrt( (coordinates[:,0]-x[p])**2. +
                        (coordinates[:,1]-y[p])**2. +
                        (coordinates[:,2]-z[p])**2. )
        

        a_norm = m[p]
        mass_gauss = a_norm * np.exp( (-1. * (dist**2.))/(2. * (hsml[p]**2.)))
        
        #any True's have to have 0 for data
        mass_gauss[wTrue] = 0 
        

        #we go ahead and blank out any values that are super
        #small. the reason is that if we end up randomly with a
        #gaussian smoothing where ALL the values are really small
        #(<<1e-100), when we try to normalize, the small number
        #calculations will screw up and later mess up the assert
        #statement; better to blank out the grid since its so small
        #and then just assign the mass of the particle of interest to
        #the nearest cell.

        mass_gauss[np.where(mass_gauss < 1e-100)] = 0

        #normalize so that the total mass == the mass of the particle.
        #just in case the distance to the nearest cell is so large
        #that no smoothing is possible (too small values in the
        #gaussian), we assign all of the mass to the nearest cell

        
        if sum(mass_gauss) != 0:
            mass_gauss /= ((sum(mass_gauss))/m[p])  

        else:
            w = np.where(dist == min(dist))
            mass_gauss[w[0]] = m[p]/len(w[0]) #/len(w) just in case a few cells fit the minimum distance criteria
           
        
        
            
        assert(np.isnan(sum(mass_gauss)) == False)
                                            
        
        mass_grid = mass_grid+mass_gauss

      
        #if (sum(mass_gauss)/m[p] > 1.01) | (sum(mass_gauss)/m[p] < 0.99): pdb.set_trace()

        assert( (sum(mass_grid)/sum(m[0:p+1]) < 1.01) & (sum(mass_grid)/sum(m[0:p+1]) > 0.99))
        assert( (sum(mass_gauss)/m[p] < 1.01) & (sum(mass_gauss)/m[p] > 0.99))

    #append a 0 to the beginning to account for the base grid
    mass_grid = np.insert(mass_grid,0,0,0)


        
 
    return mass_grid
    
    







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

