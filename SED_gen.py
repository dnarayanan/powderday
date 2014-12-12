import numpy as np
import config as cfg

from astropy.table import Table
from astropy.io import ascii
import astropy.units as u
import astropy.constants as constants
from astropy import cosmology as cosmo

import constants as const
import pdb,ipdb
import sys
from plot_generate import mass_weighted_distribution as mwd
from octree_zoom import octree_zoom_bbox_filter

import fsps 
from datetime import datetime
from datetime import timedelta
from grid_construction import stars_coordinate_boost 

from multiprocessing import Pool
import yt



class Stars:
    def __init__(self,mass,metals,positions,age,sed_bin=[-1,-1,-1],lum=-1,fsps_zmet=20):
        self.mass = mass
        self.metals = metals
        self.positions = positions
        self.age = age
        self.sed_bin = sed_bin
        self.lum = lum
        self.fsps_zmet = fsps_zmet

    def info(self):
        return(self.mass,self.metals,self.positions,self.age,self.sed_bin,self.lum,self.fsps_zmet)






def star_list_gen(boost,xcent,ycent,zcent,dx,dy,dz):
    print '[SED_gen/star_list_gen]: reading in stars particles for SPS calculation'

    fname = cfg.model.hydro_dir+cfg.model.Gadget_snap_name
    
    bbox = [[-2.*cfg.par.bbox_lim,2.*cfg.par.bbox_lim],
            [-2.*cfg.par.bbox_lim,2.*cfg.par.bbox_lim],
            [-2.*cfg.par.bbox_lim,2.*cfg.par.bbox_lim]]
 
    unit_base = {'UnitLength_in_cm'         : cfg.par.unit_length*1.e3*const.pc,
                 'UnitMass_in_g'            : cfg.par.unit_mass*const.msun,
                 'UnitVelocity_in_cm_per_s' : cfg.par.unit_velocity}

    
    if cfg.par.zoom == False:
        pf = yt.load(fname,unit_base=unit_base,bounding_box=bbox,over_refine_factor=cfg.par.oref,n_ref=cfg.par.n_ref)
    else:
        pf = octree_zoom_bbox_filter(fname,unit_base,bbox)
    pf.index

    ad = pf.all_data()


    if  ('PartType4', 'Metallicity_00') in pf.derived_field_list:
        metals = ad[('PartType4', 'Metallicity_00')]
    else:
        metals = ad[("PartType4","Metallicity")].value
    mass = ad[("PartType4","Masses")].value*cfg.par.unit_mass*const.msun
    positions = ad[("PartType4","Coordinates")].value*cfg.par.unit_length*const.pc*1.e3 #cm (as par.unit_length is kpc)


    if cfg.par.COSMOFLAG == False:
        simtime = pf.current_time.value
        simtime *= u.s
        simtime = simtime.to(u.Gyr)
        simtime = simtime.value
        age = simtime-ad[("PartType4","StellarFormationTime")].value * cfg.par.unit_age #gyr

        #make the minimum age 1 million years 
        age[np.where(age < 1.e-3)[0]] = 1.e-3


        print '\n--------------'
        print '[SED_gen/star_list_gen: ] Idealized Galaxy Simulation Assumed: Simulation time is (Gyr): ',simtime
        print '--------------\n'
    else:
        simtime = cosmo.Planck13.age(pf.current_redshift).value #what is the age of the Universe right now?
        scalefactor = ad[("PartType4","StellarFormationTime")].value
        formation_z = (1./scalefactor)-1.

        formation_time = redshift_multithread(formation_z)
        
        age = simtime - formation_time
        #make the minimum age 1 million years 
        age[np.where(age < 1.e-3)[0]] = 1.e-3

        #age = simtime-ad[("PartType4","StellarFormationTime")].value/0.7 * cfg.par.unit_age #gyr
      
        print '\n--------------'
        print '[SED_gen/star_list_gen: ] Cosmological Galaxy Simulation Assumed: Current age of Universe is (Assuming Planck13 Cosmology) is (Gyr): ',simtime
        print '--------------\n'
                

       

    median_metallicity = np.median(metals)
  
    nstars = len(age)
    print 'number of new stars =',nstars
    

    #calculate the fsps interpolated metallicity

    #if the metallicity has many fields, and not just global
    #metallicity then just extract the global metallicity
    if metals.ndim > 1:
        metals = metals[:,0]


    print '[SED_gen/star_list_gen:] Manually increasing the newstar metallicities by: ',cfg.par.Z_init
    metals += cfg.par.Z_init


    zmet = fsps_metallicity_interpolate(metals)
    #mwd(zmet,mass,'zmet_distribution.png')

    #print '[SED_gen/star_list_gen: ] fsps zmet codes:',zmet

    #create the stars_list full of Stars objects
    stars_list = []

    for i in range(nstars):
        stars_list.append(Stars(mass[i],metals[i],positions[i],age[i],fsps_zmet=zmet[i]))


    #boost stellar positions to grid center
    print 'boosting new stars to coordinate center'
    stars_list = stars_coordinate_boost(stars_list,boost)
    


   

    #if zoom is set, then pop off the stars that are not in the actual hyperion grid
    if cfg.par.zoom == True:
        orig_stars_list_len = len(stars_list)
        #we have to go through the loop in reverse so that we don't throw off the subsquent indices
        print 'Popping Stars out of Grid (potentially)...'




        for i in reversed(range(nstars)):
            
            if np.logical_or(stars_list[i].positions[0] <= xcent-dx,stars_list[i].positions[0] >= xcent+dx) == True: 
                stars_list.pop(i)
                continue

            if np.logical_or(stars_list[i].positions[1] <= ycent-dy,stars_list[i].positions[1] >= ycent+dy) == True: 
                stars_list.pop(i)
                continue

            if np.logical_or(stars_list[i].positions[2] <= zcent-dz,stars_list[i].positions[2] >= zcent+dz) == True: 
                stars_list.pop(i)
                continue

        len_difference = orig_stars_list_len-len(stars_list)
        print '\n [SED_gen]: Warning: Out of %d new stars\n'%nstars
        print '\n [SED_gen]: Warning: threw out %d stars or being outside the zoomed grid \n'%len_difference
        


 
    #ASSIGN DISK AND BULGE STARS - note, if these don't exist, it will
    #just make empty lists

   

    
    bulgestars_list = []
    diskstars_list = []


    
    #in principle, we should just be able to do the following blocks
    #if the particle types exist. the issue is that different groups
    #use PartType2 and 3 as 'filler' particle types, so they may exist
    #even if they don't correspond to disk/bulge stars.

    if cfg.par.COSMOFLAG == False:

        #Disk Stars
        if ("PartType2","Coordinates") in pf.derived_field_list:
            disk_positions = ad[("PartType2","Coordinates")].value*cfg.par.unit_length*const.pc*1.e3 #cm (as par.unit_length is kpc)
            disk_masses =  ad[("PartType2","Masses")].value*cfg.par.unit_mass*const.msun
            nstars_disk = len(disk_masses)
     
            #create the disk_list full of DiskStars objects
            for i in range(nstars_disk):
                diskstars_list.append(Stars(disk_masses[i],0.02,disk_positions[i],cfg.par.disk_stars_age))

            print 'boosting disk stars to coordinate center'    
            diskstars_list = stars_coordinate_boost(diskstars_list,boost)

            #if zoom is set, then pop off the stars that are not in the actual hyperion grid
            if cfg.par.zoom == True:
                orig_disk_stars_list_len = nstars_disk
            

                #we have to go through the loop in reverse so that we don't throw off the subsquent indices
                for i in reversed(range(nstars_disk)):
                    if np.logical_or(diskstars_list[i].positions[0] <= xcent-dx,diskstars_list[i].positions[0] >= xcent+dx) == True: 
                        diskstars_list.pop(i)
                        continue
                    if np.logical_or(diskstars_list[i].positions[1] <= ycent-dy,diskstars_list[i].positions[1] >= ycent+dy) == True: 
                        diskstars_list.pop(i)
                        continue
                    if np.logical_or(diskstars_list[i].positions[2] <= zcent-dz,diskstars_list[i].positions[2] >= zcent+dz) == True: 
                        diskstars_list.pop(i)
                        continue

                disk_len_difference = orig_disk_stars_list_len-len(diskstars_list)
                print '\n [SED_gen]: Warning: Out of %d disk stars\n'%nstars_disk
                print '\n [SED_gen]: Warning: threw out %d disk stars for being outside the zoomed grid \n'%disk_len_difference


        #Bulge Stars
        if ("PartType3","Coordinates") in pf.derived_field_list:
            bulge_positions = ad[("PartType3","Coordinates")].value*cfg.par.unit_length*const.pc*1.e3 #cm (as par.unit_length is kpc)
            bulge_masses =  ad[("PartType3","Masses")].value*cfg.par.unit_mass*const.msun
            nstars_bulge = len(bulge_masses)

            #create the bulge_list full of BulgeStars objects
            
            for i in range(nstars_bulge):
                bulgestars_list.append(Stars(bulge_masses[i],0.02,bulge_positions[i],cfg.par.bulge_stars_age))
                

            print 'boosting bulge stars to coordinate center'    
            bulgestars_list = stars_coordinate_boost(bulgestars_list,boost)

            #if zoom is set, then pop off the stars that are not in the actual hyperion grid
            if cfg.par.zoom == True:
                orig_bulge_stars_list_len = nstars_bulge

                for i in reversed(range(nstars_bulge)):
                    if np.logical_or(bulgestars_list[i].positions[0] <= xcent-dx,bulgestars_list[i].positions[0] >= xcent+dx) == True: 
                        bulgestars_list.pop(i)
                        continue
                    if np.logical_or(bulgestars_list[i].positions[1] <= ycent-dy,bulgestars_list[i].positions[1] >= ycent+dy) == True: 
                        bulgestars_list.pop(i)
                        continue
                    if np.logical_or(bulgestars_list[i].positions[2] <= zcent-dz,bulgestars_list[i].positions[2] >= zcent+dz) == True: 
                        bulgestars_list.pop(i)
                        continue


                bulge_len_difference = orig_bulge_stars_list_len-len(bulgestars_list)
                print '\n [SED_gen]: Warning: Out of %d bulge stars\n'%nstars_bulge
                print '\n [SED_gen]: Warning: threw out %d bulge stars for being outside the zoomed grid \n'%bulge_len_difference

       
    return stars_list,diskstars_list,bulgestars_list



def allstars_sed_gen(stars_list,diskstars_list,bulgestars_list):

    
    #NOTE this part is just for the gadget simulations - this will
    #eventually become obviated as it gets passed into a function to
    #populate the stars_list with objects as we start to feed in new
    #types of simulation results.

    nstars = len(stars_list)
    nstars_disk = len(diskstars_list)
    nstars_bulge = len(bulgestars_list)
    

    #get just the wavelength array
    sp = fsps.StellarPopulation(tage=stars_list[0].age,imf_type=cfg.par.imf_type,pagb = cfg.par.pagb,sfh=0,zmet=stars_list[0].fsps_zmet)
    
    spec = sp.get_spectrum(tage=stars_list[0].age,zmet=stars_list[0].fsps_zmet)
    nu = 1.e8*const.c/spec[0]
    nlam = len(nu)



    #initialize the process pool and build the chunks
    p = Pool(processes = cfg.par.n_processes)
    nchunks = cfg.par.n_processes


    chunk_start_indices = []
    chunk_start_indices.append(0) #the start index is obviously 0

    delta_chunk_indices = int(nstars / nchunks)
    print 'delta_chunk_indices = ',delta_chunk_indices
    
    for n in range(1,nchunks):
        chunk_start_indices.append(chunk_start_indices[n-1]+delta_chunk_indices)

    '''
    chunk_start_indices = list(np.fix(np.arange(0,nstars,np.fix(nstars/nchunks))))
    #because this can result in too many chunks sometimes given the number of processors:
    chunk_start_indices = chunk_start_indices[0:nchunks]
    '''
    print 'Entering Pool.map multiprocessing for Stellar SED generation'
    list_of_chunks = []
    for n in range(nchunks):
        stars_list_chunk = stars_list[chunk_start_indices[n]:chunk_start_indices[n]+delta_chunk_indices]
        #if we're on the last chunk, we might not have the full list included, so need to make sure that we have that here
        if n == nchunks-1: 
            stars_list_chunk = stars_list[chunk_start_indices[n]::]

        list_of_chunks.append(stars_list_chunk)
    
    

        
    t1=datetime.now()
    chunk_sol = p.map(newstars_gen, [arg for arg in list_of_chunks])
    
    t2=datetime.now()
    print 'Execution time for SED generation in Pool.map multiprocessing = '+str(t2-t1)

    
    stellar_fnu = np.zeros([nstars,nlam])
    star_counter=0
    for i in range(nchunks):
        fnu_list = chunk_sol[i] #this is a list of the stellar_fnu's returned by that chunk
        for j in range(len(fnu_list)):
            stellar_fnu[star_counter,:] = fnu_list[j,:]
            star_counter+=1




    p.close()
    p.terminate()
    p.join()


    stellar_nu = nu



    if cfg.par.COSMOFLAG == False:

        #calculate the SED for disk stars; note, this gets calculated
        #whether or not disk stars actually exist.  if they don't exist,
        #bogus values for the disk age and metallicity are assigned based
        #on whatever par.disk_stars_age and metallicity are.  it's no big
        #deal since these SEDs don't end up getting added to the model in
        #source_creation as long as COSMOFLAG == True.  
        
        #note, even if there are no disk/bulge stars, these are still
        #created since they're completely based on input parameters in
        #parameters_master.  they just won't get used at a later point
        #as there will be no disk/bulge star positions to add them to.

        #dust_tesc is an absolute value (not relative to min star age) as the ages of these stars are input by the user
        sp = fsps.StellarPopulation(tage = cfg.par.disk_stars_age,imf_type=cfg.par.imf_type,pagb = cfg.par.pagb,sfh=0,zmet=cfg.par.disk_stars_metals)
        spec = sp.get_spectrum(tage=cfg.par.disk_stars_age,zmet=20)
        disk_fnu = spec[1]
        
        #calculate the SED for bulge stars
        sp = fsps.StellarPopulation(tage = cfg.par.bulge_stars_age,imf_type=cfg.par.imf_type,pagb = cfg.par.pagb,sfh=0,zmet=cfg.par.bulge_stars_metals)
        spec = sp.get_spectrum(tage=cfg.par.bulge_stars_age,zmet=20)
        bulge_fnu = spec[1]
    

    else: #we have a cosmological simulation

        disk_fnu = []
        bulge_fnu = []
        

    total_lum_in_sed_gen = 0.
    for i in range(stellar_fnu.shape[0]):
        total_lum_in_sed_gen += np.absolute(np.trapz(stellar_fnu[i,:],x=nu))

    print '[SED_gen: ] total_lum_in_sed_gen = ',total_lum_in_sed_gen

    #return positions,disk_positions,bulge_positions,mass,stellar_nu,stellar_fnu,disk_masses,disk_fnu,bulge_masses,bulge_fnu
    return stellar_nu,stellar_fnu,disk_fnu,bulge_fnu


def redshift_multithread(formation_z):

        formation_z_list = formation_z
        #initialize the process pool and build the chunks
        p = Pool(processes = cfg.par.n_processes)
        nchunks = cfg.par.n_processes
        chunk_start_indices = []
        chunk_start_indices.append(0) #the start index is obviously 0
        delta_chunk_indices = int(len(formation_z_list) / nchunks)
        print 'delta_chunk_indices = ',delta_chunk_indices
        for n in range(1,nchunks):
            chunk_start_indices.append(chunk_start_indices[n-1]+delta_chunk_indices)
        list_of_chunks = []
        for n in range(nchunks):
            formation_z_chunk = formation_z_list[chunk_start_indices[n]:chunk_start_indices[n]+delta_chunk_indices]
            if n == nchunks-1: 
                formation_z_chunk = formation_z_list[chunk_start_indices[n]::]
            list_of_chunks.append(formation_z_chunk)
        print 'Entering Pool.map multiprocessing for Stellar Age calculations'
        t1=datetime.now()
        chunk_sol = p.map(redshift_gen, [arg for arg in list_of_chunks])
        

        formation_time = []
        for i in range(len(chunk_sol)):
            sub_chunk_sol = chunk_sol[i].value
            for j in range(len(sub_chunk_sol)):
                formation_time.append(sub_chunk_sol[j])

        t2=datetime.now()
        print 'Execution time for Stellar Age calculations in Pool.map multiprocessing = '+str(t2-t1)
        
      
        
        return np.array(formation_time)
        

def redshift_gen(formation_z):
    age = cosmo.Planck13.age(formation_z)
    return age

def newstars_gen(stars_list):
    #the newstars (particle type 4; so, for cosmological runs, this is all
    #stars) are calculated in a separate function with just one argument so that it is can be fed 
    #into pool.map for multithreading.
    

    #first figure out how many wavelengths there are
    sp = fsps.StellarPopulation(tage=stars_list[0].age,imf_type=cfg.par.imf_type,pagb = cfg.par.pagb,sfh=0,zmet=stars_list[0].fsps_zmet)
    spec = sp.get_spectrum(tage=stars_list[0].age,zmet=stars_list[0].fsps_zmet)
    nu = 1.e8*const.c/spec[0]
    

    nlam = len(nu)

    stellar_nu = np.zeros([nlam])
    stellar_fnu = np.zeros([len(stars_list),nlam])
    
  
    minage = 13 #Gyr
    for i in range(len(stars_list)): 
        if stars_list[i].age < minage:
            minage = stars_list[i].age

    tesc_age = np.log10((minage+cfg.par.birth_cloud_clearing_age)*1.e9)


    #calculate the SEDs for new stars
    for i in range(len(stars_list)):
        
        if cfg.par.CF_on == True:
            sp = fsps.StellarPopulation(tage=stars_list[i].age,imf_type=cfg.par.imf_type,pagb = cfg.par.pagb,sfh=0,zmet=stars_list[i].fsps_zmet,dust_type=0,dust1=1,dust2=0,dust_tesc=tesc_age)
        else:
            sp = fsps.StellarPopulation(tage=stars_list[i].age,imf_type=cfg.par.imf_type,pagb = cfg.par.pagb,sfh=0,zmet=stars_list[i].fsps_zmet)
         

        #sp = fsps.StellarPopulation(tage=stars_list[i].age,imf_type=2,sfh=0,zmet=stars_list[i].fsps_zmet)
        spec = sp.get_spectrum(tage=stars_list[i].age,zmet=stars_list[i].fsps_zmet)

        stellar_nu[:] = 1.e8*const.c/spec[0]
        stellar_fnu[i,:] = spec[1]

    return stellar_fnu






def fsps_metallicity_interpolate(metals):

    #takes a list of metallicities for star particles, and returns a
    #list of interpolated metallicities
    
    fsps_metals = np.loadtxt(cfg.par.metallicity_legend)
    nstars = len(metals)
    
    zmet = []

    for i in range(nstars):
        zmet.append(find_nearest_zmet(fsps_metals,metals[i]))
    
    return zmet
    
def find_nearest_zmet(array,value):
    #this is modified from the normal find_nearest in that it forces
    #the output to be 1 index higher than the true value since the
    #minimum zmet value fsps will take is 1 (not 0)

    idx = (np.abs(array-value)).argmin()
    
    return idx+1


