from __future__ import print_function
#source_creation.py
#purpose: to add new star sources, disks and bulge stars
import powderday.config as cfg
import numpy as np
import powderday.SED_gen as sg
from datetime import datetime
import astropy.units as units
import astropy.constants as constants
from powderday.helpers import find_nearest
from powderday.analytics import dump_AGN_SEDs


class Sed_Bins:
    def __init__(self,mass,metals,age,fsps_zmet):
        self.mass = mass
        self.metals = metals
        self.age = age
        self.fsps_zmet=fsps_zmet


def add_newstars(df_nu,stellar_nu,stellar_fnu,disk_fnu,bulge_fnu,stars_list,diskstars_list,bulgestars_list,cosmoflag,m):
    
        
    nstars = len(stars_list)
    print ('adding new stars to the grid')
    
    totallum_newstars = 0.

    for i in range(nstars):
        
        nu = stellar_nu[:]
        fnu = stellar_fnu[i,:]
    
        
        nu,fnu = wavelength_compress(nu,fnu,df_nu)
        #reverse the arrays for hyperion
        nu = nu[::-1]
        fnu = fnu[::-1]

        
        
        lum = np.absolute(np.trapz(fnu,x=nu))*stars_list[i].mass/constants.M_sun.cgs.value 
        lum *= constants.L_sun.cgs.value

        #add new stars
        totallum_newstars += lum
        #m.add_spherical_source(luminosity = lum,radius = 10.*const.rsun,spectrum = (nu,fnu),
        #position = stars_list[i].positions)
        m.add_point_source(luminosity = lum,spectrum=(nu,fnu),position = stars_list[i].positions)
        
      
                           


    print ('[source_creation/add_newstars:] totallum_newstars = ',totallum_newstars)
        
    if cosmoflag == False: add_bulge_disk_stars(df_nu,stellar_nu,stellar_fnu,disk_fnu,bulge_fnu,stars_list,diskstars_list,bulgestars_list,m)
            
    m.set_sample_sources_evenly(True)
    
    return m

def add_bulge_disk_stars(df_nu,stellar_nu,stellar_fnu,disk_fnu,bulge_fnu,stars_list,diskstars_list,bulgestars_list,m):
    print ('Non-Cosmological Simulation: Adding Disk and Bulge Stars:')
    

    nu = stellar_nu[:]
    nu,bulge_fnu = wavelength_compress(nu,bulge_fnu,df_nu)
    nu,disk_fnu = wavelength_compress(nu,disk_fnu,df_nu)


    #reverse the arrays for hyperion
    nu = nu[::-1]
    


    nstars_disk = len(diskstars_list)
    nstars_bulge = len(bulgestars_list)
    

    if nstars_disk >0: 

        fnu = disk_fnu[:]

        #reverse the arrays for hyperion
        fnu = fnu[::-1]


        print ('adding disk stars to the grid: adding as a point source collection')
        disksource = m.add_point_source_collection()
            
        
        disk_lum = np.absolute(np.trapz(fnu,x=nu))*diskstars_list[0].mass/constants.M_sun.cgs.value
        #since stellar masses are in cgs, and we need them to be in msun - we
        #multiply by mass to get the *total* luminosity of the stellar
        #cluster since int(nu,fnu) is just the luminosity of a 1 Msun single star
        disk_lum *= constants.L_sun.cgs.value
        disksource.luminosity = np.repeat(disk_lum,nstars_disk)
        
        disk_pos = np.zeros([len(diskstars_list),3])
 
        for i in range(len(diskstars_list)): disk_pos[i,:] = diskstars_list[i].positions
        disksource.position=disk_pos
        
        disksource.spectrum = (nu,fnu)
    
        print ('[source_creation/add_bulge_disk_stars:] totallum_disktars = ',disksource.luminosity[0])
    


    if nstars_bulge > 0:

        fnu = bulge_fnu[:]
        #reverse the arrays for hyperion
        fnu = fnu[::-1]



        print ('adding bulge stars to the grid: adding as a point source collection')
        bulgesource = m.add_point_source_collection()
        bulge_lum = np.absolute(np.trapz(fnu,x=nu))*bulgestars_list[0].mass/constants.M_sun.cgs.value
        bulge_lum *= constants.L_sun.cgs.value
        bulgesource.luminosity = np.repeat(bulge_lum,nstars_bulge)
        
        bulge_pos = np.zeros([len(bulgestars_list),3])
        for i in range(len(bulgestars_list)): bulge_pos[i,:] = bulgestars_list[i].positions
        bulgesource.position=bulge_pos
        bulgesource.spectrum = (nu,fnu)
        
        print ('[source_creation/add_bulge_disk_stars:] totallum_bulgetars = ',bulgesource.luminosity[0])

    

def add_binned_seds(df_nu,stars_list,diskstars_list,bulgestars_list,cosmoflag,m,sp):
    

  

    #calculate max and min ages
    minimum_age = 15 #Gyr - obviously too high of a number
    maximum_age = 0 #Gyr


    #calculate the minimum and maximum luminosity
    minimum_mass = 1e15*constants.M_sun.cgs.value #msun - some absurdly large value for a single stellar cluster
    maximum_mass = 0 #msun

    #calculate the minimum and maximum stellar metallicity
    minimum_metallicity = 1.e5 #some absurdly large metallicity
    maximum_metallicity = 0

    nstars = len(stars_list)
    for i in range(nstars):
        #if stars_list[i].metals[0] < minimum_metallicity: minimum_metallicity = stars_list[i].metals[0]
        #if stars_list[i].metals[0] > maximum_metallicity: maximum_metallicity = stars_list[i].metals[0]
        
        if stars_list[i].mass < minimum_mass: minimum_mass = stars_list[i].mass
        if stars_list[i].mass > maximum_mass: maximum_mass = stars_list[i].mass

        if stars_list[i].age < minimum_age: minimum_age = stars_list[i].age
        if stars_list[i].age > maximum_age: maximum_age = stars_list[i].age

      


    #define the metallicity bins: we do this by saying that they are the number of metallicity bins in FSPS

    fsps_metals = np.loadtxt(cfg.par.metallicity_legend)
    N_METAL_BINS = len(fsps_metals)

    #note the bins are NOT metallicity, but rather the zmet keys in
    #fsps (i.e. the zmet column in Table 1 of the fsps manual)
    metal_bins = np.arange(N_METAL_BINS)+1



    delta_age = (maximum_age-minimum_age)/cfg.par.N_STELLAR_AGE_BINS

    
    
    #define the age bins in log space so that we maximise resolution around young stars
    age_bins = 10.**(np.linspace(np.log10(minimum_age),np.log10(maximum_age),cfg.par.N_STELLAR_AGE_BINS))
    #tack on the maximum age bin
    age_bins = np.append(age_bins,age_bins[-1]+delta_age)

   
    #define the mass bins (log)
    #note - for some codes, all star particles have the same mass.  in this case, we have to have a trap:
    if minimum_mass == maximum_mass: 
        mass_bins = np.zeros(cfg.par.N_MASS_BINS+1)+minimum_mass
    else:
        delta_mass = (np.log10(maximum_mass)-np.log10(minimum_mass))/cfg.par.N_MASS_BINS
        mass_bins = np.arange(np.log10(minimum_mass),
                              np.log10(maximum_mass),
                              delta_mass)
        mass_bins = np.append(mass_bins,mass_bins[-1]+delta_mass)
        mass_bins = 10.**mass_bins
        

    print ('mass_bins = ',mass_bins)
    print ('metal_bins = ',metal_bins)
    print ('age_bins = ',age_bins)


    #has_stellar_mass is a 3D boolean array that's [wz,wa,wm] big and
    #says whether or not that bin is being used downstream for
    #creating a point source collection (i.e. that it actually has at
    #least one star cluster that falls into it)
    has_stellar_mass = np.zeros([N_METAL_BINS,cfg.par.N_STELLAR_AGE_BINS+1,cfg.par.N_MASS_BINS+1],dtype=bool)


    stars_in_bin = {} #this will be a dictionary that holds the list
    #of star particles that go in every [wz,wa,wm]
    #group.  The keys will be tuples that hold a
    #(wz,wa,wm) set that we will then use later to
    #speed up adding sources.
    
    for i in range(nstars):
        
        wz = find_nearest(metal_bins,stars_list[i].fsps_zmet)
        wa = find_nearest(age_bins,stars_list[i].age)
        wm = find_nearest(mass_bins,stars_list[i].mass)
        
        stars_list[i].sed_bin = [wz,wa,wm]

        has_stellar_mass[wz,wa,wm] = True

        if (wz,wa,wm) in stars_in_bin:
            stars_in_bin[(wz,wa,wm)].append(i)
        else:
            stars_in_bin[(wz,wa,wm)] = [i]

   



    


   
    
    print ('assigning stars to SED bins')
    sed_bins_list=[]
    sed_bins_list_has_stellar_mass = []


   
    #we loop through age bins +1 because the max values were tacked
    #onto those bin lists. but for metal bins, this isn't the case, so
    #we don't loop the extra +1
    for wz in range(N_METAL_BINS):
        for wa in range(cfg.par.N_STELLAR_AGE_BINS+1):
            for wm in range(cfg.par.N_MASS_BINS+1):
                sed_bins_list.append(Sed_Bins(mass_bins[wm],fsps_metals[wz],age_bins[wa],metal_bins[wz]))
                if has_stellar_mass[wz,wa,wm] == True:
                    sed_bins_list_has_stellar_mass.append(Sed_Bins(mass_bins[wm],fsps_metals[wz],age_bins[wa],metal_bins[wz]))

   

    #sed_bins_list is a list of Sed_Bins objects that have the
    #information about what mass bin, metal bin and age bin they
    #correspond to.  It is unnecessary, and heavy computational work
    #to re-create the SED for each of these bins - rather, we can just
    #calculate the SED for the bins that have any actual stellar mass.
            
    print ('Running SPS for Binned SEDs')
    print ('calculating the SEDs for ',len(sed_bins_list_has_stellar_mass),' bins')
    binned_stellar_nu,binned_stellar_fnu_has_stellar_mass,disk_fnu,bulge_fnu = sg.allstars_sed_gen(sed_bins_list_has_stellar_mass,cosmoflag,sp)



   


    #since the binned_stellar_fnu_has_stellar_mass is now
    #[len(sed_bins_list_has_stellar_mass),nlam)] big, we need to
    #transform it back to the a larger array.  this is an ugly loop
    #that could probably be prettier...but whatever.  this saves >an
    #order of magnitude in time in SED gen.  
    nlam = binned_stellar_nu.shape[0]
    binned_stellar_fnu = np.zeros([len(sed_bins_list),nlam])

    counter = 0
    counter_has_stellar_mass = 0
    for wz in range(N_METAL_BINS):
        for wa in range(cfg.par.N_STELLAR_AGE_BINS+1):
            for wm in range(cfg.par.N_MASS_BINS+1):
                if has_stellar_mass[wz,wa,wm] == True:
                    binned_stellar_fnu[counter,:] = binned_stellar_fnu_has_stellar_mass[counter_has_stellar_mass,:]
                    counter_has_stellar_mass += 1 
                counter+=1





    
    '''
    #DEBUG trap for nans and infs
    if np.isinf(np.sum(binned_stellar_nu)):  pdb.set_trace()
    if np.isinf(np.sum(binned_stellar_fnu)): pdb.set_trace()
    if np.isnan(np.sum(binned_stellar_nu)): pdb.set_trace()
    if np.isnan(np.sum(binned_stellar_fnu)): pdb.set_trace()
    '''

    #now binned_stellar_nu and binned_stellar_fnu are the SEDs for the bins in order of wz, wa, wm 
    
    #create the point source collections: we loop through the bins and
    #see what star particles correspond to these.  if any do, then we
    #add them to a list, and create a point source collection out of
    #these


    print ('adding point source collections')
    t1=datetime.now()


    totallum = 0 
    totalmass = 0 

    counter=0
    for wz in range(N_METAL_BINS):
        for wa in range(cfg.par.N_STELLAR_AGE_BINS+1):
            for wm in range(cfg.par.N_MASS_BINS+1):
                
                if has_stellar_mass[wz,wa,wm] == True:
                
                    source = m.add_point_source_collection()
                    
                    
                    nu = binned_stellar_nu
                    fnu = binned_stellar_fnu[counter,:]
                    nu,fnu = wavelength_compress(nu,fnu,df_nu)
                    
                    #reverse for hyperion
                    nu = nu[::-1]
                    fnu = fnu[::-1]

                    
                    #source luminosities
                    lum = np.array([stars_list[i].mass/constants.M_sun.cgs.value*constants.L_sun.cgs.value for i in stars_in_bin[(wz,wa,wm)]])
                    lum *= np.absolute(np.trapz(fnu,x=nu))
                    source.luminosity = lum
                    


                    for i in stars_in_bin[(wz,wa,wm)]:  totalmass += stars_list[i].mass
                    
                    #source positions
                    pos = np.zeros([len(stars_in_bin[(wz,wa,wm)]),3])
                    #for i in range(len(stars_in_bin[(wz,wa,wm)])): pos[i,:] = stars_list[i].positions
                    for i in range(len(stars_in_bin[(wz,wa,wm)])):
                        pos[i,:] = stars_list[stars_in_bin[(wz,wa,wm)][i]].positions

                    source.position=pos

                    #source spectrum
                    source.spectrum = (nu,fnu)
                                    
                    totallum += np.sum(source.luminosity)

                    
                    '''
                    if np.isnan(lum): 
                        print 'lum is a nan in point source collection addition. exiting now.'
                        sys.exit()
                    if np.isinf(lum): 
                        print 'lum is an inf in point source collection addition. exiting now.'
                        sys.exit()
                    '''
                counter+=1

                
    if cosmoflag == False: add_bulge_disk_stars(df_nu,binned_stellar_nu,binned_stellar_fnu,disk_fnu,bulge_fnu,stars_list,diskstars_list,bulgestars_list,m)

    m.set_sample_sources_evenly(True)

    t2=datetime.now()
    print ('[source_creation/add_binned_seds:] Execution time for point source collection adding = '+str(t2-t1))
    print ('[source_creation/add_binned_seds:] Total Luminosity of point source collection is: ',totallum)


    return m




def wavelength_compress(nu,fnu,df_nu):
    
    nu_inrange = np.logical_and(nu >= min(df_nu),nu <= max(df_nu))
    nu_inrange = np.where(nu_inrange == True)[0]
  
    compressed_nu = nu[nu_inrange]
    compressed_fnu = np.asarray(fnu)[nu_inrange]
    
  
    #get rid of all wavelengths below lyman limit
    dum_nu = compressed_nu*units.Hz
    dum_lam = constants.c.cgs/dum_nu
    dum_lam = dum_lam.to(units.angstrom)

    wll = np.where(dum_lam.value >= 912)[0] #where are lambda is above the lyman limit
    nu = nu[wll]
    fnu = fnu[wll]
   
   
    return compressed_nu,compressed_fnu



def BH_source_add(m,reg,df_nu,boost):

    print("--------------------------------\n")
    print("Adding Black Holes to Source List in source_creation\n")
    print("--------------------------------\n")
 

    try:
        nholes = reg["bhsed"].shape[0]

        #temporary wavelength compress just to get the length of the
        #compressed nu for a master array
        dumnu,dumfnu = wavelength_compress(reg["bhnu"].value,reg["bhsed"][0,:].value,df_nu)
        master_bh_fnu = np.zeros([nholes,len(dumnu)])
        
        holecounter = 0
        for i in range(nholes):  

            #don't create a BH luminsoity source if there's no luminosity since the SED will be nans/infs
            if reg["bhluminosity"][i].value > 0 :
                
                nu = reg["bhnu"].value
                fnu = reg["bhsed"][i,:].value#.tolist()
                nu,fnu = wavelength_compress(nu,fnu,df_nu)

                master_bh_fnu[i,:] = fnu
                
                if holecounter == 0:
                    fnu_compressed = np.zeros([nholes,len(nu)])
                fnu_compressed[i,:] = fnu

                #since the BH was added in a front end, we don't know
                #if the hole is in the actual cut out region of the yt
                #dataset.  so we need to filter out any holes that
                #might not be in the simulation domain.
                
                if ((reg["bhcoordinates"][i,0].in_units('kpc') <  (reg.domain_center[0].in_units('kpc')+(0.5*reg.domain_width[0].in_units('kpc'))))
                    and
                    (reg["bhcoordinates"][i,0].in_units('kpc') >  (reg.domain_center[0].in_units('kpc')-(0.5*reg.domain_width[0].in_units('kpc'))))
                    and
                    (reg["bhcoordinates"][i,1].in_units('kpc') <  (reg.domain_center[1].in_units('kpc')+(0.5*reg.domain_width[1].in_units('kpc'))))
                    and
                    (reg["bhcoordinates"][i,1].in_units('kpc') >  (reg.domain_center[1].in_units('kpc')-(0.5*reg.domain_width[1].in_units('kpc'))))
                    and
                    (reg["bhcoordinates"][i,2].in_units('kpc') <  (reg.domain_center[2].in_units('kpc')+(0.5*reg.domain_width[2].in_units('kpc'))))
                    and
                    (reg["bhcoordinates"][i,2].in_units('kpc') >  (reg.domain_center[2].in_units('kpc')-(0.5*reg.domain_width[2].in_units('kpc'))))
                ):

                    print('Boosting BH Coordinates and adding BH #%d to the source list now'%i)
                #the tolist gets rid of the array brackets
                    bh = m.add_point_source(luminosity = reg["bhluminosity"][i].value.tolist(), 
                                            spectrum = (nu,fnu),
                                            position = (reg["bhcoordinates"][i,:].in_units('cm').value-boost).tolist())
                else:
                    print('black hole #%d is not in the domain: rejecting adding it to the source list'%i)

                holecounter += 1

        dump_AGN_SEDs(nu,master_bh_fnu,reg["bhluminosity"].value)
    except:
        print('BH source creation failed.')
    #savefile = cfg.model.PD_output_dir+"/bh_sed.npz"
    #np.savez(savefile,nu = nu,fnu = master_bh_fnu,luminosity = ad["bhluminosity"].value)
    
