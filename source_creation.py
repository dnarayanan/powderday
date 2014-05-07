#source_creation.py
#purpose: to add new star sources, disks and bulge stars
#import parameters as par
import config as cfg
import pdb
import numpy as np
import constants as const
import operator
import SED_gen as sg
from datetime import datetime
from datetime import timedelta
import random


class Sed_Bins:
    def __init__(self,mass,metals,age):
        self.mass = mass
        self.metals = metals
        self.age = age



def add_super_simple_sed(stars_list,diskstars_list,bulgestars_list,m,lum,temp):

    print 'entering add_super_simple_sed function in source_creation'

    nstars = len(stars_list)
    nstars_disk = len(diskstars_list)
    nstars_bulge = len(bulgestars_list)

    for i in range(nstars):

        m.add_spherical_source(luminosity = lum, temperature = temp, radius = 10.*const.rsun,
                               position = stars_list[i].positions)
    
    
    for i in range(nstars_disk):
        
        m.add_spherical_source(luminosity = lum, temperature = temp, radius = 10.*const.rsun,
                               position = diskstars_list[i].positions)
    
    for i in range(nstars_bulge):
        
        m.add_spherical_source(luminosity = lum, temperature = temp, radius = 10.*const.rsun,
                               position = bulgestars_list[i].positions)
    


def add_newstars(df_nu,stellar_nu,stellar_fnu,disk_fnu,bulge_fnu,stars_list,diskstars_list,bulgestars_list,m):
    
        
    nstars = len(stars_list)
    print 'adding new stars to the grid'
    
    for i in range(nstars):
        
        nu = stellar_nu[:]
        fnu = stellar_fnu[i,:]
    

        nu_inrange = np.logical_and(nu >= min(df_nu),nu <= max(df_nu))
        nu_inrange = np.where(nu_inrange == True)[0]
        nu = nu[nu_inrange]
        fnu = fnu[nu_inrange]

        #reverse the arrays for hyperion
        nu = nu[::-1]
        fnu = fnu[::-1]

      

        lum = np.absolute(np.trapz(fnu,x=nu))*stars_list[i].mass/const.msun 
        lum *= const.lsun #to get in cgs

        #add new stars

        m.add_spherical_source(luminosity = lum,radius = 10.*const.rsun,spectrum = (nu,fnu),
                               position = stars_list[i].positions)
        
        
    if cfg.par.COSMOFLAG == False: add_bulge_disk_stars(df_nu,stellar_nu,stellar_fnu,disk_fnu,bulge_fnu,stars_list,diskstars_list,bulgestars_list,m)
            
    m.set_sample_sources_evenly(True)

    
def add_bulge_disk_stars(df_nu,stellar_nu,stellar_fnu,disk_fnu,bulge_fnu,stars_list,diskstars_list,bulgestars_list,m):
    print 'Non-Cosmological Simulation: Adding Disk and Bulge Stars:'
    

    nstars_disk = len(diskstars_list)
    nstars_bulge = len(bulgestars_list)
    
    fnu = disk_fnu[:]


    print 'adding disk stars to the grid: adding as a point source collection'   
    disksource = m.add_point_source_collection()
            
   
    

    nu = stellar_nu[:]

    nu_inrange = np.logical_and(nu >= min(df_nu),nu <= max(df_nu))
    nu_inrange = np.where(nu_inrange == True)[0]
    nu = nu[nu_inrange]
    fnu = fnu[nu_inrange]

    
    #reverse the arrays for hyperion
    nu = nu[::-1]
    fnu = fnu[::-1]
    


    
    disk_lum = np.absolute(np.trapz(fnu,x=nu))*diskstars_list[0].mass/const.msun
    #since stellar masses are in cgs, and we need them to be in msun - we
    #multiply by mass to get the *total* luminosity of the stellar
    #cluster since int(nu,fnu) is just the luminosity of a 1 Msun single star
    disk_lum *= const.lsun
    disksource.luminosity = np.repeat(disk_lum,nstars_disk)
    
    disk_pos = np.zeros([len(diskstars_list),3])
    for i in range(len(diskstars_list)): disk_pos[i,:] = diskstars_list[i].positions
    disksource.position=disk_pos
    
    disksource.spectrum = (nu,fnu)
            
    print 'adding bulge stars to the grid: adding as a point source collection'
    
    
    bulgesource = m.add_point_source_collection()
    
    fnu = bulge_fnu[:]
    fnu = fnu[nu_inrange]
    fnu = fnu[::-1]
    
    bulge_lum = np.absolute(np.trapz(fnu,x=nu))*bulgestars_list[0].mass/const.msun
    bulge_lum *= const.lsun
    bulgesource.luminosity = np.repeat(bulge_lum,nstars_bulge)
    
    bulge_pos = np.zeros([len(bulgestars_list),3])
    for i in range(len(bulgestars_list)): bulge_pos[i,:] = bulgestars_list[i].positions
    bulgesource.position=bulge_pos
    
    fnu = disk_fnu[:]
    fnu = fnu[::-1]
    fnu = fnu[nu_inrange]
    bulgesource.spectrum = (nu,fnu)


    

def add_binned_seds(df_nu,stars_list,diskstars_list,bulgestars_list,m):
    
    #calculate max and min ages
    minimum_age = 15 #Gyr - obviously too high of a number
    maximum_age = 0 #Gyr


    #calculate the minimum and maximum luminosity
    minimum_mass = 1e15*const.msun #msun - some absurdly large value for a single stellar cluster
    maximum_mass = 0 #msun

    #calculate the minimum and maximum stellar metallicity
    minimum_metallicity = 1.e5 #some absurdly large metallicity
    maximum_metallicity = 0

    nstars = len(stars_list)
    for i in range(nstars):
        if stars_list[i].metals[0] < minimum_metallicity: minimum_metallicity = stars_list[i].metals[0]
        if stars_list[i].metals[0] > maximum_metallicity: maximum_metallicity = stars_list[i].metals[0]
        
        if stars_list[i].mass < minimum_mass: minimum_mass = stars_list[i].mass
        if stars_list[i].mass > maximum_mass: maximum_mass = stars_list[i].mass

        if stars_list[i].age < minimum_age: minimum_age = stars_list[i].age
        if stars_list[i].age > maximum_age: maximum_age = stars_list[i].age

      


    #define the metallicity bins: we do this in log since there can be quite a spread in metallicities

    delta_metallicity =  (np.log10(maximum_metallicity)-np.log10(minimum_metallicity))/cfg.par.N_METAL_BINS
    metal_bins = np.arange(np.log10(minimum_metallicity),
                           np.log10(maximum_metallicity),
                           delta_metallicity)
    #add on the maximum metallicity bin
    metal_bins = np.append(metal_bins,metal_bins[-1]+delta_metallicity)
    metal_bins = 10.**metal_bins


    #define the age bins (linearly)
    delta_age = (maximum_age-minimum_age)/cfg.par.N_STELLAR_AGE_BINS
    age_bins = np.arange(minimum_age,maximum_age,delta_age)
    
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
        

    print 'mass_bins = ',mass_bins
    print 'metal_bins = ',metal_bins
    print 'age_bins = ',age_bins


    #has_stellar_mass is a 3D boolean array that's [wz,wa,wm] big and
    #says whether or not that bin is being used downstream for
    #creating a point source collection (i.e. that it actually has at
    #least one star cluster that falls into it)
    has_stellar_mass = np.zeros([cfg.par.N_METAL_BINS+1,cfg.par.N_STELLAR_AGE_BINS+1,cfg.par.N_MASS_BINS+1],dtype=bool)


    stars_in_bin = {} #this will be a dictionary that holds the list
    #of star particles that go in every [wz,wa,wm]
    #group.  The keys will be tuples that hold a
    #(wz,wa,wm) set that we will then use later to
    #speed up adding sources.
    
    for i in range(nstars):
        
        wz = find_nearest(metal_bins,stars_list[i].metals[0])
        wa = find_nearest(age_bins,stars_list[i].age)
        wm = find_nearest(mass_bins,stars_list[i].mass)
        
        stars_list[i].sed_bin = [wz,wa,wm]

        has_stellar_mass[wz,wa,wm] = True

        if (wz,wa,wm) in stars_in_bin:
            stars_in_bin[(wz,wa,wm)].append(i)
        else:
            stars_in_bin[(wz,wa,wm)] = [i]

   



    



    
    print 'assigning stars to SED bins'
    sed_bins_list=[]
    sed_bins_list_has_stellar_mass = []


   

    for wz in range(cfg.par.N_METAL_BINS+1):
        for wa in range(cfg.par.N_STELLAR_AGE_BINS+1):
            for wm in range(cfg.par.N_MASS_BINS+1):
                sed_bins_list.append(Sed_Bins(mass_bins[wm],metal_bins[wz],age_bins[wa]))
                if has_stellar_mass[wz,wa,wm] == True:
                    sed_bins_list_has_stellar_mass.append(Sed_Bins(mass_bins[wm],metal_bins[wz],age_bins[wa]))

   

    #sed_bins_list is a list of Sed_Bins objects that have the
    #information about what mass bin, metal bin and age bin they
    #correspond to.  It is unnecessary, and heavy computational work
    #to re-create the SED for each of these bins - rather, we can just
    #calculate the SED for the bins that have any actual stellar mass.
            
    print 'Running SPS for Binned SEDs'
    print 'calculating the SEDs for ',len(sed_bins_list_has_stellar_mass),' bins'
    binned_stellar_nu,binned_stellar_fnu_has_stellar_mass,disk_fnu,bulge_fnu = sg.allstars_sed_gen(sed_bins_list_has_stellar_mass,diskstars_list,bulgestars_list)

    #since the binned_stellar_fnu_has_stellar_mass is now
    #[len(sed_bins_list_has_stellar_mass),nlam)] big, we need to
    #transform it back to the a larger array.  this is an ugly loop
    #that could probably be prettier...but whatever.  this saves >an
    #order of magnitude in time in SED gen.  
    nlam = binned_stellar_nu.shape[0]
    binned_stellar_fnu = np.zeros([len(sed_bins_list),nlam])

    counter = 0
    counter_has_stellar_mass = 0
    for wz in range(cfg.par.N_METAL_BINS+1):
        for wa in range(cfg.par.N_STELLAR_AGE_BINS+1):
            for wm in range(cfg.par.N_MASS_BINS+1):
                if has_stellar_mass[wz,wa,wm] == True:
                    binned_stellar_fnu[counter,:] = binned_stellar_fnu_has_stellar_mass[counter_has_stellar_mass,:]
                    counter_has_stellar_mass += 1 
                counter+=1


    #now binned_stellar_nu and binned_stellar_fnu are the SEDs for the bins in order of wz, wa, wm 
    
    #create the point source collections: we loop through the bins and
    #see what star particles correspond to these.  if any do, then we
    #add them to a list, and create a point source collection out of
    #these

    nu = binned_stellar_nu
    nu_inrange = np.logical_and(nu >= min(df_nu),nu <= max(df_nu))
    nu_inrange = np.where(nu_inrange == True)[0]
    
    nu = binned_stellar_nu[nu_inrange]
    nu = nu[::-1]

    print 'adding point source collections'
    t1=datetime.now()

    counter=0
    for wz in range(cfg.par.N_METAL_BINS+1):
        for wa in range(cfg.par.N_STELLAR_AGE_BINS+1):
            for wm in range(cfg.par.N_MASS_BINS+1):
                
                if has_stellar_mass[wz,wa,wm] == True:
                
                    source = m.add_point_source_collection()
                    
                    fnu = binned_stellar_fnu[counter,:]
                    fnu = fnu[nu_inrange]
                    fnu = fnu[::-1]
                    pos = np.zeros([len(stars_in_bin[(wz,wa,wm)]),3])


                    if cfg.par.SUPER_SIMPLE_SED == False:
                        
                        lum = np.absolute(np.trapz(fnu,x=nu))*mass_bins[wm]/const.msun*const.lsun
                        source.luminosity = np.repeat(lum,len(stars_in_bin[(wz,wa,wm)]))
                        for i in range(len(stars_in_bin[(wz,wa,wm)])): pos[i,:] = stars_list[i].positions
                        source.position=pos
                        source.spectrum = (nu,fnu)

                    else:

                        lum = 1.e3*const.lsun
                        source.luminosity = np.repeat(lum,len(stars_in_bin[(wz,wa,wm)]))
                        for i in range(len(stars_in_bin[(wz,wa,wm)])): pos[i,:] = stars_list[i].positions
                        source.position=pos
                        source.temperature = 1.e4
                        print 'adding super simple SED'


                
                counter+=1

                
    if cfg.par.COSMOFLAG == False: add_bulge_disk_stars(df_nu,binned_stellar_nu,binned_stellar_fnu,disk_fnu,bulge_fnu,stars_list,diskstars_list,bulgestars_list,m)

    m.set_sample_sources_evenly(True)

    t2=datetime.now()
    print 'Execution time for point source collection adding = '+str(t2-t1)


    
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    
    return idx


