import random
import numpy as np
import pfh_readsnap
import parameters as par
from datetime import datetime
from astropy.table import Table
from astropy.io import ascii
import constants as const
import pdb
import sys

import fsps 


def new_sed_gen(sdir,snum):
    print 'reading in stars particles for SPS calculation'

    #NEW STARS
    new_stars_dict = pfh_readsnap.readsnap(sdir,snum,4)
    mass = new_stars_dict['m']*par.unit_mass*const.msun #g (as par.unit_mass is in msun)
    metals = new_stars_dict['z']
    positions = new_stars_dict['p']*par.unit_length*const.pc*1.e3 #cm (as par.unit_length is kpc)
    age = new_stars_dict['age'] #Gyr (per phopkins)

    median_metallicity = np.median(metals)
   

    if par.NEW_STARS_ONLY == False:
        #DISK STARS
        disk_stars_dict = pfh_readsnap.readsnap(sdir,snum,2)
        nstars_disk = len(disk_stars_dict['m'])

        mass = np.append(mass,disk_stars_dict['m']*par.unit_mass*const.msun) #g (as par.unit_mass is in msun)
        metals = np.append(metals,np.repeat(median_metallicity,nstars_disk))
        positions = np.concatenate((positions,disk_stars_dict['p']*par.unit_length*const.pc*1.e3)) #cm (as par.unit_length is kpc)
        age = np.append(age,np.repeat(par.disk_stars_age,nstars_disk))
        
        #BULGE STARS
        bulge_stars_dict = pfh_readsnap.readsnap(sdir,snum,3)
        nstars_bulge = len(bulge_stars_dict['m'])

        mass = np.append(mass,bulge_stars_dict['m']*par.unit_mass*const.msun) #g (as par.unit_mass is in msun)
        metals = np.append(metals,np.repeat(median_metallicity,nstars_bulge))
        positions = np.concatenate((positions,bulge_stars_dict['p']*par.unit_length*const.pc*1.e3)) #cm (as par.unit_length is kpc)
        age = np.append(age,np.repeat(par.bulge_stars_age,nstars_bulge))




    



    nstars = len(age)


   

    print 'generating stellar SEDs'

    #first figure out how many wavelengths there are
    sp = fsps.StellarPopulation(tage=age[0],imf_type=1,sfh=0)
    spec = sp.get_spectrum(tage=age[0])
    nu = 1.e8*const.c/spec[0]
    fnu = spec[1]

    nlam = len(nu)

    
    stellar_nu = np.zeros([nstars,nlam])
    stellar_fnu = np.zeros([nstars,nlam])
    
    
    
    #DEBUG DEBUG DEBUG 
    #for now we don't have a metallicity in the sps calculations
    print '========================='
    print 'WARNING: METALLICITIES NOT ACCOUNTED FOR IN STELLAR SEDS'
    print '========================'

    for i in range(nstars):
        
      



        #do sp? at the python line to find out what all possible
        #parameters there are
        sp = fsps.StellarPopulation(tage=age[i],imf_type=1,sfh=0)
        spec = sp.get_spectrum(tage=age[i])

        stellar_nu[i,:] = 1.e8*const.c/spec[0]
        stellar_fnu[i,:] = spec[1]

    

    return positions,mass,stellar_nu,stellar_fnu
