from __future__ import print_function
import numpy as np
import powderday.config as cfg
import pdb

import astropy.units as u
import astropy.constants as constants
from astropy import cosmology as cosmo

import fsps 
from datetime import datetime
from powderday.grid_construction import stars_coordinate_boost

from multiprocessing import Pool
from functools import partial
from scipy.integrate import simps
from powderday.nebular_emission.cloudy_tools import calc_LogU
from powderday.nebular_emission.cloudy_tools import cmdf
from powderday.analytics import logu_diagnostic,dump_emlines
from powderday.nebular_emission.cloudy_model import get_nebular

#this is required to keep the reg as a strong reference.  for some
#reason in the star_list.append in star_list_gen, reg otherwise gets
#garbage collected.
import gc
gc.set_threshold(0)

# Lazily initialize FSPS
sp = None

class Stars:
    def __init__(self,mass,metals,positions,age,sed_bin=[-1,-1,-1],lum=-1,fsps_zmet=20,all_metals=[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]):
        self.mass = mass
        self.metals = metals
        self.positions = positions
        self.age = age
        self.sed_bin = sed_bin
        self.lum = lum
        self.fsps_zmet = fsps_zmet
        self.all_metals = all_metals

    def info(self):
        return(self.mass,self.metals,self.positions,self.age,self.sed_bin,self.lum,self.fsps_zmet)


def star_list_gen(boost,dx,dy,dz,reg,ds):
    print ('[SED_gen/star_list_gen]: reading in stars particles for SPS calculation')
    mass = reg["starmasses"].value
    positions = reg["starcoordinates"].value
    age = reg["stellarages"].value
    nstars = len(reg["stellarages"].value) 
    el = ['He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe' ]

    try:                                                                                                                                                                
        metals = np.zeros((nstars,11))-10.0
        for i in range(11):
            if i == 0:
                el_str = ""
            else:
                el_str = "_"+el[i-1]
            metals[:, i] = reg["starmetals"+el_str].value
    except:
        metals = reg["starmetals"].value
    print ('number of new stars =',nstars)
    
    #calculate the fsps interpolated metallicity

    #if the metallicity has many fields, and not just global
    #metallicity then just extract the global metallicity
    if metals.ndim > 1:
        metals_tot = metals[:,0]
    else:
        metals_tot = metals

    print ('[SED_gen/star_list_gen:] Manually increasing the newstar metallicities by: ',cfg.par.Z_init)
    metals_tot += cfg.par.Z_init
    
    #ADVANCED FEATURE - if force_stellar_metallcities or force_stellar_ages are set, then we set to those values
    if cfg.par.FORCE_STELLAR_AGES:
        print ("[SED_GEN/stars_list_gen:]  FORCE_STELLAR_AGES is set to True: setting all stars to age: %e Gyr"%cfg.par.FORCE_STELLAR_AGES_VALUE)
        age = np.repeat(cfg.par.FORCE_STELLAR_AGES_VALUE,nstars)

    if cfg.par.FORCE_STELLAR_METALLICITIES:
        print ("[SED_GEN/stars_list_gen:]  FORCE_STELLAR_METALLICITIES is set to True: setting all stars to metallicity: %e "%cfg.par.FORCE_STELLAR_METALLICITIES_VALUE)
        metals_tot = np.repeat(cfg.par.FORCE_STELLAR_METALLICITIES_VALUE,nstars)



    zmet = fsps_metallicity_interpolate(metals_tot)
    #mwd(zmet,mass,'zmet_distribution.png')

    #print '[SED_gen/star_list_gen: ] fsps zmet codes:',zmet

    #create the stars_list full of Stars objects
    stars_list = []
    
    if metals.ndim > 1:
        for i in range(nstars):
            stars_list.append(Stars(mass[i],metals_tot[i],positions[i],age[i],fsps_zmet=zmet[i],all_metals = metals[i]))
    else:
        for i in range(nstars):
            stars_list.append(Stars(mass[i],metals_tot[i],positions[i],age[i],fsps_zmet=zmet[i]))
    
    #boost stellar positions to grid center
    print ('boosting new stars to coordinate center')
    stars_list = stars_coordinate_boost(stars_list,boost)

    #orig_stars_list_len = len(stars_list)
    
    #ASSIGN DISK AND BULGE STARS - note, if these don't exist, it will
    #just make empty lists

   

    
    bulgestars_list = []
    diskstars_list = []

    
    #in principle, we should just be able to do the following blocks
    #if the particle types exist. the issue is that different groups
    #use PartType2 and 3 as 'filler' particle types, so they may exist
    #even if they don't correspond to disk/bulge stars.


    if ds.cosmological_simulation == False:

        #Disk Stars

        if ("diskstarcoordinates") in ds.derived_field_list:
            
            disk_positions = reg[("diskstarcoordinates")].value
            disk_masses =  reg[("diskstarmasses")].value
            nstars_disk = len(disk_masses)
     
            #create the disk_list full of DiskStars objects
            for i in range(nstars_disk):
                diskstars_list.append(Stars(disk_masses[i],cfg.par.solar,disk_positions[i],cfg.par.disk_stars_age))

            print ('boosting disk stars to coordinate center')    
            diskstars_list = stars_coordinate_boost(diskstars_list,boost)

        #orig_disk_stars_list_len = nstars_disk
            
       
        #Bulge Stars


        if ("bulgestarcoordinates") in ds.derived_field_list:
            bulge_positions = reg[("bulgestarcoordinates")].value
            bulge_masses =  reg[("bulgestarmasses")].value
            nstars_bulge = len(bulge_masses)
            
            #create the bulge_list full of BulgeStars objects
            
            for i in range(nstars_bulge):
                bulgestars_list.append(Stars(bulge_masses[i],cfg.par.solar,bulge_positions[i],cfg.par.bulge_stars_age))
                

            print ('boosting bulge stars to coordinate center')
            bulgestars_list = stars_coordinate_boost(bulgestars_list,boost)


    #EXPERIMENTAL FEATURES
    if cfg.par.SOURCES_IN_CENTER == True:
        for i in range(nstars):
            stars_list[i].positions[:] =  np.array([0,0,0])
        if ("bulgestarcoordinates") in ds.derived_field_list:
            for i in range(nstars_bulge):
                bulgestars_list[i].positions[:] =  np.array([0,0,0])
            for i in range(nstars_disk):
                diskstars_list[i].positions[:] = np.array([0,0,0])

    if cfg.par.SOURCES_RANDOM_POSITIONS == True:
        print ("================================")
        print ("SETTING SOURCES TO RANDOM POSITIONS")
        print ("================================")
        for i in range(nstars):
            xpos,ypos,zpos = np.random.uniform(-0.9*dx/2.,0.9*dx/2.),np.random.uniform(-0.9*dy/2.,0.9*dy/2.),np.random.uniform(-0.9*dz/2.,0.9*dz/2.)
            stars_list[i].positions[:] = np.array([xpos,ypos,zpos])

        if ("bulgestarcoordinates") in ds.derived_field_list:
            for i in range(nstars_bulge):
                xpos,ypos,zpos = np.random.uniform(-0.9*dx/2.,0.9*dx/2.),np.random.uniform(-0.9*dy/2.,0.9*dy/2.),np.random.uniform(-0.9*dz/2.,0.9*dz/2.)
                bulgestars_list[i].positions[:] = np.array([xpos,ypos,zpos])
            for i in range(nstars_disk):
                xpos,ypos,zpos = np.random.uniform(-0.9*dx/2.,0.9*dx/2.),np.random.uniform(-0.9*dy/2.,0.9*dy/2.),np.random.uniform(-0.9*dz/2.,0.9*dz/2.)
                diskstars_list[i].positions[:] = np.array([xpos,ypos,zpos])

    return stars_list,diskstars_list,bulgestars_list,reg



def allstars_sed_gen(stars_list,cosmoflag,sp):


    #NOTE this part is just for the gadget simulations - this will
    #eventually become obviated as it gets passed into a function to
    #populate the stars_list with objects as we start to feed in new
    #types of simulation results.

    nstars = len(stars_list)
    
    #get just the wavelength array
    sp.params["tage"] = stars_list[0].age
    sp.params["imf_type"] = cfg.par.imf_type
    sp.params["pagb"] = cfg.par.pagb
    sp.params["sfh"] = 0
    sp.params["zmet"] = stars_list[0].fsps_zmet
    sp.params["add_neb_emission"] = cfg.par.add_neb_emission
    sp.params["add_agb_dust_model"] = cfg.par.add_agb_dust_model
    sp.params['gas_logu'] = cfg.par.gas_logu
    if cfg.par.FORCE_gas_logz == False:
        sp.params['gas_logz'] = np.log10(stars_list[0].metals/cfg.par.solar)
    else:
        sp.params['gas_logz'] = cfg.par.gas_logz

        '''
    sp = fsps.StellarPopulation(tage=stars_list[0].age,imf_type=cfg.par.imf_type,pagb = cfg.par.pagb,sfh=0,zmet=stars_list[0].fsps_zmet,
                                add_neb_emission = cfg.par.add_neb_emission, add_agb_dust_model=cfg.par.add_agb_dust_model)
                                '''
    spec = sp.get_spectrum(tage=stars_list[0].age,zmet=stars_list[0].fsps_zmet)
    nu = 1.e8*constants.c.cgs.value/spec[0]
    nlam = len(nu)

    nprocesses = np.min([cfg.par.n_processes,len(stars_list)]) #the pool.map will barf if there are less star bins than process threads

    #initialize the process pool and build the chunks
    p = Pool(processes = nprocesses)
    nchunks = nprocesses


    chunk_start_indices = []
    chunk_start_indices.append(0) #the start index is obviously 0


    #this should just be int(nstars/nchunks) but in case nstars < nchunks, we need to ensure that this is at least  1
    delta_chunk_indices = np.max([int(nstars / nchunks),1]) 
    print ('delta_chunk_indices = ',delta_chunk_indices)
    
    for n in range(1,nchunks):
        chunk_start_indices.append(chunk_start_indices[n-1]+delta_chunk_indices)

    '''
    chunk_start_indices = list(np.fix(np.arange(0,nstars,np.fix(nstars/nchunks))))
    #because this can result in too many chunks sometimes given the number of processors:
    chunk_start_indices = chunk_start_indices[0:nchunks]
    '''
    print ('Entering Pool.map multiprocessing for Stellar SED generation')
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
    print ('Execution time for SED generation in Pool.map multiprocessing = '+str(t2-t1))

    
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



    if cosmoflag == False:

        #calculate the SED for disk stars; note, this gets calculated
        #whether or not disk stars actually exist.  if they don't exist,
        #bogus values for the disk age and metallicity are assigned based
        #on whatever par.disk_stars_age and metallicity are.  it's no big
        #deal since these SEDs don't end up getting added to the model in
        #source_creation.  
        
        #note, even if there are no disk/bulge stars, these are still
        #created since they're completely based on input parameters in
        #parameters_master.  they just won't get used at a later point
        #as there will be no disk/bulge star positions to add them to.

        #dust_tesc is an absolute value (not relative to min star age) as the ages of these stars are input by the user

        # Load in the metallicity legend
        fsps_metals = np.loadtxt(cfg.par.metallicity_legend)

        sp.params["tage"] = cfg.par.disk_stars_age
        sp.params["imf_type"] = cfg.par.imf_type
        sp.params["imf1"] = cfg.par.imf1
        sp.params["imf2"] = cfg.par.imf2
        sp.params["imf3"] = cfg.par.imf3
        sp.params["pagb"] = cfg.par.pagb
        sp.params["sfh"] = 0
        sp.params["zmet"] = cfg.par.disk_stars_metals
        sp.params["add_neb_emission"] = cfg.par.add_neb_emission
        sp.params["add_agb_dust_model"] = cfg.par.add_agb_dust_model
        sp.params['gas_logu'] = cfg.par.gas_logu
        if cfg.par.FORCE_gas_logz == False:
            sp.params['gas_logz'] = np.log10(fsps_metals[cfg.par.disk_stars_metals]/cfg.par.solar)
        else:
            sp.params['gas_logz'] = cfg.par.gas_logz

        spec = sp.get_spectrum(tage=cfg.par.disk_stars_age,zmet=cfg.par.disk_stars_metals)
        disk_fnu = spec[1]
        
        #calculate the SED for bulge stars
        sp.params["tage"] = cfg.par.bulge_stars_age
        sp.params["imf_type"] = cfg.par.imf_type
        sp.params["imf1"] = cfg.par.imf1
        sp.params["imf2"] = cfg.par.imf2
        sp.params["imf3"] = cfg.par.imf3
        sp.params["pagb"] = cfg.par.pagb
        sp.params["sfh"] = 0
        sp.params["zmet"] = cfg.par.bulge_stars_metals
        sp.params["add_neb_emission"] = cfg.par.add_neb_emission
        sp.params["add_agb_dust_model"] = cfg.par.add_agb_dust_model
        sp.params['gas_logu'] = cfg.par.gas_logu
        if cfg.par.FORCE_gas_logz == False:
            sp.params['gas_logz'] = np.log10(fsps_metals[cfg.par.bulge_stars_metals]/cfg.par.solar)
        else:
            sp.params['gas_logz'] = cfg.par.gas_logz


        spec = sp.get_spectrum(tage=cfg.par.bulge_stars_age,zmet=cfg.par.bulge_stars_metals)
        bulge_fnu = spec[1]
    

    else: #we have a cosmological simulation

        disk_fnu = []
        bulge_fnu = []
        

    total_lum_in_sed_gen = 0.
    for i in range(stellar_fnu.shape[0]):
        total_lum_in_sed_gen += np.absolute(np.trapz(stellar_fnu[i,:],x=nu))

    print ('[SED_gen: ] total_lum_in_sed_gen = ',total_lum_in_sed_gen)

    #return positions,disk_positions,bulge_positions,mass,stellar_nu,stellar_fnu,disk_masses,disk_fnu,bulge_masses,bulge_fnu
    return stellar_nu,stellar_fnu,disk_fnu,bulge_fnu


def newstars_gen(stars_list):
    global sp
    if sp is None:
        sp = fsps.StellarPopulation()

    #the newstars (particle type 4; so, for cosmological runs, this is all
    #stars) are calculated in a separate function with just one argument so that it is can be fed 
    #into pool.map for multithreading.
    #sp = fsps.StellarPopulation()
    sp.params["tage"] = stars_list[0].age
    sp.params["imf_type"] = cfg.par.imf_type
    sp.params["pagb"] = cfg.par.pagb
    sp.params["sfh"] = 0
    sp.params["zmet"] = stars_list[0].fsps_zmet
    sp.params["add_neb_emission"] = False
    sp.params["add_agb_dust_model"] = cfg.par.add_agb_dust_model
    sp.params['gas_logu'] = cfg.par.gas_logu

    if cfg.par.FORCE_gas_logz == False:
        sp.params['gas_logz'] = np.log10(stars_list[0].metals/cfg.par.solar)
    else:
        sp.params['gas_logz'] = cfg.par.gas_logz

    #first figure out how many wavelengths there are
    
    spec = sp.get_spectrum(tage=stars_list[0].age,zmet=stars_list[0].fsps_zmet)
    nu = 1.e8*constants.c.cgs.value/spec[0]

    nlam = len(nu)

    stellar_nu = np.zeros([nlam])
    stellar_fnu = np.zeros([len(stars_list),nlam])
    
  
    minage = 13 #Gyr
    for i in range(len(stars_list)): 
        if stars_list[i].age < minage:
            minage = stars_list[i].age

    tesc_age = np.log10((minage+cfg.par.birth_cloud_clearing_age)*1.e9)


    # Get the number of ionizing photons from SED

    #calculate the SEDs for new stars
    for i in range(len(stars_list)):
        
        sp.params["tage"] = stars_list[i].age
        sp.params["imf_type"] = cfg.par.imf_type
        sp.params["imf1"] = cfg.par.imf1
        sp.params["imf2"] = cfg.par.imf2
        sp.params["imf3"] = cfg.par.imf3
        sp.params["pagb"] = cfg.par.pagb
        sp.params["sfh"] = 0
        sp.params["zmet"] = stars_list[i].fsps_zmet
        sp.params["add_neb_emission"] = False
        sp.params["add_agb_dust_model"] = cfg.par.add_agb_dust_model

        if cfg.par.FORCE_gas_logz == False:
            LogZ = np.log10(stars_list[i].metals/cfg.par.solar)
        else:
            LogZ = cfg.par.gas_logz
        
        if cfg.par.CF_on == True:
            sp.params["dust_type"] = 0
            sp.params["dust1"] = 1
            sp.params["dust2"] = 0
            sp.params["dust_tesc"] = tesc_age


        #sp = fsps.StellarPopulation(tage=stars_list[i].age,imf_type=2,sfh=0,zmet=stars_list[i].fsps_zmet)
        spec = sp.get_spectrum(tage=stars_list[i].age,zmet=stars_list[i].fsps_zmet)
        f = spec[1]
        
        pagb = cfg.par.add_pagb_stars and stars_list[i].age >= cfg.par.PAGB_min_age and stars_list[i].age <= cfg.par.PAGB_max_age
        
        #Only including particles below the maximum age limit for calulating nebular emission
        if (cfg.par.add_neb_emission or cfg.par.use_cmdf) and (stars_list[i].age <= cfg.par.HII_max_age or pagb):
            
            cluster_mass, num_clusters = cmdf(stars_list[i].mass/constants.M_sun.cgs.value,int(cfg.par.cmdf_bins),cfg.par.cmdf_min_mass,
                                                cfg.par.cmdf_max_mass, cfg.par.cmdf_beta)
            f = np.zeros(nlam)
            line_em = np.zeros([cfg.par.cloudy_nlam])
            for j in range(len(cluster_mass)):
                num_HII_clusters = num_clusters[j]
                neb_file_output = cfg.par.NEB_DEBUG
            
                sp.params["add_neb_emission"] = False
                
                if cfg.par.add_neb_emission:
                    spec = sp.get_spectrum(tage=stars_list[i].age,zmet=stars_list[i].fsps_zmet)
                    print (j, 10**cluster_mass[j], num_HII_clusters, stars_list[i].age, LogZ)
                    if cfg.par.FORCE_gas_logu:
                        alpha = 2.5e-13*((cfg.par.HII_T/(10**4))**(-0.85))
                        LogU = cfg.par.gas_logu
                        LogQ = np.log10((10 ** (3*LogU))*(36*np.pi*(constants.c.cgs.value**3))/((alpha**2)*cfg.par.HII_nh))
                        Rs = ((3*(10 ** LogQ))/(4*np.pi*(cfg.par.HII_nh**2)*alpha))**(1./3.)
                    else:
                        LogQ, Rs, LogU = calc_LogU(1.e8*constants.c.cgs.value/spec[0], spec[1]*constants.L_sun.cgs.value,
                                                cfg.par.HII_nh, cfg.par.HII_T, cfg.par.HII_escape_fraction, mstar=10**cluster_mass[j])
                
                    alpha = 2.5e-13*((cfg.par.HII_T/(10**4))**(-0.85))
                    LogU = LogU + cfg.par.gas_logu_init
                    LogQ = np.log10((10 ** (3*LogU))*(36*np.pi*(constants.c.cgs.value**3))/((alpha**2)*cfg.par.HII_nh))
                    Rs = ((3*(10 ** LogQ))/(4*np.pi*(cfg.par.HII_nh**2)*alpha))**(1./3.)
                    
                    Rin = cfg.par.Rinner_per_Rs*Rs
                    
                    if cfg.par.FORCE_logq:
                        LogQ = cfg.par.source_logq

                    if cfg.par.FORCE_inner_radius:
                        Rin = cfg.par.inner_radius
               
                    if neb_file_output:
                        if cfg.par.use_cloudy_tables:
                            Rin = 1.e19   # Rinner is fixed at 1.e19 cm for lookup tables
                        logu_diagnostic(LogQ, LogU, LogZ, Rin, 10**cluster_mass[j], num_HII_clusters, stars_list[i].age, append=True)
                        neb_file_output = False

                        sp.params['gas_logu'] = LogU
                        sp.params['gas_logz'] = LogZ
                        sp.params["add_neb_emission"] = True 
                
                    if cfg.par.use_cloudy_tables:
                        lam_neb, spec_neb = sp.get_spectrum(tage=stars_list[i].age, zmet=stars_list[i].fsps_zmet)
                        line_lum = sp.emline_luminosity
                        wave_line = sp.emline_wavelengths
                    else:
                        try:
                            # Calculating ionizing photons again but for 1 Msun in order to scale the output for FSPS
                            LogQ_1, Rin_1, LogU_1 = calc_LogU(1.e8 * constants.c.cgs.value / spec[0],
                                                            spec[1] * constants.L_sun.cgs.value, cfg.par.HII_nh,
                                                            cfg.par.HII_T, cfg.par.HII_escape_fraction)
                               
                            
                            spec_neb, wave_line, line_lum = get_nebular(spec[0], spec[1], cfg.par.HII_nh, LogQ, Rin, LogU, LogZ, LogQ_1, stars_list[i].all_metals, 
                                                    Dust=cfg.par.neb_dust, abund=cfg.par.neb_abund, useq = cfg.par.use_Q, Pagb=pagb, clean_up = cfg.par.cloudy_cleanup)
                        except ValueError as err:
                            print ("WARNING: CLOUDY run was unsucessful. Using lookup tables for nebular emission")
                            lam_neb, spec_neb = sp.get_spectrum(tage=stars_list[i].age, zmet=stars_list[i].fsps_zmet)
                            line_lum = sp.emline_luminosity
                            wave_line = sp.emline_wavelengths
                
                else:
                    lam_neb, spec_neb = sp.get_spectrum(tage=stars_list[i].age, zmet=stars_list[i].fsps_zmet)

                weight = num_HII_clusters*(10**cluster_mass[j])/(stars_list[i].mass/constants.M_sun.cgs.value)
                f = f + spec_neb*weight
                if cfg.par.add_neb_emission and cfg.par.dump_emlines:
                    line_em = line_em + line_lum*weight
        
            if cfg.par.add_neb_emission and cfg.par.dump_emlines:
                #the stellar population returns the calculation in units of Lsun/1 Msun: https://github.com/dfm/python-fsps/issues/117#issuecomment-546513619
                line_em = line_em * ((stars_list[i].mass*u.g).to(u.Msun).value) * (3.839e33)  # Units: ergs/s
                line_em = np.append(line_em, stars_list[i].age)
                dump_emlines(wave_line, line_em)

        stellar_nu[:] = 1.e8*constants.c.cgs.value/spec[0]
        stellar_fnu[i,:] = f

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
