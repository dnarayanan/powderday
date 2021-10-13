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
from itertools import repeat
from scipy.integrate import simps

from powderday.nebular_emission.cloudy_tools import calc_LogQ, cmdf, get_nearest,convert_metals
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
    mass = reg["star","masses"].value
    positions = reg["star","coordinates"].value
    age = reg["stellar","ages"].value
    nstars = len(reg["stellar","ages"].value) 
    el = ['He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe' ]

    try:                                                                                                                                                                
        metals = np.zeros((nstars,11))-10.0
        for i in range(11):
            if i == 0:
                el_str = ""
            else:
                el_str = "_"+el[i-1]
            metals[:, i] = reg["star","metals"+el_str].value
    except:
        metals = reg["star","metals"].value
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

        if ("diskstar","coordinates") in ds.derived_field_list:
            
            disk_positions = reg[("diskstar","coordinates")].value
            disk_masses =  reg[("diskstar","masses")].value
            nstars_disk = len(disk_masses)
     
            #create the disk_list full of DiskStars objects
            for i in range(nstars_disk):
                diskstars_list.append(Stars(disk_masses[i],cfg.par.solar,disk_positions[i],cfg.par.disk_stars_age))

            print ('boosting disk stars to coordinate center')    
            diskstars_list = stars_coordinate_boost(diskstars_list,boost)

        #orig_disk_stars_list_len = nstars_disk
            
       
        #Bulge Stars


        if ("bulgestarcoordinates") in ds.derived_field_list:
            bulge_positions = reg[("bulgestar","coordinates")].value
            bulge_masses =  reg[("bulgestar","masses")].value
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
        if ("bulgestar","coordinates") in ds.derived_field_list:
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

        if ("bulgestar","coordinates") in ds.derived_field_list:
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
    mfrac = np.zeros(nstars)
    # see newstars_gen() for info on mfrac
    star_counter=0
    for i in range(nchunks):
        fnu_list = chunk_sol[i][0] #this is a list of the stellar_fnu's returned by that chunk
        #chunk_sol now returns two things for each star/star bin: the spectrum and the associated surviving stellar mass fraction for that SSP 
        mfrac_list = chunk_sol[i][1]
        for j in range(len(fnu_list)):
            mfrac[star_counter] = mfrac_list[j] 
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
    return stellar_nu,stellar_fnu,disk_fnu,bulge_fnu, mfrac


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

    #first figure out how many wavelengths there are
    
    spec = sp.get_spectrum(tage=stars_list[0].age,zmet=stars_list[0].fsps_zmet)
    nu = 1.e8*constants.c.cgs.value/spec[0]

    nlam = len(nu)

    stellar_nu = np.zeros([nlam])
    stellar_fnu = np.zeros([len(stars_list),nlam])
    mfrac = np.zeros([len(stars_list)])
  
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

        if cfg.par.CF_on == True:
            sp.params["dust_type"] = 0
            sp.params["dust1"] = 1
            sp.params["dust2"] = 0
            sp.params["dust_tesc"] = tesc_age

        spec_noneb = sp.get_spectrum(tage=stars_list[i].age,zmet=stars_list[i].fsps_zmet)
        f = spec_noneb[1]
        #NOTE: FSPS SSP/CSP spectra are scaled by *formed* mass, not current mass. i.e., the SFHs of the SSP/CSP are normalized such that 1 solar mass
        #is formed over the history. This means that the stellar spectra are normalized by the integral of the SFH =/= current 
        #(surviving, observed, etc.) stellar mass. In simulations, we only know the current star particle mass. To get the formed mass for an SSP,
        #we generate the surviving mass fraction (sp.stellar_mass) to extrapolate the initial mass from the current mass, metallicity, and age
        #this 'mfrac' is used to scale the FSPS SSP luminosities in source_creation
        mfrac[i] = sp.stellar_mass

        pagb = cfg.par.add_pagb_stars and cfg.par.PAGB_min_age <= stars_list[i].age <= cfg.par.PAGB_max_age
        young_star = cfg.par.add_young_stars and cfg.par.HII_min_age <= stars_list[i].age <= cfg.par.HII_max_age

        if (cfg.par.add_neb_emission or cfg.par.use_cmdf) and (young_star or pagb):

            # Cluster Mass Distribution Funtion is used only when the star particle's mass is gretaer than the maximum cluster mass and use_cmdf is True. 

            if stars_list[i].mass/constants.M_sun.cgs.value > 10**cfg.par.cmdf_max_mass and cfg.par.use_cmdf:
                cluster_mass, num_clusters = cmdf(stars_list[i].mass/constants.M_sun.cgs.value,int(cfg.par.cmdf_bins),cfg.par.cmdf_min_mass,
                        cfg.par.cmdf_max_mass, cfg.par.cmdf_beta, rescale_masses=cfg.par.cmdf_rescale)
            
            else:
                cluster_mass = [np.log10(stars_list[i].mass/constants.M_sun.cgs.value)]
                num_clusters = [1]

            f = np.zeros(nlam)
            cloudy_nlam = len(np.genfromtxt(cfg.par.pd_source_dir + "/powderday/nebular_emission/data/refLines.dat", delimiter=','))
            line_em = np.zeros([cloudy_nlam])
            
            for j in range(len(cluster_mass)):
                num_HII_clusters = num_clusters[j]
                neb_file_output = cfg.par.NEB_DEBUG

                sp.params["add_neb_emission"] = False
                if cfg.par.add_neb_emission:
                    # id_val = 0, 1, 2 for young stars, Post-AGB star and AGNs respectively.
                    if young_star:
                        id_val = 0
                        Rinner_per_Rs = cfg.par.HII_Rinner_per_Rs
                        nh = cfg.par.HII_nh       
                        escape_fraction  = cfg.par.HII_escape_fraction
                    
                    elif pagb:
                        id_val = 1
                        Rinner_per_Rs = cfg.par.PAGB_Rinner_per_Rs
                        nh = cfg.par.PAGB_nh    
                        escape_fraction  = cfg.par.PAGB_escape_fraction

                    if cfg.par.HII_alpha_enhance: #Setting Zstar based on Fe/H 
                        Fe = stars_list[i].all_metals[-1]

                        # Gizmo metallicity structure, photospheric abundances from Asplund et al. 2009:
                        # Photospheric mass fraction of H = 0.7381
                        # Photospheric mass fraction of Fe = 1.31e-3

                        # Converting from mass fraction to atomic fraction of Fe
                        # Taking atmoic mass of H = 1.008u
                        # Taking atmoic mass of Fe = 55.845u
                        FeH = (Fe/0.7381)*(1.008/55.845)

                        # Solar atomic fraction of Fe. Calculated by substituting Fe = 1.31e-3 in the previous equation
                        FeH_sol = 3.22580645e-5

                        Logzsol = np.log10(FeH/FeH_sol)
                        
                        sp1 = fsps.StellarPopulation(zcontinuous=1)
                        sp1.params["tage"] = stars_list[i].age
                        sp1.params["imf_type"] = cfg.par.imf_type
                        sp1.params["imf1"] = cfg.par.imf1
                        sp1.params["imf2"] = cfg.par.imf2
                        sp1.params["imf3"] = cfg.par.imf3
                        sp1.params["pagb"] = cfg.par.pagb
                        sp1.params["sfh"] = 0
                        sp1.params["zmet"] = stars_list[i].fsps_zmet
                        sp1.params["add_neb_emission"] = False
                        sp1.params["add_agb_dust_model"] = cfg.par.add_agb_dust_model
                        sp1.params["logzsol"] = Logzsol

                        if cfg.par.CF_on == True:
                            sp1.params["dust_type"] = 0
                            sp1.params["dust1"] = 1
                            sp1.params["dust2"] = 0
                            sp1.params["dust_tesc"] = tesc_age

                        spec = sp1.get_spectrum(tage=stars_list[i].age)

                    else:
                        spec = sp.get_spectrum(tage=stars_list[i].age,zmet=stars_list[i].fsps_zmet)

                    alpha = 2.5e-13 # Recombination Rate (assuming T = 10^4 K)

                    if cfg.par.FORCE_gas_logu[id_val]:
                        LogU = cfg.par.gas_logu[id_val]
                        LogQ = np.log10((10 ** (3*LogU))*(36*np.pi*(constants.c.cgs.value**3))/((alpha**2)*nh))
                        Rs = ((3*(10 ** LogQ))/(4*np.pi*(nh**2)*alpha))**(1./3.)
                    
                    elif cfg.par.FORCE_logq[id_val]:
                        LogQ = cfg.par.source_logq[id_val]
                        Rs = ((3*(10 ** LogQ))/(4*np.pi*(nh**2)*alpha))**(1./3.)
                        LogU = np.log10((10**LogQ)/(4*np.pi*Rs*Rs*nh*constants.c.cgs.value))

                    else:
                        LogQ = calc_LogQ(1.e8*constants.c.cgs.value/spec[0], spec[1]*constants.L_sun.cgs.value
                                , efrac=escape_fraction, mstar=10**cluster_mass[j])   
                        Rs = ((3*(10 ** LogQ))/(4*np.pi*(nh**2)*alpha))**(1./3.)
                        LogU = np.log10((10**LogQ)/(4*np.pi*Rs*Rs*nh*constants.c.cgs.value))+cfg.par.gas_logu_init[id_val]
                        LogQ = np.log10((10 ** (3*LogU))*(36*np.pi*(constants.c.cgs.value**3))/((alpha**2)*nh))
                        Rs = ((3*(10 ** LogQ))/(4*np.pi*(nh**2)*alpha))**(1./3.)

                    if cfg.par.FORCE_inner_radius[id_val]:
                        Rin = cfg.par.inner_radius[id_val]

                    else:
                        Rin = Rinner_per_Rs*Rs

                    if cfg.par.FORCE_gas_logz[id_val]:
                        LogZ = cfg.par.gas_logz[id_val]
                    
                    else:
                        LogZ = np.log10(stars_list[i].metals/cfg.par.solar)


                    if neb_file_output:
                        if cfg.par.use_cloudy_tables:
                            Rin = 1.e19   # Rinner is fixed at 1.e19 cm for lookup tables
                        
                        if cfg.par.FORCE_inner_radius[id_val]:
                            Rin = cfg.par.inner_radius[id_val]
                        
                        LogU = np.log10((10**LogQ)/(4*np.pi*Rin*Rin*nh*constants.c.cgs.value))

                        logu_diagnostic(LogQ, LogU, LogZ, Rs, 10**cluster_mass[j], num_HII_clusters, stars_list[i].age, append=True)
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
                            LogQ_1 = calc_LogQ(1.e8 * constants.c.cgs.value / spec[0], spec[1] * constants.L_sun.cgs.value,
                                    efrac=escape_fraction) 
                            #LogQ_1 = LogQ_1 + cfg.par.gas_logu_init[id_val]
                               
                            spec_neb, wave_line, line_lum = get_nebular(spec[0], spec[1], nh, stars_list[i].all_metals, logq = LogQ, radius = Rin, 
                                                    logu = LogU, logz = LogZ, logq_1 = LogQ_1, Dust=False, abund=cfg.par.neb_abund[id_val], 
                                                    clean_up = cfg.par.cloudy_cleanup, index=id_val)
                        except ValueError as err:
                            # If the CLOUDY run crashes we switch to using lookup tables for young stars but throw an error for post-AGB stars.
                            if  young_star:
                                print ("WARNING: Switching to using lookup tables pre-packed with FSPS to calculate nebular emission for this particle.") 
                                print ("WARNING: The emission line fluxes repoted may not be accurate if the particle lies outside the range of the lookup table paramters.")
                                lam_neb, spec_neb = sp.get_spectrum(tage=stars_list[i].age, zmet=stars_list[i].fsps_zmet)
                                line_lum = sp.emline_luminosity
                                wave_line = sp.emline_wavelengths
                            else:
                                print ("ERROR: Can't switch to using lookup tables.")
                                print ("ERROR: Please check the CLOUDY output file to figure out why the run was unsuccessful" )
                                raise ValueError('CLOUDY run was unsucessful')
                
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
                dump_emlines(wave_line, line_em, id_val)

        stellar_nu[:] = 1.e8*constants.c.cgs.value/spec[0]
        stellar_fnu[i,:] = f
    return stellar_fnu, mfrac


def get_gas_metals(ngas):
    # This function outputs the metallicity (total as well as all the 10 elements tracked by the simulation)
    # for all the gas particles
    reg = reg_gl
    el = ['He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe']
    metals = np.zeros((ngas,11))-10.0
    try:
        for i in range(11):
            if i == 0:
                el_str = ""
            else:
                el_str = "_"+el[i-1]
            metals[:, i] = reg["gas","metals"+el_str].value
    except:
        metals[:,0] = reg["gas","metals"].value
    
    return metals


def get_nearest_gas_metals(all_gas_coordinates, particle_coordinates):

    all_gas_metals = get_gas_metals(len(all_gas_coordinates))
    
    # Getting N nearest gas particles to the AGN where N is defined by cfg.par.AGN_num_gas
    nearest_gas_dist, nearest_gas_id = get_nearest(all_gas_coordinates, particle_coordinates, num=cfg.par.AGN_num_gas)
    
    metals_avg = []
    
    # We take the distance weighted avearge of the metallicity of nearest N gas particles.
    # This is used as input to the CLOUDY model.
    for q in range(11):
        nearest_gas_metals = np.array(all_gas_metals[:,q][nearest_gas_id])
        metals_avg.append(np.sum(nearest_gas_metals*nearest_gas_dist)/np.sum(nearest_gas_dist))

    return metals_avg
                                                                                                                

def get_agn_seds(agn_ids, reg):
  
    print ('Starting AGN SED generation')
    
    if cfg.par.dump_emlines:
        dump_emlines(None, None, 2, add_header=True)

    t1 = datetime.now()
    nprocesses = np.min([cfg.par.n_processes,len(agn_ids)])
    p = Pool(processes = nprocesses)

    # Pre-calculating average metallicity and fluxes for all the BH particles
    nu = []
    fnu_in = []
    metals_avg = []
    for agn_id in agn_ids:
        fnu_in.append(reg["bh","sed"][agn_id,:].in_units("Lsun").value/reg["bhnu"].in_units("Hz").value)
        nu.append(reg["bh","nu"].in_units("Hz").value)
        all_gas_coordinates = reg["gas","coordinates"].in_units('kpc').value
        agn_coordinates = reg["bh","coordinates"][agn_id].in_units('kpc').value
        metals_avg.append(get_nearest_gas_metals(all_gas_coordinates, agn_coordinates))

    z = zip(agn_ids, nu, fnu_in, metals_avg)
    fnu_out = p.starmap(agn_sed, z)
    fnu_out = np.atleast_2d(fnu_out)
    p.close()
    p.terminate()
    p.join()
    t2 = datetime.now()

    print ('Execution time for AGN SED generation = '+str(t2-t1))
    
    return fnu_out


def agn_sed(agn_id, nu, fnu, metals_avg):
    agn_id = int(agn_id)

    # Hopkins model returns the nu and fnu in reversed order. 
    # Thus, we have to reverse it back to make it compatible.
    if cfg.par.BH_model != 'Nenkova': 
        nu = nu[::-1]
        fnu = fnu[::-1]

    if cfg.par.add_neb_emission and cfg.par.add_AGN_neb:
        
        id_val = 2
            
        tot_metals = metals_avg[0]
        metals = metals_avg
        
        if cfg.par.FORCE_gas_logz[id_val]:
            LogZ = gas_logz[id_val]
        
        else:
            LogZ = tot_metals/cfg.par.solar

        Rin = cfg.par.inner_radius[id_val]

        if cfg.par.FORCE_gas_logu[id_val]:
            LogU = cfg.par.gas_logu[id_val]
            LogQ = np.log10((10**LogU)*(4*np.pi*Rin*Rin*cfg.par.AGN_nh*constants.c.cgs.value))
        
        elif cfg.par.FORCE_logq[id_val]:
            LogQ = cfg.par.source_logq[id_val]
            LogU = np.log10((10**LogQ)/(4*np.pi*Rin*Rin*cfg.par.AGN_nh*constants.c.cgs.value))

        else:
            LogQ = calc_LogQ(nu, fnu*constants.L_sun.cgs.value)
            LogU = np.log10((10**LogQ)/(4*np.pi*Rin*Rin*cfg.par.AGN_nh*constants.c.cgs.value))+cfg.par.gas_logu_init[id_val]
            LogQ = np.log10((10**LogU)*(4*np.pi*Rin*Rin*cfg.par.AGN_nh*constants.c.cgs.value))

        spec, wave_line, line_lum = get_nebular(1.e8 * constants.c.cgs.value / nu, fnu, cfg.par.AGN_nh, metals, logq = LogQ, radius = Rin, logu = LogU, 
                                        logz = LogZ, logq_1 = LogQ, Dust=False, abund=cfg.par.neb_abund[id_val],clean_up = cfg.par.cloudy_cleanup, index=id_val)

        if cfg.par.dump_emlines:
        # The stellar population returns the calculation in units of Lsun
            line_em = line_lum * 3.839e33  # Units: ergs/s
            # The last column in dump_emlines is reserved for age of the star particle. For AGN we just set it to -1 as a place holder.
            line_em = np.append(line_em, -1.0)
            dump_emlines(wave_line, line_em, id_val)

    else:
        spec = fnu

    return spec


def get_dig_seds(factors, cell_widths, metals):

    print ('Starting DIG SED generation')
    
    t1 = datetime.now()
    nprocesses = np.min([cfg.par.n_processes,len(cell_widths)])
    p = Pool(processes = nprocesses)

    z = zip(factors, cell_widths, metals)
    fnu_out = p.starmap(dig_sed, z)
    fnu_out = np.atleast_2d(fnu_out)
    
    p.close()
    p.terminate()
    p.join()
    t2 = datetime.now()

    print ('Execution time for DIG SED generation = '+str(t2-t1))

    return fnu_out


def dig_sed(factor, cell_width, metal):    
    id_val = 3
    
    dat = np.load(cfg.par.pd_source_dir + "/powderday/nebular_emission/data/black_1987.npz")
    spec_lam = dat["lam"]
    sspi = dat["sed"]*(cell_width**2)*factor # Lsun/Hz
    
    spec, wave_line, line_lum = get_nebular(spec_lam, sspi, cfg.par.DIG_nh, metal, Factor=factor, Cell_width=cell_width, Dust=False,
                                            abund=cfg.par.neb_abund[id_val], clean_up = cfg.par.cloudy_cleanup, index=id_val)

    if cfg.par.dump_emlines:
        # The stellar population returns the calculation in units of Lsun
        line_em = line_lum * 3.839e33  # Units: ergs/s
        # The last column in dump_emlines is reserved for age of the star particle. For DIG we just set it to -1 as a place holder.
        line_em = np.append(line_em, -1.0)
        dump_emlines(wave_line, line_em, id_val)
        
    return spec


def fsps_metallicity_interpolate(metals):

    # takes a list of metallicities for star particles, and returns a
    # list of interpolated metallicities
    
    fsps_metals = np.loadtxt(cfg.par.metallicity_legend)
    nstars = len(metals)
    
    zmet = []

    for i in range(nstars):
        zmet.append(find_nearest_zmet(fsps_metals,metals[i]))
    
    return zmet
    
def find_nearest_zmet(array,value):
    # this is modified from the normal find_nearest in that it forces
    # the output to be 1 index higher than the true value since the
    # minimum zmet value fsps will take is 1 (not 0)

    idx = (np.abs(array-value)).argmin()
    
    return idx+1      
