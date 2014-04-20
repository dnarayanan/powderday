#source_creation.py
#purpose: to add new star sources, disks and bulge stars
import parameters as par
import pdb
import numpy as np
import constants as const

def add_newstars(df_nu,stellar_nu,stellar_fnu,disk_fnu,bulge_fnu,stars_list,diskstars_list,bulgestars_list,m):
    
    if par.COSMOFLAG == False:
        
        nstars = len(stars_list)
        nstars_disk = len(diskstars_list)
        nstars_bulge = len(bulgestars_list)

        print 'adding new stars to the grid'
  
        for i in range(nstars):
            
            nu = stellar_nu[:]
            fnu = stellar_fnu[i,:]
    

            nu_inrange = np.logical_and(nu >= min(df_nu),nu <= max(df_nu))
            nu_inrange = np.where(nu_inrange == True)[0]
            nu = nu[nu_inrange]

            #reverse the arrays for hyperion
            nu = nu[::-1]
            fnu = fnu[::-1]

            fnu = fnu[nu_inrange]

            lum = np.absolute(np.trapz(fnu,x=nu))*stars_list[i].mass/const.msun 
            lum *= const.lsun #to get in cgs

            #add new stars

            m.add_spherical_source(luminosity = lum,radius = 10.*const.rsun,spectrum = (nu,fnu),
                                   position = stars_list[i].positions)



        print 'Non-Cosmological Simulation: Adding Disk and Bulge Stars:'
            
        print 'adding disk stars to the grid: adding as a point source collection'   
        disksource = m.add_point_source_collection()
       
    
        fnu = disk_fnu[:]
        fnu = fnu[::-1]
        fnu = fnu[nu_inrange]

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
        fnu = fnu[::-1]
        fnu = fnu[nu_inrange]

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
    
        m.set_sample_sources_evenly(True)
