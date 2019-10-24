from __future__ import print_function
import numpy as np
import yt
import matplotlib
matplotlib.use('Agg')
import powderday.config as cfg
from astropy import constants
import astropy.units as u
from hyperion.model import ModelOutput


def proj_plots(pf):
    print ('\n[analytics/proj_plots] Saving Diagnostic Projection Plots \n')
    p = yt.ProjectionPlot(pf,"x",("gas","density"),width=(cfg.par.zoom_box_len,'kpc'))
    p.save(cfg.model.PD_output_dir+'/proj_plot_x.png')
    p = yt.ProjectionPlot(pf,"y",("gas","density"),width=(cfg.par.zoom_box_len,'kpc'))
    p.save(cfg.model.PD_output_dir+'/proj_plot_y.png')
    p = yt.ProjectionPlot(pf,"z",("gas","density"),width=(cfg.par.zoom_box_len,'kpc'))
    p.save(cfg.model.PD_output_dir+'/proj_plot_z.png')
    
 
    return None


def stellar_sed_write(m):
    
    totallum = 0
    nsources = len(m.sources)
    
    for i in range(nsources):
        tempnu = m.sources[i].spectrum["nu"]
        tempfnu = m.sources[i].spectrum["fnu"]

        if i == 0: fnu = np.zeros(len(tempnu))
        

        #now we need to scale this because the spectrum is just in
        #terms of an SSP, and we need to scale by the total luminosity
        #that wen t into the model (i.e. by the actual stellar mass
        #used in powderday).
        ssp_lum = np.absolute(np.trapz(tempnu,tempfnu))*constants.L_sun.cgs
        lum_scale = np.sum(m.sources[i].luminosity)/ssp_lum #we have to do np.sum in case the sources were in a collection
        tempfnu *= lum_scale.value
        



        for i in range(len(fnu)):
            fnu[i] += tempfnu[i]

    nu = tempnu 

    #go ahead and calculate lambda and flambda
    lam = (constants.c/(nu*u.Hz)).cgs
    flam = nu*u.Hz*(fnu*u.Lsun/u.Hz)/lam
    flam = flam.to(u.Lsun/u.micron)
    lam = lam.to(u.micron)


    README = "Note: nu is in Hz, and fnu is in Lsun/Hz; lam is in micron and flam is in Lsun/micron"
    #saving: nu is in Hz and fnu is in Lsun/Hz
    outfile = cfg.model.PD_output_dir+"stellar_seds."+cfg.model.snapnum_str+".npz"



    np.savez(outfile,nu=nu,fnu=fnu,lam = lam.value, flam = flam.value, README=README)
   

def dump_cell_info(refined,fc1,fw1,xmin,xmax,ymin,ymax,zmin,zmax):
    outfile = cfg.model.PD_output_dir+"cell_info."+cfg.model.snapnum_str+".npz"
    np.savez(outfile,refined=refined,fc1=fc1,fw1=fw1,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,zmin=zmin,zmax=zmax)

def dump_data(pf,model):
    ad = pf.all_data()
    particle_fh2 = ad["gasfh2"]
    particle_fh1 = np.ones(len(particle_fh2))-particle_fh2
    particle_gas_mass = ad["gasmasses"]
    particle_star_mass = ad["starmasses"]
    particle_star_metallicity = ad["starmetals"]
    particle_stellar_formation_time = ad["starformationtime"]
    particle_sfr = ad['gassfr'].in_units('Msun/yr')

    #these are in try/excepts in case we're not dealing with gadget and yt 3.x
    try: grid_gas_mass = ad["gassmoothedmasses"]
    except: grid_gas_mass = -1
    try: grid_gas_metallicity = ad["gassmoothedmetals"]
    except: grid_gas_metallicity = -1
    try: grid_star_mass = ad["starsmoothedmasses"]
    except: grid_star_mass = -1

    #get tdust
    m = ModelOutput(model.outputfile+'.sed')
    oct = m.get_quantities()
    tdust_pf = oct.to_yt()
    tdust_ad = tdust_pf.all_data()
    tdust = tdust_ad[ ('gas', 'temperature')]


    try: outfile = cfg.model.PD_output_dir+"grid_physical_properties."+cfg.model.snapnum_str+'_galaxy'+cfg.model.galaxy_num_str+".npz"
    except:
        outfile = cfg.model.PD_output_dir+"grid_physical_properties."+cfg.model.snapnum_str+".npz"

    np.savez(outfile,particle_fh2=particle_fh2,particle_fh1 = particle_fh1,particle_gas_mass = particle_gas_mass,particle_star_mass = particle_star_mass,particle_star_metallicity = particle_star_metallicity,particle_stellar_formation_time = particle_stellar_formation_time,grid_gas_metallicity = grid_gas_metallicity,grid_gas_mass = grid_gas_mass,grid_star_mass = grid_star_mass,particle_sfr = particle_sfr,tdust = tdust)


def SKIRT_data_dump(pf,ad,m,stars_list,hsml_in_pc):

    #create stars file.  this assumes the 'extragalactic [length in pc, distance in Mpc]' units for SKIRT

    spos_x = ad["starcoordinates"][:,0].in_units('pc').value
    spos_y = ad["starcoordinates"][:,1].in_units('pc').value
    spos_z = ad["starcoordinates"][:,2].in_units('pc').value
    smasses = ad["starmasses"].in_units('Msun').value

    try:
        disk_x = ad["diskstarcoordinates"][:,0].in_units('pc').value
        disk_y = ad["diskstarcoordinates"][:,1].in_units('pc').value
        disk_z = ad["diskstarcoordinates"][:,2].in_units('pc').value
        diskmasses = ad["diskstarmasses"].in_units('Msun').value
    except:
        disk_x, disk_y, disk_z, diskmasses = (np.array([]),)*4

    try:
        bulge_x = ad["bulgestarcoordinates"][:,0].in_units('pc').value
        bulge_y = ad["bulgestarcoordinates"][:,1].in_units('pc').value
        bulge_z = ad["bulgestarcoordinates"][:,2].in_units('pc').value
        bulgemasses = ad["bulgestarmasses"].in_units('Msun').value
    except:
        bulge_x, bulge_y, bulge_z, bulgemasses = (np.array([]),)*4

    spos_x = np.concatenate((spos_x, disk_x, bulge_x))
    spos_y = np.concatenate((spos_y, disk_y, bulge_y))
    spos_z = np.concatenate((spos_z, disk_z, bulge_z))
    smasses = np.concatenate((smasses, diskmasses, bulgemasses))

    fsps_metals = np.loadtxt(cfg.par.metallicity_legend)

    dmet = [fsps_metals[cfg.par.disk_stars_metals]]*len(diskmasses)
    dage = [(cfg.par.disk_stars_age*u.Gyr).to(u.yr).value]*len(diskmasses)

    bmet = [fsps_metals[cfg.par.bulge_stars_metals]]*len(bulgemasses)
    bage = [(cfg.par.bulge_stars_age*u.Gyr).to(u.yr).value]*len(bulgemasses)

    #ages and metallicities need to come from the stars list in case
    #we do something in parameters master to change the values
    smetallicity = [stars.metals for stars in stars_list] + dmet + bmet
    sage = [(stars.age*u.Gyr).to(u.yr).value for stars in stars_list] + dage + bage
    shsml = np.repeat(hsml_in_pc,len(sage))

    try: outfile = cfg.model.PD_output_dir+"SKIRT."+cfg.model.snapnum_str+'_galaxy'+cfg.model.galaxy_num_str+".stars.particles.txt"
    except: outfile = cfg.model.PD_output_dir+"SKIRT."+cfg.model.snapnum_str+".stars.particles.txt"
    np.savetxt(outfile, np.column_stack((spos_x,spos_y,spos_z,shsml,smasses,smetallicity,sage)))
    
    #create the gas file.  this assumes the 'extragalactic [length in pc, distance in Mpc]' units for SKIRT
    gpos_x = ad["gascoordinates"][:,0].in_units('pc').value
    gpos_y = ad["gascoordinates"][:,1].in_units('pc').value
    gpos_z = ad["gascoordinates"][:,2].in_units('pc').value
    gmass = ad["gasmasses"].in_units('Msun').value
    ghsml = np.repeat(hsml_in_pc,len(gpos_x))
    gmetallicity = ad["gasmetals"].value

    try: outfile = cfg.model.PD_output_dir+"SKIRT."+cfg.model.snapnum_str+'_galaxy'+cfg.model.galaxy_num_str+".gas.particles.txt"
    except: outfile = cfg.model.PD_output_dir+"SKIRT."+cfg.model.snapnum_str+".gas.particles.txt"
    np.savetxt(outfile, np.column_stack((gpos_x,gpos_y,gpos_z,ghsml,gmass,gmetallicity)))
    
# Saves logU, Q and other related parameters in a file (seperate file is created for each galaxy)
def logu_diagnostic(logU, Q, mstar, age, zmet,append = True):
    if append == False:
        try: outfile = cfg.model.PD_output_dir + "nebular_properties_galaxy" + cfg.model.galaxy_num_str + ".txt"
        except: outfile = cfg.model.PD_output_dir + "nebular_properties_galaxy.txt"
        f = open(outfile, 'w+')
        f.close()
    else:
        try:outfile = cfg.model.PD_output_dir + "nebular_properties_galaxy" + cfg.model.galaxy_num_str + ".txt"
        except: outfile = cfg.model.PD_output_dir + "nebular_properties_galaxy.txt"
        f = open(outfile, 'a+')
        f.write(str(logU) + "\t" + str(Q) + "\t" + str(mstar) + "\t"+ str(age) + "\t" + str(zmet) + "\n")
        f.close()

    # Dumps AGN SEDs
def dump_AGN_SEDs(nu,fnu,luminosity):
    savefile = cfg.model.PD_output_dir+"/bh_sed.npz"
    np.savez(savefile,nu = nu,fnu = fnu, luminosity = luminosity)
                      
