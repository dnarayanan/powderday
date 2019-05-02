from __future__ import print_function
import numpy as np
import yt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import config as cfg
import pdb
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



def mass_weighted_distribution(array,masses,fileout='mwd.png',nbins=100):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.hist(array,bins=nbins,weights=masses,log=True,normed=True)
    
    fig.savefig(fileout)


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


def dust_histograms(refined,dust_smoothed_dtm,dust_smoothed_remy_ruyer):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    wFalse = np.where(np.array(refined) == False)[0]
    d_dtm = dust_smoothed_dtm[wFalse]
    d_dtm = d_dtm[d_dtm>0]
    d_rr = dust_smoothed_remy_ruyer[wFalse]
    d_rr = d_rr[d_rr>0]
    histvals,binvals,patches = plt.hist(np.log10(d_dtm),bins=100)
    histvals,binvals,patches = plt.hist(np.log10(d_rr),bins=100)
    ax.set_xlabel('dust density')
    ax.set_ylabel('N')
    fig.savefig(cfg.model.PD_output_dir+'dust_density.'+cfg.model.snapnum_str+'.png',dpi=300)

def SKIRT_data_dump(pf,ad,m,stars_list,hsml_in_pc):

    #create stars file.  this assumes the 'extragalactic [length in pc, distance in Mpc]' units for SKIRT

    #ages and metallicities need to come from the stars list in case
    #we do something in parameters master to change the values
    smetallicity = [stars.metals for stars in stars_list]
    sage = [(stars.age*u.Gyr).to(u.yr).value for stars in stars_list] #to get in yr
    shsml = np.repeat(hsml_in_pc,len(sage))

    spos_x = ad["starcoordinates"][:,0].in_units('pc').value
    spos_y = ad["starcoordinates"][:,1].in_units('pc').value
    spos_z = ad["starcoordinates"][:,2].in_units('pc').value
    smasses = ad["starmasses"].in_units('Msun').value

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
    

    
                      
