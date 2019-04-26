from __future__ import print_function
import numpy as np
import yt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import config as cfg
import pdb,ipdb
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

def dump_data(ad,model):

    try: particle_fh2 = ad[('PartType0', 'FractionH2')]
    except: particle_fh2 = np.zeros(len(ad[('PartType0','Masses')]))
    particle_fh1 = np.ones(len(particle_fh2))-particle_fh2
    particle_gas_mass = ad[('PartType0','Masses')].in_units('Msun')
    particle_star_mass = ad[('PartType4','Masses')].in_units('Msun')
    particle_star_metallicity = ad[('PartType4','Metallicity_00')]
    particle_stellar_formation_time = ad[('PartType4', 'StellarFormationTime')]
    grid_gas_mass = ad[('deposit', 'PartType0_mass')].in_units('Msun')
    grid_gas_metallicity = ad[('deposit', 'PartType0_smoothed_metallicity')]
    grid_star_mass = ad[ ('deposit', 'PartType4_mass')].in_units('Msun')
    grid_star_metallicity = ad[('PartType4','Metallicity_00')]

    #get tdust

    m = ModelOutput(model.outputfile+'.sed')
    oct = m.get_quantities()
    tdust_pf = oct.to_yt()
    tdust_ad = tdust_pf.all_data()
    tdust = tdust_ad[ ('gas', 'temperature')]


    try: outfile = cfg.model.PD_output_dir+"grid_physical_properties."+cfg.model.snapnum_str+'_galaxy'+cfg.model.galaxy_num_str+".npz"
    except:
        outfile = cfg.model.PD_output_dir+"grid_physical_properties."+cfg.model.snapnum_str+".npz"

    np.savez(outfile,particle_fh2=particle_fh2,particle_fh1 = particle_fh1,particle_gas_mass = particle_gas_mass,particle_star_mass = particle_star_mass,particle_star_metallicity = particle_star_metallicity,particle_stellar_formation_time = particle_stellar_formation_time,grid_gas_metallicity = grid_gas_metallicity,grid_gas_mass = grid_gas_mass,grid_star_mass = grid_star_mass,grid_star_metallicity = grid_star_metallicity,tdust=tdust)


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

def SKIRT_data_dump(pf,m,stars_list,hsml_in_pc):

    #create stars file.  this assumes SI units for SKIRT
    stars_pos_x = [pf.quan(stars.positions[0],'code_length').in_units('pc').value for stars in stars_list]
    stars_pos_y = [pf.quan(stars.positions[1],'code_length').in_units('pc').value for stars in stars_list]
    stars_pos_z = [pf.quan(stars.positions[2],'code_length').in_units('pc').value for stars in stars_list]

    metallicity = [stars.metals for stars in stars_list]

    masses = [(stars.mass*u.g).to(u.Msun).value for stars in stars_list]
    age = [(stars.age*u.Gyr).to(u.yr).value for stars in stars_list] #to get in yr
    hsml = np.repeat(hsml_in_pc,len(age))



    try: outfile = cfg.model.PD_output_dir+"SKIRT."+cfg.model.snapnum_str+'_galaxy'+cfg.model.galaxy_num_str+".particles.txt"
    except: outfile = cfg.model.PD_output_dir+"SKIRT."+cfg.model.snapnum_str+".particles.txt"
    np.savetxt(outfile, np.column_stack((stars_pos_x,stars_pos_y,stars_pos_z,hsml,masses,metallicity,age)))
    
    #create wavelength grid
    '''
    nu = m.sources[0].spectrum["nu"]*u.Hz
    lam = ((constants.c)/nu).to(u.micron).value
    try: outfile = cfg.model.PD_output_dir+"SKIRT."+cfg.model.snapnum_str+'_galaxy'+cfg.model.galaxy_num_str+".wavelengths.txt"
    except: outfile = cfg.model.PD_output_dir+"SKIRT."+cfg.model.snapnum_str+".wavelengths.txt"
    np.savetxt(outfile,lam[::-1])
    '''
