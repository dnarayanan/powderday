import numpy as np
import yt
import matplotlib.pyplot as plt
import config as cfg
import pdb,ipdb
from astropy import constants
import astropy.units as u

def proj_plots(pf):
    print '\n[analytics/proj_plots] Saving Diagnostic Projection Plots \n'
    p = yt.ProjectionPlot(pf,"x",("gas","density"),width=(cfg.par.zoom_box_len,'kpc'))
    p.save(cfg.model.PD_output_dir+'/proj_plot_x.png')
    p = yt.ProjectionPlot(pf,"y",("gas","density"),width=(cfg.par.zoom_box_len,'kpc'))
    p.save(cfg.model.PD_output_dir+'/proj_plot_y.png')
    p = yt.ProjectionPlot(pf,"z",("gas","density"),width=(cfg.par.zoom_box_len,'kpc'))
    p.save(cfg.model.PD_output_dir+'/proj_plot_z.png')
    
 
    return None



def mass_weighted_distribution(array,masses,fileout='junk.png',nbins=100):
    
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


    README = "Note: nu is in Hz, and fnu is in Lsun/Hz; lam is in micron and flam is in Lsun/micron; if monochromatic_idx == -1,then there was no monochromatic photons; else, monochromatic_idx refers to the indices in in the nu/fnu array that the monochromatic  photons were run at."
    #saving: nu is in Hz and fnu is in Lsun/Hz
    outfile = cfg.model.PD_output_dir+"stellar_seds."+cfg.model.snapnum_str+".npz"

    #save monochromatic indexes if there are any
    monochromatic_lam = -1 #default
    monochromatic_idx = -1 #default
    if cfg.par.FIX_SED_MONOCHROMATIC_WAVELENGTHS == True:
        monochromatic_nu = m.sources[0].spectrum['nu']*u.Hz
        monochromatic_lam = (constants.c/monochromatic_nu).to(u.micron).value[::-1]
        monochromatic_idx = np.round(np.linspace(np.min(np.where(monochromatic_lam > cfg.par.SED_MONOCHROMATIC_min_lam)[0]),\
                                   np.max(np.where(monochromatic_lam < cfg.par.SED_MONOCHROMATIC_max_lam)[0]),\
                                   cfg.par.SED_MONOCHROMATIC_nlam))

        monochromatic_lam = np.take(monochromatic_lam,list(monochromatic_idx))

    np.savez(outfile,nu=nu,fnu=fnu,lam = lam.value, flam = flam.value, README=README,monochromatic_lam = monochromatic_lam,monochromatic_idx = monochromatic_idx)
   
