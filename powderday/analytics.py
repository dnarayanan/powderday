from __future__ import print_function
import numpy as np
import yt
import matplotlib
matplotlib.use('Agg')
import powderday.config as cfg
from astropy import constants
import astropy.units as u
from hyperion.model import ModelOutput
import os,pdb

def proj_plots(ds):
    print ('\n[analytics/proj_plots] Saving Diagnostic Projection Plots \n')
    p = yt.ProjectionPlot(ds,"x",("gas","density"),width=(cfg.par.zoom_box_len,'kpc'))
    p.save(cfg.model.PD_output_dir+'/proj_plot_x.png')
    p = yt.ProjectionPlot(ds,"y",("gas","density"),width=(cfg.par.zoom_box_len,'kpc'))
    p.save(cfg.model.PD_output_dir+'/proj_plot_y.png')
    p = yt.ProjectionPlot(ds,"z",("gas","density"),width=(cfg.par.zoom_box_len,'kpc'))
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

    
    try: outfile = cfg.model.PD_output_dir+"/stellar_seds."+cfg.model.snapnum_str+'_galaxy'+cfg.model.galaxy_num_str+".npz"
    except:
        outfile = cfg.model.PD_output_dir+"/stellar_seds."+cfg.model.snapnum_str+".npz"

    np.savez(outfile,nu=nu,fnu=fnu,lam = lam.value, flam = flam.value, README=README)
   

def dump_cell_info(refined,fc1,fw1,xmin,xmax,ymin,ymax,zmin,zmax):
    outfile = cfg.model.PD_output_dir+"cell_info."+cfg.model.snapnum_str+"_"+cfg.model.galaxy_num_str+".npz"
    np.savez(outfile,refined=refined,fc1=fc1,fw1=fw1,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,zmin=zmin,zmax=zmax)

def dump_data(reg,model):
    
    particle_fh2 = reg["gas","fh2"]
    particle_fh1 = np.ones(len(particle_fh2))-particle_fh2
    particle_gas_mass = reg["gas","masses"]
    particle_star_mass = reg["star","masses"]
    particle_star_metallicity = reg["star","metals"]
    #particle_stellar_formation_time = reg["starformationtime"]
    particle_stellar_formation_time = reg["stellar","ages"]
    particle_sfr = reg['gas','sfr'].in_units('Msun/yr')
    

    try: cell_size = reg.parameters['cell_size']
    except: cell_size = -1

    #save dustmasses.  for particle type codes where the particles are
    #projected onto an octree (gizmo/gadget) this will be the
    #smoothedmasses field; else, it will just be the ('dust','mass')
    #tuple which is the mesh values 
    try: grid_dustmass = reg['dust','smoothedmasses'].in_units('Msun')
    except: grid_dustmass = reg["dust","mass"].in_units('Msun')  #for arepo/AMR codes,
                                                 #this is the grid
                                                 #values (and also
                                                 #parttype0 for
                                                 #arepo).  for
                                                 #gizmo/gadget, this
                                                 #is particle
                                                 #information which is
                                                 #why it has to be in
                                                 #the except
                                                 #statement, and not
                                                 #the try statemetn in logical flow.

    #if we have separate particle dust masses, save them. 
    try: particle_dustmass = reg['particle_dust','mass'].in_units('Msun')
    except: particle_dustmass = -1


    #these are in try/excepts in case we're not dealing with gadget and yt 3.x
    
    if reg.parameters['dataset_type'] in ['gadget_hdf5','tipsy']:
        try: 
            grid_gas_mass = reg["gas","smoothedmasses"]
        except: grid_gas_mass = -1
    else:
        try:
            grid_gas_mass = reg["gas","masses"]
        except: grid_gas_mass = -1
        
    #for sph we have a different nomenclature thanks to the smoothing onto an octree
    if reg.parameters['dataset_type'] in ['gadget_hdf5','tipsy']:

        try: 
            grid_gas_metallicity = []
            grid_gas_metallicity.append(reg["gas","smoothedmetals"].value)
            abund_el = ['He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe']
            for i in abund_el:
                grid_gas_metallicity.append(reg["gas","smoothedmetals_"+str(i)].value)
                
        except: grid_gas_metallicity = reg["gas","smoothedmetals"].value
    
    else:
        try:
            grid_gas_metallicity = []
            grid_gas_metallicity.append(reg["gas","metals"].value)
            abund_el = ['He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe']
            for i in abund_el:
                grid_gas_metallicity.append(reg["gas","metals_"+str(i)].value)

        except: grid_gas_metallicity = reg["gas","metals"].value

    try: grid_star_mass = reg["star","smoothedmasses"]
    except: grid_star_mass = -1
    
    #PAH information
    try: grid_PAH_luminosity = reg.parameters['grid_PAH_luminosity']
    except: grid_PAH_luminosity = -1
    try: grid_neutral_PAH_luminosity = reg.parameters['grid_neutral_PAH_luminosity']
    except: grid_neutral_PAH_luminosity = -1
    try: grid_ion_PAH_luminosity = reg.parameters['grid_ion_PAH_luminosity']
    except: grid_ion_PAH_luminosity = -1


    try: PAH_lam = reg.parameters['PAH_lam']
    except: PAH_lam = -1
    try: total_PAH_luminosity = reg.parameters['total_PAH_luminosity']
    except: total_PAH_luminosity = -1
    try: total_neutral_PAH_luminosity = reg.parameters['total_neutral_PAH_luminosity']
    except: total_neutral_PAH_luminosity = -1
    try: total_ion_PAH_luminosity = reg.parameters['total_ion_PAH_luminosity']
    except: total_ion_PAH_luminosity = -1


    try: integrated_grid_PAH_luminosity = reg.parameters['integrated_grid_PAH_luminosity']
    except: integrated_grid_PAH_luminosity = -1
    try: integrated_grid_neutral_PAH_luminosity = reg.parameters['integrated_grid_neutral_PAH_luminosity']
    except: integrated_grid_neutral_PAH_luminosity = -1
    try: integrated_grid_ion_PAH_luminosity = reg.parameters['integrated_grid_ion_PAH_luminosity']
    except: integrated_grid_ion_PAH_luminosity = -1


    try: m_pah = reg.parameters['m_pah']
    except: m_pah = -1
    try: q_pah = reg.parameters['q_pah']
    except: q_pah = -1
    try: particle_mass_weighted_gsd = reg.parameters['particle_mass_weighted_gsd']
    except: particle_mass_weighted_gsd = -1
    try: grid_mass_weighted_gsd = reg.parameters['grid_mass_weighted_gsd']
    except: grid_mass_weighted_gsd = -1
    try: simulation_sizes = reg.parameters['simulation_sizes']
    except: simulation_sizes = -1
    try: beta_nnls = reg.parameters['beta_nnls']
    except: beta_nnls = -1

    #get tdust
    #m = ModelOutput(model.outputfile+'.sed')
    #oct = m.get_quantities()
    #tdust_ds = oct.to_yt()
    #tdust_ad = tdust_ds.all_data()
    #tdust = tdust_ad[ ('gas', 'temperature')]


    try: outfile = cfg.model.PD_output_dir+"/grid_physical_properties."+cfg.model.snapnum_str+'_galaxy'+cfg.model.galaxy_num_str+".npz"
    except:
        outfile = cfg.model.PD_output_dir+"/grid_physical_properties."+cfg.model.snapnum_str+".npz"

    np.savez(outfile,particle_fh2=particle_fh2,particle_fh1 = particle_fh1,particle_gas_mass = particle_gas_mass,particle_star_mass = particle_star_mass,particle_star_metallicity = particle_star_metallicity,particle_stellar_formation_time = particle_stellar_formation_time,grid_gas_metallicity = grid_gas_metallicity,grid_gas_mass = grid_gas_mass,grid_star_mass = grid_star_mass,particle_sfr = particle_sfr,particle_dustmass = particle_dustmass,grid_dustmass=grid_dustmass,grid_PAH_luminosity = grid_PAH_luminosity,grid_neutral_PAH_luminosity=grid_neutral_PAH_luminosity,grid_ion_PAH_luminosity=grid_ion_PAH_luminosity,PAH_lam=PAH_lam,total_PAH_luminosity = total_PAH_luminosity,total_neutral_PAH_luminosity=total_neutral_PAH_luminosity,total_ion_PAH_luminosity=total_ion_PAH_luminosity,integrated_grid_PAH_luminosity = integrated_grid_PAH_luminosity,integrated_grid_neutral_PAH_luminosity=integrated_grid_neutral_PAH_luminosity,integrated_grid_ion_PAH_luminosity=integrated_grid_ion_PAH_luminosity,m_pah = m_pah, q_pah=q_pah,particle_mass_weighted_gsd = particle_mass_weighted_gsd,grid_mass_weighted_gsd = grid_mass_weighted_gsd,simulation_sizes=simulation_sizes,cell_size=cell_size,beta_nnls=beta_nnls)#,tdust = tdust)


def SKIRT_data_dump(reg,ds,m,stars_list,bulgestars_list,diskstars_list,ds_type,sp,hsml_in_pc = 10):
    
    #the work flow for this function is: for all dataset types, we
    #dump stars in the same manner (since we don't allow for mappings
    #in our skirt dumps, all stars are "old stars").  #for gas,
    #however, we separate based on SPH type outputs or arepo-types
    #since they require different formats.

    #create stars file.  this assumes the 'extragalactic [length in pc, distance in Mpc]' units for SKIRT

    spos_x = ([stars.positions[0] for stars in stars_list]*u.cm).to(u.pc).value
    spos_y = ([stars.positions[1] for stars in stars_list]*u.cm).to(u.pc).value
    spos_z = ([stars.positions[2] for stars in stars_list]*u.cm).to(u.pc).value
    smasses = ([stars.mass for stars in stars_list]*u.g).to(u.Msun).value

    try:
        disk_x = ([diskstars.positions[0] for diskstars in diskstars_list]*u.cm).to(u.pc).value
        disk_y = ([diskstars.positions[1] for diskstars in diskstars_list]*u.cm).to(u.pc).value
        disk_z = ([diskstars.positions[2] for diskstars in diskstars_list]*u.cm).to(u.pc).value
        diskmasses = ([diskstars.mass for diskstars in diskstars_list]*u.g).to(u.Msun).value

    except: 
        disk_x, disk_y, disk_z, diskmasses = (np.array([]),)*4

    try:
        bulge_x = ([bulgestars.positions[0] for bulgestars in bulgestars_list]*u.cm).to(u.pc).value
        bulge_y = ([bulgestars.positions[1] for bulgestars in bulgestars_list]*u.cm).to(u.pc).value
        bulge_z = ([bulgestars.positions[2] for bulgestars in bulgestars_list]*u.cm).to(u.pc).value
        bulgemasses = ([bulgestars.mass for bulgstars in bulgestars_list]*u.g).to(u.Msun).value
    except:
        bulge_x, bulge_y, bulge_z, bulgemasses = (np.array([]),)*4



    spos_x = np.concatenate((spos_x, disk_x, bulge_x))
    spos_y = np.concatenate((spos_y, disk_y, bulge_y))
    spos_z = np.concatenate((spos_z, disk_z, bulge_z))
    smasses = np.concatenate((smasses, diskmasses, bulgemasses))

    fsps_metals = np.array(sp.zlegend)
    
    if ds.cosmological_simulation:
        dmet = [0.]*len(diskmasses)
        dage = [0.]*len(diskmasses)
        bmet = [0.]*len(bulgemasses)
        bage = [0.]*len(bulgemasses)

    else:
        dmet = [fsps_metals[cfg.par.disk_stars_metals]]*len(diskmasses)
        dage = [(cfg.par.disk_stars_age*u.Gyr).to(u.yr).value]*len(diskmasses)
        bmet = [fsps_metals[cfg.par.bulge_stars_metals]]*len(bulgemasses)
        bage = [(cfg.par.bulge_stars_age*u.Gyr).to(u.yr).value]*len(bulgemasses)


    #ages and metallicities need to come from the stars list in case
    #we do something in parameters master to change the values
    smetallicity = [stars.metals for stars in stars_list] + dmet + bmet
    sage = [(stars.age*u.Gyr).to(u.yr).value for stars in stars_list] + dage + bage
    shsml = np.repeat(hsml_in_pc,len(sage))

    #create the gas file for SPH-oids.  this assumes the 'extragalactic [length in pc, distance in Mpc]' units for SKIRT

    gpos_x = reg["gas","coordinates"][:,0].in_units('pc').value
    gpos_y = reg["gas","coordinates"][:,1].in_units('pc').value
    gpos_z = reg["gas","coordinates"][:,2].in_units('pc').value
    gmass = reg["gas","masses"].in_units('Msun').value
    gmetallicity = reg["gas","metals"].value
    grho = (reg["gas","density"]*reg["gas","metals"].value*cfg.par.dusttometals_ratio).in_units('Msun/cm**3').value



    #set the smoothing lengths. see if we have one defined from the front ends.  if not, then we just use a constant value
    try:
        ghsml = reg["gas","smoothinglength"].in_units('pc').value
    except:
        ghsml = np.repeat(hsml_in_pc,len(gpos_x))

    #file I/O
    #stars output
    try: outfile_stars = cfg.model.PD_output_dir+"SKIRT."+cfg.model.snapnum_str+'_galaxy'+cfg.model.galaxy_num_str+".stars.particles.txt"
    except: outfile_stars = cfg.model.PD_output_dir+"SKIRT."+cfg.model.snapnum_str+".stars.particles.txt"
    np.savetxt(outfile_stars, np.column_stack((spos_x,spos_y,spos_z,shsml,smasses,smetallicity,sage)))

    #gas output
    try: outfile_gas = cfg.model.PD_output_dir+"SKIRT."+cfg.model.snapnum_str+'_galaxy'+cfg.model.galaxy_num_str+".gas.particles.txt"
    except: outfile_gas = cfg.model.PD_output_dir+"SKIRT."+cfg.model.snapnum_str+".gas.particles.txt"

    
    if ds_type in ['gadget_hdf5','tipsy']:
        np.savetxt(outfile_gas, np.column_stack((gpos_x,gpos_y,gpos_z,ghsml,gmass,gmetallicity)))
    else:
    #if we ever get the arepo SKIRT ski files working, this is
        #actually the line we need. but since we are currently running
        #SKIRT for arepo in SPH/octree mode, we have to
        # save as though it's an octree..
        #np.savetxt(outfile,np.column_stack((gpos_x,gpos_y,gpos_z,grho)))
        np.savetxt(outfile_gas, np.column_stack((gpos_x,gpos_y,gpos_z,ghsml,gmass,gmetallicity)))


# Saves logU, Q and other related parameters in a file (seperate file is created for each galaxy)
def logu_diagnostic(logQ, LogU, LogZ, Rin, cluster_mass, num_cluster, age, append = True):
    if append == False:
        try: outfile = cfg.model.PD_output_dir + "nebular_properties_galaxy" + cfg.model.galaxy_num_str + ".txt"
        except: outfile = cfg.model.PD_output_dir + "nebular_properties_galaxy.txt"
        f = open(outfile, 'w+')
        f.close()
    else:
        try:outfile = cfg.model.PD_output_dir + "nebular_properties_galaxy" + cfg.model.galaxy_num_str + ".txt"
        except: outfile = cfg.model.PD_output_dir + "nebular_properties_galaxy.txt"
        f = open(outfile, 'a+')
        f.write(str(logQ) + "\t" + str(LogU) + "\t" + str(LogZ) + "\t" + str(Rin) + "\t"+ str(cluster_mass) + "\t" + str(num_cluster) + "\t" + str(age) + "\n")
        f.close()


# Dumps emission lines
def dump_emlines(line_em, append=True):
    if hasattr(cfg.model, 'galaxy_num_str'):
        outfile_lines = cfg.model.PD_output_dir + "emlines.galaxy" + cfg.model.galaxy_num_str + ".txt"
    else:
        outfile_lines = cfg.model.PD_output_dir + "emlines.galaxy.txt"

    if append == False:
        refline_file = cfg.par.pd_source_dir + "/powderday/nebular_emission/data/refLines.dat"
        f = open(outfile_lines,'w')
        wdat = np.genfromtxt(refline_file, delimiter=',')
        wl = np.array([dat[0] for dat in wdat])
        sinds = np.argsort(wl)
        line_wav = wl[sinds]  
        np.savetxt(f,np.expand_dims(line_wav,axis=0))
        f.close()
    else:
        f = open(outfile_lines,'a+')
        np.savetxt(f,np.expand_dims(line_em,axis=0))
        
        f.close()

# Dumps AGN SEDs
def dump_AGN_SEDs(nu,fnu,luminosity):
    
    if hasattr(cfg.model,'galaxy_num_str'):
        outfile_bh = cfg.model.PD_output_dir + "bh_sed." + cfg.model.galaxy_num_str+".npz"
    else:
        outfile_bh = cfg.model.PD_output_dir+"/bh_sed.npz"

    np.savez(outfile_bh,nu = nu,fnu = fnu, luminosity = luminosity)
                      

def dump_NEB_SEDs(nu_arr, fnu_arr, pos_arr, append=True, clean_up=False):
    
    if hasattr(cfg.model,'galaxy_num_str'):
        outfile = cfg.model.PD_output_dir+"/neb_seds_galaxy_"+cfg.model.galaxy_num_str+".npz"
    else:
        outfile = cfg.model.PD_output_dir+"/neb_seds.npz"

    if clean_up:
        os.remove(outfile)
        return
        
    # If append is False then just save an empty npz file.
    if not append:
        np.savez(outfile)
        return

    data = np.load(outfile)
    
    # If the npz file is not empty then the new data is appeneded 
    # to the previous arrays
    if (len(list(data.keys())) == 0):
        fnu_arr_tot = np.atleast_2d(fnu_arr)
        pos_arr_tot = np.atleast_2d(pos_arr)
        
    else:
        fnu_old = np.atleast_2d(data["fnu"])
        pos_old = np.atleast_2d(data["positions"])
        pos_arr = np.atleast_2d(pos_arr)
        fnu_arr = np.atleast_2d(fnu_arr)
        fnu_arr_tot = np.append(fnu_old, fnu_arr, axis=0)
        pos_arr_tot = np.append(pos_old, pos_arr, axis=0)

    np.savez(outfile, nu=nu_arr, fnu=fnu_arr_tot, positions=pos_arr_tot)



'''
#def dump_emline(emline_wavelengths,emline_luminosity,append=True):
def dump_emline(emline_wave,emline_lum,append=True):

    if hasattr(cfg.model, 'galaxy_num_str'):
        outfile_lines = cfg.model.PD_output_dir + "emlines.galaxy" + cfg.model.galaxy_num_str + ".txt"
    else:
        outfile_lines = cfg.model.PD_output_dir + "emlines.galaxy.txt"

    if append == False:
        f = open(outfile_lines,'w')
        #np.expand_dims sets up a dummy dimension so it saves as a row
        np.savetxt(f,np.expand_dims(emline_wave,axis=0))
        f.close()
    else:
        f = open(outfile_lines,'a+')
        np.savetxt(f,emline_lum)#,fmt='%.18g', delimiter=' ', newline=os.linesep)
        f.close()
'''
