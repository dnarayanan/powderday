# coding: utf-8
#Code:  pd_front_end.py

# =========================================================
# IMPORT STATEMENTS
# =========================================================
from __future__ import print_function
from powderday.source_creation import direct_add_stars, add_binned_seds, BH_source_add
from powderday.analytics import stellar_sed_write, dump_data, SKIRT_data_dump, logu_diagnostic
from astropy import constants
import fsps
from powderday.image_processing import add_transmission_filters, convolve
from powderday.m_control_tools import m_control_sph, m_control_enzo, m_control_arepo
import powderday.backwards_compatibility as bc
import powderday.error_handling as eh
import powderday.SED_gen as sg
from powderday.front_ends.front_end_controller import stream
from astropy import units as u
import powderday.config as cfg
import h5py
import matplotlib as mpl
import copy
import numpy as np
import sys
import gc
import random

gc.set_threshold(0)

script, pardir, parfile, modelfile = sys.argv


mpl.use('Agg')


sys.path.insert(0, pardir)
par = __import__(parfile)
model = __import__(modelfile)


cfg.par = par  # re-write cfg.par for all modules that read this in now
cfg.model = model


# =========================================================
# CHECK FOR THE EXISTENCE OF A FEW CRUCIAL FILES FIRST
# =========================================================
eh.file_exist(model.hydro_dir+model.snapshot_name)
eh.file_exist(par.dustdir+par.dustfile)


# =========================================================
# Enforce Backwards Compatibility for Non-Critical Variables
# =========================================================
cfg.par.FORCE_RANDOM_SEED, cfg.par.max_age_direct, cfg.par.BH_SED, cfg.par.IMAGING, cfg.par.SED, cfg.par.IMAGING_TRANSMISSION_FILTER, cfg.par.SED_MONOCHROMATIC, cfg.par.SKIP_RT, cfg.par.FIX_SED_MONOCHROMATIC_WAVELENGTHS, cfg.par.n_MPI_processes, cfg.par.SOURCES_RANDOM_POSITIONS, cfg.par.FORCE_gas_logu, cfg.par.gas_logu, cfg.par.gas_logu_init, cfg.par.gas_logz, cfg.par.FORCE_gas_logz, cfg.par.source_logq, cfg.par.FORCE_logq, cfg.par.FORCE_inner_radius, cfg.par.inner_radius, cfg.par.use_Q, cfg.par.neb_dust, cfg.par.cmdf_min_mass, cfg.par.cmdf_max_mass, cfg.par.cmdf_bins, cfg.par.cmdf_beta, cfg.par.SUBLIMATION, cfg.par.SUBLIMATION_TEMPERATURE, cfg.model.TCMB, cfg.model.THETA, cfg.model.PHI, cfg.par.MANUAL_ORIENTATION, cfg.par.solar, cfg.par.dust_grid_type, cfg.par.BH_model, cfg.par.BH_modelfile, cfg.par.BH_var, cfg.par.FORCE_STELLAR_AGES, cfg.par.FORCE_STELLAR_AGES_VALUE, cfg.par.FORCE_STELLAR_METALLICITIES, cfg.par.FORCE_STELLAR_METALLICITIES_VALUE, cfg.par.HII_T, cfg.par.HII_nh, cfg.par.HII_max_age, cfg.par.HII_escape_fraction, cfg.par.neb_abund, cfg.par.use_cloudy_tables, cfg.par.cloudy_cleanup, cfg.par.neb_file_output, cfg.par.stellar_cluster_mass, cfg.par.filterdir, cfg.par.filterfiles, cfg.par.PAH_frac, cfg.par.FORCE_BINNED = bc.variable_set()

# =========================================================
# GRIDDING
# =========================================================

fname = cfg.model.hydro_dir+cfg.model.snapshot_name
field_add, ds = stream(fname)

# figure out which tributary we're going to

ds_type = ds.dataset_type
# define the options dictionary
options = {'gadget_hdf5': m_control_sph,
           'tipsy': m_control_sph,
           'enzo_packed_3d': m_control_enzo,
           'arepo_hdf5': m_control_arepo}

m_gen = options[ds_type]()
m, xcent, ycent, zcent, dx, dy, dz, reg, ds, boost = m_gen(fname, field_add)


sp = fsps.StellarPopulation()


# Get dust wavelengths. This needs to preceed the generation of sources
# for hyperion since the wavelengths of the SEDs need to fit in the
# dust opacities.

df = h5py.File(cfg.par.dustdir+cfg.par.dustfile, 'r')
o = df['optical_properties']
df_nu = o['nu']
df_chi = o['chi']

df.close()


# add sources to hyperion
stars_list, diskstars_list, bulgestars_list, reg = sg.star_list_gen(boost, dx, dy, dz, reg, ds)
nstars = len(stars_list)

if cfg.par.BH_SED == True:
    BH_source_add(m, reg, df_nu, boost)


# figure out N_METAL_BINS:
fsps_metals = np.loadtxt(cfg.par.metallicity_legend)
N_METAL_BINS = len(fsps_metals)


#initializing the logU file newly
if cfg.par.add_neb_emission: logu_diagnostic(None,None,None,None,None,None,None,None,append=False)

if cfg.par.FORCE_BINNED == False:
    m = direct_add_stars(df_nu, stars_list, diskstars_list, bulgestars_list, ds.cosmological_simulation, m, sp)

# note - the generation of the SEDs is called within
# add_binned_seds itself, unlike add_newstars, which requires
# that sg.allstars_sed_gen() be called first.
m = add_binned_seds(df_nu, stars_list, diskstars_list,bulgestars_list, ds.cosmological_simulation, m, sp)



#set the random seets
if cfg.par.FORCE_RANDOM_SEED == False:
    m.set_seed(random.randrange(0,10000)*-1)
else:
    m.set_seed(cfg.par.seed)

# save SEDs
# stars and black holes can't both be in the sim and write stellar SEDs to a file becuase they have different wavelength sizes
if (par.STELLAR_SED_WRITE == True) and not (par.BH_SED):
    stellar_sed_write(m)


if ds_type in ['gadget_hdf5','tipsy','arepo_hdf5']:
    SKIRT_data_dump(reg, ds, m, stars_list, ds_type, 10)


nstars = len(stars_list)
nstars_disk = len(diskstars_list)
nstars_bulge = len(bulgestars_list)


'''
#EXPERIMENTAL FEATURES
if par.SOURCES_IN_CENTER == True:
    for i in range(nstars):
        stars_list[i].positions[:] =  np.array([xcent,ycent,zcent])
    for i in range(nstars_bulge):
        bulgestars_list[i].positions[:] =  np.array([xcent,ycent,zcent])
    for i in range(nstars_disk):
        diskstars_list[i].positions[:] = np.array([xcent,ycent,zcent])

if par.SOURCES_RANDOM_POSITIONS == True:
    print "================================"
    print "SETTING SOURCES TO RANDOM POSITIONS"
    print "================================"
    for i in range(nstars):
        xpos,ypos,zpos = np.random.uniform(-dx,dx),np.random.uniform(-dy,dy),np.random.uniform(-dz,dz)
        stars_list[i].positions[:] = np.array([xpos,ypos,zpos])
    for i in range(nstars_bulge):
        xpos,ypos,zpos = np.random.uniform(-dx,dx),np.random.uniform(-dy,dy),np.random.uniform(-dz,dz)
        bulgestars_list[i].positions[:] = np.array([xpos,ypos,zpos])
    for i in range(nstars_disk):
        xpos,ypos,zpos = np.random.uniform(-dx,dx),np.random.uniform(-dy,dy),np.random.uniform(-dz,dz)
        diskstars_list[i].positions[:] = np.array([xpos,ypos,zpos])
'''


print('Done adding Sources')


# set up the CMB field -- place holder to put in haardt/madau eventually
'''
cmb = m.add_external_box_source()
cmb.temperature = cfg.model.TCMB
cmb_box_len = ds.quan(cfg.par.zoom_box_len,'kpc').in_units('cm').value
cmb.bounds = [[-cmb_box_len,cmb_box_len],[-cmb_box_len,cmb_box_len],[-cmb_box_len,cmb_box_len]]
pdb.set_trace()
L_CMB = (constants.sigma_sb*(cfg.model.TCMB*u.K)**4.).to(u.erg/u.cm**2/u.s)*4*(cmb_box_len*u.cm)**2 #get_J_CMB()
cmb.luminosity = L_CMB.cgs.value
'''

'''
energy_density_absorbed=energy_density_absorbed_by_CMB()
m.add_density_grid(density, dust, specific_energy=energy_density_absorbed)
m.set_specific_energy_type('additional')
'''


print('Setting up Model')
m_imaging = copy.deepcopy(m)

if cfg.par.SED == True:
    # set up the SEDs and images

    if cfg.par.SED_MONOCHROMATIC == True:

         # since all sources have the same spectrum just take the nu
         # from the input SED from the first source

        monochromatic_nu = m.sources[0].spectrum['nu']*u.Hz
        monochromatic_lam = (constants.c/monochromatic_nu).to(u.micron).value[::-1]

        if cfg.par.FIX_SED_MONOCHROMATIC_WAVELENGTHS == True:
            # idx = np.round(np.linspace(np.min(np.where(monochromatic_lam > cfg.par.SED_MONOCHROMATIC_min_lam)[0]),\
            ##                           np.max(np.where(monochromatic_lam < cfg.par.SED_MONOCHROMATIC_max_lam)[0]),\
            #                           cfg.par.SED_MONOCHROMATIC_nlam))

            idx = np.where((monochromatic_lam > cfg.par.SED_MONOCHROMATIC_min_lam) & (monochromatic_lam < cfg.par.SED_MONOCHROMATIC_max_lam))[0]
            monochromatic_lam = np.take(monochromatic_lam, list(idx))
        m.set_monochromatic(True, wavelengths=monochromatic_lam)
        m.set_raytracing(True)
        m.set_n_photons(initial=par.n_photons_initial,
                        imaging_sources=par.n_photons_imaging,
                        imaging_dust=par.n_photons_imaging,
                        raytracing_sources=par.n_photons_raytracing_sources,
                        raytracing_dust=par.n_photons_raytracing_dust)

        m.set_n_initial_iterations(3)
        m.set_convergence(True, percentile=99., absolute=1.01, relative=1.01)
        sed = m.add_peeled_images(sed=True, image=False)

        if cfg.par.MANUAL_ORIENTATION == True:
            sed.set_viewing_angles(np.array(cfg.model.THETA), np.array(cfg.model.PHI))

        else:
            sed.set_viewing_angles(np.linspace(0, 90, par.NTHETA).tolist()*par.NPHI, np.repeat(np.linspace(0, 90, par.NPHI), par.NPHI))
        sed.set_track_origin('basic')

        if cfg.par.SKIP_RT == False:
            m.write(model.inputfile+'.sed', overwrite=True)
            m.run(model.outputfile+'.sed', mpi=True, n_processes=par.n_MPI_processes, overwrite=True)

        print('[pd_front_end]: Beginning RT Stage: Calculating SED using a monochromatic spectrum equal to the input SED')

    else:

        m.set_raytracing(True)
        m.set_n_photons(initial=par.n_photons_initial, imaging=par.n_photons_imaging, raytracing_sources=par.n_photons_raytracing_sources, raytracing_dust=par.n_photons_raytracing_dust)
        m.set_n_initial_iterations(7)
        m.set_convergence(True, percentile=99., absolute=1.01, relative=1.01)

        sed = m.add_peeled_images(sed=True, image=False)
        sed.set_wavelength_range(2500, 0.001, 1000.)

        if cfg.par.MANUAL_ORIENTATION == True:
            sed.set_viewing_angles(np.array(cfg.model.THETA), np.array(cfg.model.PHI))
        else:
            sed.set_viewing_angles(np.linspace(0, 90, par.NTHETA).tolist(
            )*par.NPHI, np.repeat(np.linspace(0, 90, par.NPHI), par.NPHI))
        sed.set_track_origin('basic')

        print('[pd_front_end]: Beginning RT Stage: Calculating SED using a binned spectrum')

        # Run the Model
        if cfg.par.SKIP_RT == False:
            m.write(model.inputfile+'.sed', overwrite=True)
            m.run(model.outputfile+'.sed', mpi=True,
                  n_processes=par.n_MPI_processes, overwrite=True)

# see if the variable exists to make code backwards compatible

if cfg.par.IMAGING == True:

    print("Beginning Monochromatic Imaging RT")

    if cfg.par.IMAGING_TRANSMISSION_FILTER == False:

        # read in the filters file
        try:
            filter_data = [np.loadtxt(par.filterdir+f) for f in par.filterfiles]
        except:
            raise ValueError("Filters not found. You may be running above changeset 'f1f16eb' with an outdated parameters_master file. Please update to the most recent parameters_master format or ensure that the 'filterdir' and 'filterfiles' parameters are set properly.")

        # Extract and flatten all wavelengths in the filter files
        wavs = [wav[0] for single_filter in filter_data for wav in single_filter]
        wavs = list(set(wavs))      # Remove duplicates, if they exist

        m_imaging.set_monochromatic(True, wavelengths=wavs)
        m_imaging.set_raytracing(True)
        m_imaging.set_n_photons(initial=par.n_photons_initial,
                                imaging_sources=par.n_photons_imaging,
                                imaging_dust=par.n_photons_imaging,
                                raytracing_sources=par.n_photons_raytracing_sources,
                                raytracing_dust=par.n_photons_raytracing_dust)
    else:
        m_imaging.set_n_photons(initial=par.n_photons_initial, imaging=par.n_photons_imaging)

    m_imaging.set_n_initial_iterations(7)
    m_imaging.set_convergence(True, percentile=99., absolute=1.01, relative=1.01)

    image = m_imaging.add_peeled_images(sed=True, image=True)

    if cfg.par.IMAGING_TRANSMISSION_FILTER == True:
        add_transmission_filters(image)

    if cfg.par.MANUAL_ORIENTATION == True:
        image.set_viewing_angles(np.array(cfg.model.THETA), np.array(cfg.model.PHI))
    else:
        image.set_viewing_angles(np.linspace(0, 90, par.NTHETA).tolist()*par.NPHI, np.repeat(np.linspace(0, 90, par.NPHI), par.NPHI))

    image.set_track_origin('basic')
    image.set_image_size(cfg.par.npix_x, cfg.par.npix_y)
    image.set_image_limits(-dx/2., dx/2., -dy/2., dy/2.)

    if cfg.par.SKIP_RT == False:
        m_imaging.write(model.inputfile+'.image', overwrite=True)
        m_imaging.run(model.outputfile+'.image', mpi=True, n_processes=par.n_MPI_processes, overwrite=True)

        convolve(model.outputfile+'.image', par.filterfiles, filter_data)
    else:
        # Print a message in case that skip_rt debugging flag is set:
        print('++++++++++++++++++++++++++++++++++++')
        print('WARNING: SKIP RT is set in the parameters_master file - this is why your code didnt run')
        print('++++++++++++++++++++++++++++++++++++')

if ds_type in ['gadget_hdf5','tipsy','arepo_hdf5']:
    dump_data(reg, model)
