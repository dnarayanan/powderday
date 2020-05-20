from __future__ import print_function
import numpy as np
import yt
from yt.fields.particle_fields import add_volume_weighted_smoothed_field
import powderday.config as cfg
from powderday.mlt.dgr_extrarandomtree_part import dgr_ert

def gadget_field_add(fname, bounding_box=None, ds=None,add_smoothed_quantities=True):

    def _starmetals_00(field, data):
        el_dict = {'He': '01',
                   'C': '02',
                   'N': '03',
                   'O': '04',
                   'Ne': '05',
                   'Mg': '06',
                   'Si': '07',
                   'S': '08',
                   'Ca': '09',
                   'Fe': '10'}
        el_str = field.name[1]
        if '_' in el_str:
            el_name = field.name[1][field.name[1].find('_')+1:]
            el_num = el_dict[el_name]
        else:
            el_num = '00'
        return data[('PartType4', 'Metallicity_'+el_num)]

    def _starmetals(field, data):
        return data[('PartType4', 'Metallicity')]

    def _starcoordinates(field, data):
        return data[('PartType4', 'Coordinates')]

    def _starformationtime(field, data):
        return data[('PartType4', 'StellarFormationTime')]

    def _starmasses(field, data):
        return data[("PartType4", "Masses")]

    def _diskstarcoordinates(field, data):
        return data[('PartType2', 'Coordinates')]

    def _diskstarmasses(field, data):
        return data[("PartType2", "Masses")]

    def _bulgestarcoordinates(field, data):
        return data[('PartType3', 'Coordinates')]

    def _bulgestarmasses(field, data):
        return data[("PartType3", "Masses")]

    def _gasdensity(field, data):
        return data[('PartType0', 'Density')]

    def _gasmetals_00(field, data):
        return data[('PartType0', 'Metallicity_00')]

    def _gasmetals(field, data):
        return data[('PartType0', 'Metallicity')]

    def _gascoordinates(field, data):
        return data[('PartType0', 'Coordinates')]

    def _gasmasses(field, data):
        return data[('PartType0', 'Masses')]

    def _gasfh2(field, data):
        try: return data[('PartType0', 'FractionH2')]
        except: return data[('PartType0', 'metallicity')]*0.

    def _gassfr(field, data):
        return data[('PartType0', 'StarFormationRate')]

    def _gassmootheddensity(field, data):
        if  yt.__version__ == '4.0.dev0':
            return data.ds.parameters['octree'][('PartType0', 'density')]
        else:
            return data[("deposit","PartType0_smoothed_density")]

    def _gassmoothedmetals(field, data):
        if  yt.__version__ == '4.0.dev0':
            return data.ds.parameters['octree'][('PartType0', 'metallicity')]
        else:
            return data[("deposit","PartType0_smoothed_metallicity")]

    def _gassmoothedmasses(field, data):
        if  yt.__version__ == '4.0.dev0':
            return data.ds.parameters['octree'][('PartType0', 'Masses')]
        else:
            return data[('deposit', 'PartType0_mass')]

    def _metaldens_00(field, data):
        return (data["PartType0", "Density"]*data["PartType0", "Metallicity_00"])

    def _metaldens(field, data):
        return (data["PartType0", "Density"]*data["PartType0", "Metallicity"])

    def _metalmass_00(field, data):
        return (data["PartType0", "Masses"]*(data["PartType0", "Metallicity_00"].value))

    def _metalmass(field, data):
        return (data["PartType0", "Masses"]*(data["PartType0", "Metallicity"].value))

    def _metalsmoothedmasses(field, data):
        if  yt.__version__ == '4.0.dev0':
            return (data.ds.parameters['octree'][('PartType0', 'Masses')]* data.ds.parameters['octree'][('PartType0','metallicity')])
        else:
            return (data[('deposit', 'PartType0_smoothed_metalmass')].value)

        
    def _dustmass_manual(field, data):
        return (data.ds.arr(data[("PartType0", "Dust_Masses")].value, 'code_mass'))

    def _dustmass_dtm(field,data):
        return (data["PartType0","metalmass"]*cfg.par.dusttometals_ratio)

    def _li_ml_dustmass(field,data):
        return (data.ds.arr(data.ds.parameters['li_ml_dustmass'].value,'code_mass'))

    def _dustsmoothedmasses(field, data):
        if yt.__version__ == '4.0.dev0':
            return (data.ds.parameters['octree'][('PartType0','Dust_Masses')])

        else:
            return (data.ds.arr(data[("deposit", "PartType0_sum_Dust_Masses")].value, 'code_mass'))


    def _li_ml_dustsmoothedmasses(field,data):
        if yt.__version__ == '4.0.dev0':
            return (data.ds.parameters['octree'][('PartType0', 'li_ml_dustmass')])
        else:
            return (data.ds.arr(data[("deposit", "PartType0_sum_li_ml_dustmass")].value, 'code_mass'))
        



    def _stellarages(field, data):
        ad = data.ds.all_data()
        if data.ds.cosmological_simulation == False:
            simtime = data.ds.current_time.in_units('Gyr')
            simtime = simtime.value

            age = simtime-data.ds.arr(ad[('PartType4', 'StellarFormationTime')],'Gyr').value
            # make the minimum age 1 million years
            age[np.where(age < 1.e-3)[0]] = 1.e-3

            print('\n--------------')
            print(
                '[gadget2pd: ] Idealized Galaxy Simulation Assumed: Simulation time is (Gyr): ', simtime)
            print('--------------\n')
        else:
            yt_cosmo = yt.utilities.cosmology.Cosmology(hubble_constant=data.ds.hubble_constant,
                                                        omega_matter=data.ds.omega_matter,
                                                        omega_lambda=data.ds.omega_lambda)
            simtime = yt_cosmo.t_from_z(ds.current_redshift).in_units('Gyr').value # Current age of the universe
            scalefactor = data[('PartType4', 'StellarFormationTime')].value
            formation_z = (1./scalefactor)-1.
            formation_time = yt_cosmo.t_from_z(formation_z).in_units('Gyr').value
            age = simtime - formation_time
            # Minimum age is set to 1 Myr (FSPS doesn't work properly for ages below 1 Myr)
            age[np.where(age < 1.e-3)[0]] = 1.e-3

            print('\n--------------')
            print('[gadget2pd: ] Cosmological Galaxy Simulation Assumed: Current age of Universe is (Gyr): ', simtime)
            print('--------------\n')

        age = data.ds.arr(age, 'Gyr')
        return age
        
    def _starsmoothedmasses(field, data):
        return data[('deposit', 'PartType4_mass')]

    def _bhluminosity(field, data):
        ad = data.ds.all_data()
        mdot = ad[("PartType5", "BH_Mdot")]
        # give it a unit since usually these are dimensionless in yt

        mdot = data.ds.arr(mdot, "code_mass/code_time")


        c = yt.utilities.physical_constants.speed_of_light_cgs
        bhluminosity = (cfg.par.BH_eta * mdot * c**2.).in_units("erg/s")

        print("[front_ends/gadget2pd:] Generating the black hole luminosity")
        if cfg.par.BH_var:
            return bhluminosity * cfg.par.bhlfrac
        else:
            return bhluminosity

    def _bhcoordinates(field, data):
        return data["PartType5", "Coordinates"]

    def _bhsed_nu(field, data):
        bhluminosity = data["bhluminosity"]
        log_lum_lsun = np.log10(bhluminosity[0].in_units("Lsun"))
        nu, bhlum = agn_spectrum(log_lum_lsun)
        # the last 4 numbers aren't part of the SED
        nu = nu[0:-4]
        nu = 10.**nu
        nu = yt.YTArray(nu, "Hz")

        return nu

    def _bhsed_sed(field, data):
        bhluminosity = data["bhluminosity"]
        nholes = len(bhluminosity)

        # get len of nu just for the 0th hole so we know how long the vector is
        log_lum_lsun = np.log10(bhluminosity[0].in_units("Lsun"))
        nu, l_band_vec = agn_spectrum(log_lum_lsun)
        nu = nu[0:-4]
        n_nu = len(nu)

        bh_sed = np.zeros([nholes, n_nu])

        for i in range(nholes):

            log_lum_lsun = np.log10(bhluminosity[i].in_units("Lsun"))
            nu, l_band_vec = agn_spectrum(log_lum_lsun)

            l_band_vec = 10.**l_band_vec
            l_band_vec = l_band_vec[0:-4]
            for l in range(len(l_band_vec)):
                l_band_vec[l] = data.ds.quan(l_band_vec[l], "erg/s")

            bh_sed[i, :] = l_band_vec
        bh_sed = yt.YTArray(bh_sed, "erg/s")
        return bh_sed

    # load the ds (but only if this is our first passthrough and we pass in fname)
    if fname != None:
        if  yt.__version__ == '4.0.dev0':
            ds = yt.load(fname)
            ds.index
            ad = ds.all_data()
        else:
            ds = yt.load(fname,bounding_box=bounding_box,over_refine_factor=cfg.par.oref,n_ref=cfg.par.n_ref)
            ds.index
            ad = ds.all_data()
    
    #if we're in the 4.x branch of yt, load up the octree for smoothing
    if  yt.__version__ == '4.0.dev0':
        left = np.array([pos[0] for pos in bounding_box])
        right = np.array([pos[1] for pos in bounding_box])
        #octree = ds.octree(left, right, over_refine_factor=cfg.par.oref, n_ref=cfg.par.n_ref, force_build=True)
        octree = ds.octree(left,right,n_ref=cfg.par.n_ref)
        ds.parameters['octree'] = octree
        

    # for the metal fields have a few options since gadget can have different nomenclatures
    ad = ds.all_data()
    if ('PartType4', 'Metallicity_00') in ds.derived_field_list:
        try:
            ds.add_field(('starmetals'), function=_starmetals_00, units="code_metallicity", particle_type=True)
            ds.add_field(('starmetals_He'), function=_starmetals_00, units="code_metallicity", particle_type=True)
            ds.add_field(('starmetals_C'), function=_starmetals_00, units="code_metallicity", particle_type=True)
            ds.add_field(('starmetals_N'), function=_starmetals_00, units="code_metallicity", particle_type=True)
            ds.add_field(('starmetals_O'), function=_starmetals_00, units="code_metallicity", particle_type=True)
            ds.add_field(('starmetals_Ne'), function=_starmetals_00, units="code_metallicity", particle_type=True)
            ds.add_field(('starmetals_Mg'), function=_starmetals_00, units="code_metallicity", particle_type=True)
            ds.add_field(('starmetals_Si'), function=_starmetals_00, units="code_metallicity", particle_type=True)
            ds.add_field(('starmetals_S'), function=_starmetals_00, units="code_metallicity", particle_type=True)
            ds.add_field(('starmetals_Ca'), function=_starmetals_00, units="code_metallicity", particle_type=True)
            ds.add_field(('starmetals_Fe'), function=_starmetals_00, units="code_metallicity", particle_type=True)
        except:
            ds.add_field(('starmetals'), function=_starmetals_00, units="code_metallicity", particle_type=True)
    else:
        ds.add_field(('starmetals'), function=_starmetals, units="code_metallicity", particle_type=True)

    if ('PartType0', 'Metallicity_00') in ds.derived_field_list:
        ds.add_field(('gasmetals'), function=_gasmetals_00, units="code_metallicity", particle_type=True)
        ds.add_field(('metaldens'), function=_metaldens_00, units="g/cm**3", particle_type=True)
        # we add this as part type 0 (a non-general name) as it gets
        # smoothed immediately and that's all we end up using downstream
        ds.add_field(('PartType0', 'metalmass'), function=_metalmass_00, units="g", particle_type=True)

    else:
        ds.add_field(('gasmetals'), function=_gasmetals, units="code_metallicity", particle_type=True)
        ds.add_field(('metaldens'), function=_metaldens, units="g/cm**3", particle_type=True)
        ds.add_field(('PartType0', 'metalmass'), function=_metalmass, units="g", particle_type=True)

    metalmass_fn = add_volume_weighted_smoothed_field("PartType0", "Coordinates", "Masses",
                                                 "SmoothingLength", "Density", "metalmass",
                                                 ds.field_info)

    if add_smoothed_quantities == True: ds.add_field(('metalsmoothedmasses'), function=_metalsmoothedmasses, units='code_metallicity', particle_type=True)

    # get the dust mass

    if cfg.par.dust_grid_type == 'dtm':
        ds.add_field(('dustmass'), function=_dustmass_dtm,units='code_mass',particle_type=True)
    if cfg.par.dust_grid_type == 'manual':
        #if ('PartType0', 'Dust_Masses') in ds.derived_field_list:
        ds.add_field(('dustmass'), function=_dustmass_manual, units='code_mass', particle_type=True)
        ds.add_deposited_particle_field(("PartType0", "Dust_Masses"), "sum")
        if add_smoothed_quantities == True: ds.add_field(('dustsmoothedmasses'), function=_dustsmoothedmasses, units='code_mass', particle_type=True)

    #if we have the Li, Narayanan & Dave 2019 Extreme Randomized Trees
    #dust model in place, create a field for these so that
    #dust_grid_gen can use these dust masses
    if cfg.par.dust_grid_type == 'li_ml':
        #get the dust to gas ratio
        ad = ds.all_data()
        li_ml_dgr = dgr_ert(ad["PartType0","Metallicity_00"],ad["PartType0","StarFormationRate"],ad["PartType0","Masses"])
        li_ml_dustmass = ((10.**li_ml_dgr)*ad["PartType0","Masses"]).in_units('code_mass')
        #this is an icky way to pass this to the function for ds.add_field in the next line. but such is life.
        ds.parameters['li_ml_dustmass'] = li_ml_dustmass
        ds.add_field(('PartType0','li_ml_dustmass'),function=_li_ml_dustmass,units='code_mass',particle_type=True)
        ds.add_deposited_particle_field(("PartType0","li_ml_dustmass"),"sum")
        if add_smoothed_quantities == True: 
            ds.add_field(("li_ml_dustsmoothedmasses"), function=_li_ml_dustsmoothedmasses, units='code_mass',particle_type=True)

    ds.add_field(('starmasses'), function=_starmasses, units='g', particle_type=True)
    ds.add_field(('starcoordinates'), function=_starcoordinates, units='cm', particle_type=True)
    ds.add_field(('starformationtime'), function=_starformationtime, units='dimensionless', particle_type=True)
    
    ds.add_field(('stellarages'),function=_stellarages,units='Gyr',particle_type=True)
   
    if ('PartType2', 'Masses') in ds.derived_field_list:
        ds.add_field(('diskstarmasses'), function=_diskstarmasses, units='g', particle_type=True)
        ds.add_field(('diskstarcoordinates'), function=_diskstarcoordinates, units='cm', particle_type=True)

    if ('PartType3', 'Masses') in ds.derived_field_list:
        ds.add_field(('bulgestarmasses'), function=_bulgestarmasses, units='g', particle_type=True)
        ds.add_field(('bulgestarcoordinates'), function=_bulgestarcoordinates, units='cm', particle_type=True)

    if add_smoothed_quantities == True: ds.add_field(('starsmoothedmasses'), function=_starsmoothedmasses, units='g', particle_type=True)

    ds.add_field(('gasdensity'), function=_gasdensity, units='g/cm**3', particle_type=True)
    # Gas Coordinates need to be in Comoving/h as they'll get converted later.
    ds.add_field(('gascoordinates'), function=_gascoordinates, units='cm', particle_type=True)
    if add_smoothed_quantities == True:
        ds.add_field(('gassmootheddensity'), function=_gassmootheddensity, units='g/cm**3', particle_type=True)
        ds.add_field(('gassmoothedmetals'), function=_gassmoothedmetals, units='code_metallicity', particle_type=True)
        ds.add_field(('gassmoothedmasses'), function=_gassmoothedmasses, units='g', particle_type=True)



    ds.add_field(('gasmasses'), function=_gasmasses, units='g', particle_type=True)
    ds.add_field(('gasfh2'), function=_gasfh2, units='dimensionless', particle_type=True)
    ds.add_field(('gassfr'), function=_gassfr, units='g/s', particle_type=True)

    if cfg.par.BH_SED == True:
        if ('PartType5','BH_Mass') in ds.derived_field_list:
            nholes = len(ds.all_data()[('PartType5', 'BH_Mass')])
            print("The number of black holes is:",nholes)
            if nholes > 0:
                if cfg.par.BH_model == 'Nenkova':
                    from powderday.agn_models.nenkova import Nenkova2008
                    agn_spectrum = Nenkova2008(*cfg.par.nenkova_params).agn_spectrum
                else:
                    from powderday.agn_models.hopkins import agn_spectrum

                if cfg.par.BH_var:
                    from powderday.agn_models.hickox import vary_bhluminosity
                    cfg.par.bhlfrac = vary_bhluminosity(nholes)

                ds.add_field(("bhluminosity"),function=_bhluminosity,units='erg/s',particle_type=True)
                ds.add_field(("bhcoordinates"),function=_bhcoordinates,units="cm",particle_type=True)
                ds.add_field(("bhnu"),function=_bhsed_nu,units='Hz',particle_type=True)
                ds.add_field(("bhsed"),function=_bhsed_sed,units="erg/s",particle_type=True)
                


            else:
                print('No black holes found (length of BH_Mass field is 0)')


    return ds
