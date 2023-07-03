#FIELD NOMENCLATURE

#For the arepo front end, we only pull in particle information.  This
#is because in arepo_tributary, we end up re-building the Voronoi
#tesslation via the built in voro++ calls in Hyperion.  As a result,
#there's (for example) only one field nomenclature that looks like
#[('dust','mass')] which represents the particle information.

from __future__ import print_function
import numpy as np
import yt
#from yt.fields.particle_fields import add_volume_weighted_smoothed_field
import powderday.config as cfg
from powderday.mlt.dgr_extrarandomtree_part import dgr_ert
#from yt.data_objects.particle_filters import add_particle_filter

def arepo_field_add(fname, bounding_box=None, ds=None):

    def _gassmoothinglength(field,data):
        return data[('PartType0', 'smoothing_length')].in_units('pc')

    def _starmetals(field, data):
        return data[('newstars', 'GFM_Metallicity')]

    def _starcoordinates(field, data):
        return data[('newstars', 'Coordinates')]

    def _starformationtime(field, data):
        return data[('newstars', 'GFM_StellarFormationTime')]

    def _starmasses(field, data):
        return data[("newstars", "Masses")]

#    def _diskstarcoordinates(field, data):
#        return data[('PartType2', 'Coordinates')]

#    def _diskstarmasses(field, data):
#        return data[("PartType2", "Masses")]

#    def _bulgestarcoordinates(field, data):
#        return data[('PartType3', 'Coordinates')]

#    def _bulgestarmasses(field, data):
#        return data[("PartType3", "Masses")]

    def _gasdensity(field, data):
        return data[('PartType0', 'density')]

    def _gasmetals(field, data):
        return data[('PartType0', 'GFM_Metallicity')]

    def _gascoordinates(field, data):
        return data[('PartType0', 'Coordinates')]

    def _gasmasses(field, data):
        return data[('PartType0', 'Masses')]

    def _gasfh2(field, data):
        try: return data[('PartType0', 'FractionH2')]
        except: return data[('PartType0', 'GFM_Metallicity')]*0. #just some dimensionless array

    def _gassfr(field, data):
        return data[('PartType0', 'StarFormationRate')]

    def _metaldens(field, data):
        return (data["PartType0", "density"]*data["PartType0", "GFM_Metallicity"])

    def _metalmass(field, data):
        return (data["PartType0", "Masses"]*(data["PartType0", "GFM_Metallicity"].value))


    def _dustcoordinates(field, data):
        return data[('PartType3', 'Coordinates')]
        
    def _dustmass_manual(field, data):
        if cfg.par.otf_extinction == True:
            dust_mass = (data[('PartType3','Masses')]).in_units('code_mass')
            
            return dust_mass
        else:
            return (data.ds.arr(data[("PartType0", "Dust_Masses")].value, 'code_mass'))

    def _dustmass_dtm(field,data):
        return (data["PartType0","metalmass"]*cfg.par.dusttometals_ratio)

    #def _dust_numgrains(field,data):
    #    return (data[('PartType0', 'NumGrains')])


    def _li_ml_dustmass(field,data):
        li_ml_dgr = dgr_ert(data["PartType0","GFM_Metallicity"],data["PartType0","StarFormationRate"],data["PartType0","Masses"])
        li_ml_dustmass = ((10.**li_ml_dgr)*data["PartType0","Masses"]).in_units('code_mass')
        
        #ds.parameters['li_ml_dustmass'] = li_ml_dustmass
        ##return (data.ds.arr(data.ds.parameters['li_ml_dustmass'].value,'code_mass'))
        return li_ml_dustmass


    def _dustmass_rr(field,data):
        #hard coded values from remy-ruyer table 1
        a = 2.21
        alpha = 2.02
        x_sun = 8.69

        x = 12.+np.log10(data["PartType0","GFM_Metallicity"]/cfg.par.solar * 10.**(x_sun-12.) )
        y = a + alpha*(x_sun-np.asarray(x))
        gas_to_dust_ratio = 10.**(y)
        dust_to_gas_ratio = 1./gas_to_dust_ratio
        return dust_to_gas_ratio * data["PartType0","Masses"]

    def _dustmass_li_bestfit(field,data):
        log_dust_to_gas_ratio = (2.445*np.log10(data["PartType0","GFM_Metallicity"]/cfg.par.solar))-(2.029)
        dust_to_gas_ratio = 10.**(log_dust_to_gas_ratio)
        return dust_to_gas_ratio * data["PartType0","Masses"]


    def _particle_dust_carbon_fraction(field,data):
        return data['PartType3','Dust_MetalFractions'][:,2]

        #you need to be careful: flattening the ds via ad =
        #ds.all_data() and then returning this value (vs just getting
        #it from data[] as above) ends up in a bug where zoom does not
        #take only the cut out data for this field for reg downstream.

        #ad = data.ds.all_data() return
        #(ad['PartType3','Dust_MetalFractions'][:,2])

    def _particle_dust_numgrains(field,data):
        return data['PartType3','Dust_NumGrains']
        #ad = data.ds.all_data()
        #return (ad['PartType3','Dust_NumGrains'])

    def _particle_dust_mass(field,data):
        dust_mass = (data[('PartType3','Masses')]).in_units('code_mass')
        return dust_mass
        #ad = data.ds.all_data()
        #return (ad['PartType3','Masses'])

    def _particle_dust_coordinates(field,data):
        return data['PartType3','Coordinates']
        #ad = data.ds.all_data()
        #return (ad['PartType3','Coordinates'])

    def _stellarages(field, data):
        ad = data.ds.all_data()
        if data.ds.cosmological_simulation == False:

            simtime = data.ds.current_time.in_units('Gyr')
            simtime = simtime.value

            print(" ")
            print("------------------------------------------------------------------")
            print("WARNING WARNING WARNING:")
            print("Assuming units in stellar ages are s*kpc/km")
            print("if this is not true - please edit _stellarages in front_ends/arepo2pd.py right under this warning message")
            print("------------------------------------------------------------------")

            age = simtime-(data.ds.arr(ad[("newstars","GFM_StellarFormationTime")],'s*kpc/km').in_units('Gyr')).value
            # make the minimum age 1 million years
            age[np.where(age < 1.e-3)[0]] = 1.e-3

            print('\n--------------')
            print(
                '[arepo2pd: ] Idealized Galaxy Simulation Assumed: Simulation time is (Gyr): ', simtime)
            print('--------------\n')
        else:
            yt_cosmo = yt.utilities.cosmology.Cosmology(hubble_constant=data.ds.hubble_constant,
                                                        omega_matter=data.ds.omega_matter,
                                                        omega_lambda=data.ds.omega_lambda)
            simtime = yt_cosmo.t_from_z(ds.current_redshift).in_units('Gyr').value # Current age of the universe
            scalefactor = data[("newstars","GFM_StellarFormationTime")].value
            formation_z = (1./scalefactor)-1.
            formation_time = yt_cosmo.t_from_z(formation_z).in_units('Gyr').value
            age = simtime - formation_time
            # Minimum age is set to 1 Myr (FSPS doesn't work properly for ages below 1 Myr)
            age[np.where(age < 1.e-3)[0]] = 1.e-3

            print('\n--------------')
            print('[arepo2pd: ] Cosmological Galaxy Simulation Assumed: Current age of Universe is (Gyr): ', simtime)
            print('--------------\n')

        age = data.ds.arr(age, 'Gyr')
        return age
        
    '''
    def _bhluminosity(field, data):
        ad = data.ds.all_data()
        mdot = ad[("PartType5", "BH_Mdot")]
        # give it a unit since usually these are dimensionless in yt

        mdot = data.ds.arr(mdot, "code_mass/code_time")


        c = yt.utilities.physical_constants.speed_of_light_cgs
        bhluminosity = (cfg.par.BH_eta * mdot * c**2.).in_units("erg/s")
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
    '''

    # load the ds (but only if this is our first passthrough and we pass in fname)
    if fname != None:
        try:
            float(yt.__version__[0:3]) >= 4
            ds = yt.load(fname)
            ds.index
            ad = ds.all_data()
        except:
            raise ValueError("It appears as though you are running in yt3.x  The vornoi mesh cannot be read in yt3.x.  Please update to yt4.x following the instructions here: https://powderday.readthedocs.io/en/latest/installation.html#yt-4-x-configuration-wip")
    

    #set up particle_filters to figure out which particles are stars.
    #we'll call particles that have ages > 0 newstars.

    def _newstars(pfilter,data):
        filter = data[(pfilter.filtered_type, "GFM_StellarFormationTime")] > 0
        return filter

    yt.add_particle_filter("newstars",function=_newstars,filtered_type='PartType4')
    ds.add_particle_filter("newstars")
    
    ds.add_field(('star','metals'), function=_starmetals, sampling_type='particle',units="code_metallicity", particle_type=True)

    ds.add_field(('gas','metals'), function=_gasmetals, sampling_type='particle', units="code_metallicity", particle_type=True)
    ds.add_field(('metal','dens'), function=_metaldens,  sampling_type='particle',units="g/cm**3", particle_type=True)
    ds.add_field(('PartType0', 'metalmass'), function=_metalmass,  sampling_type='particle', units="g", particle_type=True)

    ds.add_field(('gas','density'), function=_gasdensity, sampling_type='particle', units='g',particle_type=True)
    ds.add_field(('gas','masses'), function=_gasmasses,  sampling_type='particle', units='g', particle_type=True)
    ds.add_field(('gas','fh2'), function=_gasfh2,  sampling_type='particle', units='dimensionless', particle_type=True)
    ds.add_field(('gas','sfr'), function=_gassfr,  sampling_type='particle', units='g/s', particle_type=True)
    ds.add_field(('gas','smoothinglength'),function=_gassmoothinglength, sampling_type='particle', units='pc',particle_type=True)

    # get the dust mass

    if cfg.par.dust_grid_type == 'dtm':
        ds.add_field(('dust','mass'), function=_dustmass_dtm,  sampling_type='particle', units='code_mass',particle_type=True)
    if cfg.par.dust_grid_type == 'manual':
        ds.add_field(('dust','mass'), function=_dustmass_manual,  sampling_type='particle', units='code_mass', particle_type=True)
                
        #ds.add_deposited_particle_field(("PartType0", "Dust_Masses"), "sum")
        #if add_smoothed_quantities == True: ds.add_field(('dustsmoothedmasses'), function=_dustsmoothedmasses,  sampling_type='particle', units='code_mass', particle_type=True)


    if cfg.par.dust_grid_type == 'rr':
        ds.add_field(("dust","mass"),function=_dustmass_rr, sampling_type='particle', units='code_mass',particle_type=True)
    if cfg.par.dust_grid_type == 'li_bestfit':
        ds.add_field(("dust","mass"),function=_dustmass_li_bestfit, sampling_type='particle', units='code_mass',particle_type=True)


    #if we have the Li, Narayanan & Dave 2019 Extreme Randomized Trees
    #dust model in place, create a field for these so that
    #dust_grid_gen can use these dust masses
    if cfg.par.dust_grid_type == 'li_ml':
        #get the dust to gas ratio
        #ad = ds.all_data()
        #li_ml_dgr = dgr_ert(ad["gasmetals"],ad["PartType0","StarFormationRate"],ad["PartType0","Masses"])
        #li_ml_dustmass = ((10.**li_ml_dgr)*ad["PartType0","Masses"]).in_units('code_mass')
        #this is an icky way to pass this to the function for ds.add_field in the next line. but such is life.
        #ds.parameters['li_ml_dustmass'] = li_ml_dustmass
        #ds.add_field(('li_ml_dustmass'),function=_li_ml_dustmass, sampling_type='particle', units='code_mass',particle_type=True)
        ds.add_field(("dust","mass"),function=_li_ml_dustmass, sampling_type='particle', units='code_mass',particle_type=True)



    #add the numgrains (in case we're in OTF extinction mode)
    if cfg.par.otf_extinction:
        ds.add_field(("dust","coordinates"),function=_dustcoordinates, sampling_type='particle',units='code_length',particle_type=True)
        #ds.add_field(('dust','numgrains'),function=_dust_numgrains,units='dimensionless',sampling_type='particle',particle_type=True)
        ds.add_field(('particle_dust','carbon_fraction'),function=_particle_dust_carbon_fraction,units='dimensionless',sampling_type='particle',particle_type=True)
        ds.add_field(('particle_dust','numgrains'),function=_particle_dust_numgrains,units='dimensionless',sampling_type='particle',particle_type=True)
        ds.add_field(('particle_dust','mass'),function=_particle_dust_mass,units='code_mass',sampling_type='particle',particle_type=True)
        ds.add_field(('particle_dust','coordinates'),function=_particle_dust_coordinates,units='code_length',sampling_type='particle',particle_type=True)


    ds.add_field(('star','masses'), function=_starmasses,  sampling_type='particle', units='g', particle_type=True)
    ds.add_field(('star','coordinates'), function=_starcoordinates,  sampling_type='particle', units='cm', particle_type=True)
    ds.add_field(('star','formationtime'), function=_starformationtime,  sampling_type='particle', units='dimensionless', particle_type=True)
    ds.add_field(('stellar','ages'),function=_stellarages, sampling_type='particle', units='Gyr',particle_type=True)


    ds.add_field(('gas','density'), function=_gasdensity, sampling_type='particle',units='g/cm**3', particle_type=True)
    # Gas Coordinates need to be in Comoving/h as they'll get converted later.
    ds.add_field(('gas','coordinates'), function=_gascoordinates, sampling_type='particle',units='cm', particle_type=True)


    if cfg.par.BH_SED == True:
        try:
            nholes = len(ds.all_data()[('PartType5', 'BH_Mass')])
            if nholes > 0:
                if cfg.par.BH_model == 'Nenkova':
                    from powderday.agn_models.nenkova import Nenkova2008
                    agn_spectrum = Nenkova2008(*cfg.par.nenkova_params).agn_spectrum
                else:
                    from powderday.agn_models.hopkins import agn_spectrum

                if cfg.par.BH_var:
                    from powderday.agn_models.hickox import vary_bhluminosity
                    cfg.par.bhlfrac = vary_bhluminosity(nholes)

                ds.add_field(("bh","luminosity"),function=_bhluminosity,sampling_type='particle',units='erg/s',particle_type=True)
                ds.add_field(("bh","coordinates"),function=_bhcoordinates,sampling_type='particle',units="cm",particle_type=True)
                ds.add_field(("bh","nu"),function=_bhsed_nu,sampling_type='particle',units='Hz',particle_type=True)
                ds.add_field(("bh","sed"),function=_bhsed_sed,sampling_type='particle',units="erg/s",particle_type=True)
            else:
                print('No black holes found (length of BH_Mass field is 0)')
        except:
            print('Unable to find field "BH_Mass" in snapshot. Skipping.')



    return ds
