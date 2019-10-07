from __future__ import print_function
import numpy as np
import yt
from yt.fields.particle_fields import add_volume_weighted_smoothed_field
import powderday.config as cfg

from astropy.cosmology import Planck13
import astropy.units as u
from powderday.front_ends.redshift_multithread import redshift_vectorized


def gadget_field_add(fname,bounding_box = None,ds=None,starages=False):
    
    
    def _starcoordinates(field,data):
        return data[('PartType4', 'Coordinates')]

    def _starformationtime(field,data):
        return data[('PartType4', 'StellarFormationTime')]

    def _starmasses(field,data):
        return data[("PartType4","Masses")]

    def _starmetals_00_CS(field,data):
        ad = ds.all_data()
        starmass = ad[("PartType4","Masses")].value
        metals = np.zeros(len(ad[("PartType4","Metallicity_00")]))
        for i in range(1,6):
            metals += ad[("PartType4","Metallicity_0%s"%i)]
            #we avoid 6 because this is the hydrogen mass
        for i in range(7,10):
            metals += ad[("PartType4","Metallicity_0%s"%i)]
        for i in (10,11):
            metals += ad[("PartType4","Metallicity_%s"%i)]

        #the CS metals are actually masses in code units, so to get
        #metallicity, divide by particle gas mass in code units
        
        metals /= starmass
        return data.ds.arr(metals,'code_metallicity')

    def _diskstarcoordinates(field,data):
        return data[('PartType2','Coordinates')]
    
    def _diskstarmasses(field,data):
        return data[("PartType2","Masses")]
    
    def _bulgestarcoordinates(field,data):
        return data[('PartType3','Coordinates')]

    def _bulgestarmasses(field,data):
        return data[("PartType3","Masses")]

    def _gasdensity(field,data):
        return data[('PartType0', 'Density')]

    def _gasmetals_00_CS(field,data):
        ad = ds.all_data()
        gasmass = ad[("PartType0","Masses")].value
        metals = np.zeros(len(ad[("PartType0","Metallicity_00")]))
        for i in range(1,6):
            metals += ad[("PartType0","Metallicity_0%s"%i)]
        #we avoid 6 because this is the hydrogen mass
        for i in range(7,10):
            metals += ad[("PartType0","Metallicity_0%s"%i)]
        for i in (10,11):
            metals += ad[("PartType0","Metallicity_%s"%i)]

        #the CS metals are actually masses in code units, so to get
        #metallicity, divide by particle gas mass in code units
        
        metals /= gasmass
        return data.ds.arr(metals,'code_metallicity')
        
    def _gasmasses(field,data):
        return data[('PartType0','Masses')]
        
    def _gasfh2(field,data):
        try: return data[('PartType0','FractionH2')]
        except: return np.zeros(len(data[('PartType0','Masses')]))

    def _gassfr(field,data):
        return data[('PartType0','StarFormationRate')]
        
    def _gascoordinates(field,data):
        return data[('PartType0','Coordinates')]

    def _gassmootheddensity(field,data):
        return data[("deposit","PartType0_smoothed_density")]

    def _gassmoothedmetals(field,data):
        return data[('deposit', 'PartType0_smoothed_gasmetals')]


    def _gassmoothedmasses(field,data):
        return data[('deposit', 'PartType0_mass')]
    
    def _metaldens_00_CS(field,data):
        return (data["PartType0","Density"]*data["gasmetals"].value)
    
    def _metalmass_00_CS(field,data):
        return (data["PartType0","Masses"]*(data["gasmetals"].value))
    
    def _metalsmoothedmasses(field,data):
        return (data[('deposit', 'PartType0_smoothed_metalmass')].value)

    def _stellarages(field,data):
        ad = data.ds.all_data()
        if data.ds.cosmological_simulation == False:

            simtime = data.ds.current_time.in_units('Gyr')
            simtime = simtime.value

            age = simtime-ad[("starformationtime")].value #Gyr (assumes that ad["starformationtime"] is in Gyr for Gadget)
            #make the minimum age 1 million years 
            age[np.where(age < 1.e-3)[0]] = 1.e-3
            
            
            print ('\n--------------')
            print ('[SED_gen/star_list_gen: ] Idealized Galaxy Simulation Assumed: Simulation time is (Gyr): ',simtime)
            print ('--------------\n')
        else:
            simtime = Planck13.age(data.ds.current_redshift).to(u.Gyr).value #what is the age of the Universe right now?
            
            scalefactor = ad[("starformationtime")].value
            formation_z = (1./scalefactor)-1.
            
            formation_time = redshift_vectorized(formation_z)
            #drop the Gyr unit
            formation_time = np.asarray([formation_time[i].value for i in range(len(formation_time))])
           
            age = simtime - formation_time
            #make the minimum age 1 million years 
            age[np.where(age < 1.e-3)[0]] = 1.e-3

        
            print ('\n--------------')
            print ('[SED_gen/star_list_gen: ] Cosmological Galaxy Simulation Assumed: Current age of Universe is (Assuming Planck13 Cosmology) is (Gyr): ',simtime)
            print ('--------------\n')
                
        age = data.ds.arr(age,'Gyr')
        return age

    def _starsmoothedmasses(field,data):
            return data[('deposit','PartType4_mass')]

        
    def _bhluminosity(field,data):
        ad = data.ds.all_data()
        mdot = ad[("PartType5","BH_Mdot")]
        #give it a unit since usually these are dimensionless in yt
        
        mdot = data.ds.arr(mdot,"code_mass/code_time")

        
        c = yt.utilities.physical_constants.speed_of_light_cgs
        bhluminosity = (cfg.par.BH_eta * mdot * c**2.).in_units("erg/s")
        if cfg.par.BH_var:
            from powderday.agn_models.hickox import vary_bhluminosity
            return vary_bhluminosity(bhluminosity)
        else:
            return bhluminosity
        
    def _bhcoordinates(field,data):
        return data["PartType5","Coordinates"]


    def _bhsed_nu(field,data):
        bhluminosity = data["bhluminosity"]
        log_lum_lsun = np.log10(bhluminosity[0].in_units("Lsun"))
        nu,bhlum = agn_spectrum(log_lum_lsun)
        #the last 4 numbers aren't part of the SED
        nu = nu[0:-4]
        nu = 10.**nu
        nu = yt.YTArray(nu,"Hz")
        
        return nu
        
    def _bhsed_sed(field,data):
        bhluminosity = data["bhluminosity"]
        nholes = len(bhluminosity)
        
        #get len of nu just for the 0th hole so we know how long the vector is
        log_lum_lsun = np.log10(bhluminosity[0].in_units("Lsun"))
        nu,l_band_vec = agn_spectrum(log_lum_lsun)
        nu = nu[0:-4]
        n_nu = len(nu)
        
        bh_sed = np.zeros([nholes,n_nu])

        for i in range(nholes):
            
            log_lum_lsun = np.log10(bhluminosity[i].in_units("Lsun"))
            nu,l_band_vec = agn_spectrum(log_lum_lsun)
            
            l_band_vec = 10.**l_band_vec
            l_band_vec = l_band_vec[0:-4]
            for l in range(len(l_band_vec)):
                l_band_vec[l] = data.ds.quan(l_band_vec[l],"erg/s")

            bh_sed[i,:] = l_band_vec
        bh_sed = yt.YTArray(bh_sed,"erg/s")
        return bh_sed

   #load the ds
    if fname != None:
        ds = yt.load(fname,bounding_box=bounding_box,over_refine_factor=cfg.par.oref,n_ref=cfg.par.n_ref)
        ds.index

    #------------------
    #Munich Group Gadget Metallicity Fields based on Cecila Scannapieco's Metallicity Implementation
    #------------------

                
    ds.add_field(('gasmetals'),function=_gasmetals_00_CS,units="code_metallicity",particle_type=True)
    ds.add_field(('starmetals'),function=_starmetals_00_CS,units="code_metallicity",particle_type=True)
    ds.add_field(('metaldens'),function=_metaldens_00_CS,units="g/cm**3", particle_type=True)
    ds.add_field(('PartType0','metalmass'),function=_metalmass_00_CS,units="g", particle_type=True)

    #this is exactly the same as 'gasmetals' but in 'parttype0'
    #notation so we can project it with
    #add_volume_weighted_smoothed_field
    ds.add_field((('PartType0','gasmetals')),function=_gasmetals_00_CS,units="code_metallicity",particle_type=True)








    metalmass_fn = add_volume_weighted_smoothed_field("PartType0", "Coordinates", "Masses",
                                                  "SmoothingLength", "Density","metalmass",
                                                  ds.field_info)

    metallicity_smoothed_fn = add_volume_weighted_smoothed_field("PartType0", "Coordinates", "Masses",
                                                  "SmoothingLength", "Density","gasmetals",
                                                  ds.field_info)

    
    ds.add_field(('metalsmoothedmasses'),function=_metalsmoothedmasses,units='code_metallicity',particle_type=True)



    ds.add_field(('starmasses'),function=_starmasses,units='g',particle_type=True)
    ds.add_field(('starcoordinates'),function=_starcoordinates,units='cm',particle_type=True)
    ds.add_field(('starformationtime'),function=_starformationtime,units='dimensionless',particle_type=True)
    ds.add_field(('starsmoothedmasses'),function=_starsmoothedmasses,units='g',particle_type=True)

    if ('PartType2','Masses') in ds.derived_field_list:
        ds.add_field(('diskstarmasses'),function=_diskstarmasses,units='g',particle_type=True)
        ds.add_field(('diskstarcoordinates'),function=_diskstarcoordinates,units='cm',particle_type=True)

    if ('PartType3','Masses') in ds.derived_field_list:
        ds.add_field(('bulgestarmasses'),function=_bulgestarmasses,units='g',particle_type=True)
        ds.add_field(('bulgestarcoordinates'),function=_bulgestarcoordinates,units='cm',particle_type=True)
    
    ds.add_field(('gasdensity'),function=_gasdensity,units='g/cm**3',particle_type=True)
    #Gas Coordinates need to be in Comoving/h as they'll get converted later.
    ds.add_field(('gascoordinates'),function=_gascoordinates,units='cm',particle_type=True)
    ds.add_field(('gassmootheddensity'),function=_gassmootheddensity,units='g/cm**3',particle_type=True)
    ds.add_field(('gassmoothedmetals'),function=_gassmoothedmetals,units='code_metallicity',particle_type=True)
    ds.add_field(('gassmoothedmasses'),function=_gassmoothedmasses,units='g',particle_type=True)
    ds.add_field(('gasmasses'),function=_gasmasses,units='g',particle_type=True)
    ds.add_field(('gasfh2'),function=_gasfh2,units='dimensionless',particle_type=True)
    ds.add_field(('gassfr'),function=_gassfr,units='g/s',particle_type=True)

    if cfg.par.BH_SED == True:
        try:
            if len(ds.all_data()[('PartType5', 'BH_Mass')]) > 0:
                if cfg.par.BH_model == 'Nenkova':
                    from powderday.agn_models.nenkova import Nenkova2008
                    try:
                        model = Nenkova2008(cfg.par.nenkova_params)
                    except:
                        model = Nenkova2008
                    agn_spectrum = model.agn_spectrum
                else:
                    from powderday.agn_models.hopkins import agn_spectrum

                ds.add_field(("bhluminosity"),function=_bhluminosity,units='erg/s',particle_type=True)
                ds.add_field(("bhcoordinates"),function=_bhcoordinates,units="cm",particle_type=True)
                ds.add_field(("bhnu"),function=_bhsed_nu,units='Hz',particle_type=True)
                ds.add_field(("bhsed"),function=_bhsed_sed,units="erg/s",particle_type=True)
            else:
                print('No black holes found (length of BH_Mass field is 0)')
        except:
            print('Unable to find field "BH_Mass" in snapshot. Skipping.')

    if starages == True:
        ds.add_field(('stellarages'),function=_stellarages,units='Gyr',particle_type=True)


    return ds
