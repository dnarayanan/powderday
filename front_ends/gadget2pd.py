import numpy as np
import yt
from yt import derived_field
from yt.fields.particle_fields import add_volume_weighted_smoothed_field
import pdb,ipdb
import config as cfg

from astropy.cosmology import Planck13
import astropy.units as u
from redshift_multithread import *



#need - 
#1. PartType4_Metallicity (newstar metals)
#2. PartType4_coordinates
#3. Parttype4_StelarFormationTiome
#4. PartType0_Density
#5. PartType0_Metallicity
#6. Parttype0_Coordinates
#7. Parttype0_Smoothed_Density


def gadget_field_add(fname,bounding_box = None,ds=None,starages=False):
    
    
    def _starmetals_00(field,data):
        return data[('PartType4', 'Metallicity_00')]
        
    def _starmetals(field,data):
        return data[('PartType4', 'Metallicity')]

    def _starcoordinates(field,data):
        return data[('PartType4', 'Coordinates')]

    def _starformationtime(field,data):
        return data[('PartType4', 'StellarFormationTime')]

    def _starmasses(field,data):
        return data[("PartType4","Masses")]

    def _diskstarcoordinates(field,data):
        return data[('PartType2','Coordinates')]
    
    def _diskstarmasses(field,data):
        return data[("PartType2","Masses")]
    
    def _bulgestarcoordinates(field,data):
        return data[('PartType3','Coordinates')]

    def _bulgestarmasses(field,data):
        return data[("PartType3","Masses")]

    def _gasdensity(fied,data):
        return data[('PartType0', 'Density')]

    def _gasmetals_00(field,data):
        return data[('PartType0', 'Metallicity_00')]
    
    def _gasmetals(field,data):
        return data[('PartType0', 'Metallicity')]

    def _gascoordinates(field,data):
        return data[('PartType0','Coordinates')]

    def _gassmootheddensity(field,data):
        return data[("deposit","PartType0_smoothed_density")]

    def _gassmoothedmetals(field,data):
        return data[("deposit","PartType0_smoothed_metallicity")]

    def _gassmoothedmasses(field,data):
        return data[('deposit', 'PartType0_mass')]
    
    def _metaldens_00(field,data):
        return (data["PartType0","Density"]*data["PartType0","Metallicity_00"])

    def _metaldens(field,data):
        return (data["PartType0","Density"]*data["PartType0","Metallicity"])

    def _metalmass_00(field,data):
        return (data["PartType0","Masses"]*(data["PartType0","Metallicity_00"].value))
    
    def _metalmass(field,data):
        return (data["PartType0","Masses"]*(data["PartType0","Metallicity"].value))

    def _metalsmoothedmasses(field,data):
        return (data[('deposit', 'PartType0_smoothed_metalmass')].value)

    def _stellarages(field,data):
        ad = data.ds.all_data()
        if cfg.par.COSMOFLAG == False:

            simtime = data.ds.current_time.in_units('Gyr')
            simtime = simtime.value
            
            age = simtime-ad[("starformationtime")].value #Gyr (assumes that ad["starformationtime"] is in Gyr for Gadget)
            #make the minimum age 1 million years 
            age[np.where(age < 1.e-3)[0]] = 1.e-3
            
            
            print '\n--------------'
            print '[SED_gen/star_list_gen: ] Idealized Galaxy Simulation Assumed: Simulation time is (Gyr): ',simtime
            print '--------------\n'
        else:
            simtime = Planck13.age(data.ds.current_redshift).to(u.Gyr).value #what is the age of the Universe right now?
            
            scalefactor = ad[("starformationtime")].value
            formation_z = (1./scalefactor)-1.
            
            formation_time = redshift_multithread(formation_z)
            
            age = simtime - formation_time
            #make the minimum age 1 million years 
            age[np.where(age < 1.e-3)[0]] = 1.e-3

        
            print '\n--------------'
            print '[SED_gen/star_list_gen: ] Cosmological Galaxy Simulation Assumed: Current age of Universe is (Assuming Planck13 Cosmology) is (Gyr): ',simtime
            print '--------------\n'
                
        return age
        
    #load the ds
    if fname != None:
        ds = yt.load(fname,bounding_box=bounding_box,over_refine_factor=cfg.par.oref,n_ref=cfg.par.n_ref)
        ds.index


    #for the metal fields have a few options since gadget can have different nomenclatures
    if  ('PartType4', 'Metallicity_00') in ds.derived_field_list:
        ds.add_field(('starmetals'),function=_starmetals_00,units="code_metallicity",particle_type=True)
    else:
        ds.add_field(('starmetals'),function=_starmetals,units="code_metallicity",particle_type=True)

        
    if  ('PartType0', 'Metallicity_00') in ds.derived_field_list:
        ds.add_field(('gasmetals'),function=_gasmetals_00,units="code_metallicity",particle_type=True)
        ds.add_field(('metaldens'),function=_metaldens_00,units="g/cm**3", particle_type=True)

        #we add this as part type 0 (a non-general name) as it gets
        #smoothed immediately and that's all we end up using downstream
        ds.add_field(('PartType0','metalmass'),function=_metalmass_00,units="g", particle_type=True)

    else:
        ds.add_field(('gasmetals'),function=_gasmetals,units="code_metallicity",particle_type=True)
        ds.add_field(('metaldens'),function=_metaldens,units="g/cm**3", particle_type=True)
        ds.add_field(('PartType0','metalmass'),function=_metalmass,units="g", particle_type=True)


    metalmass_fn = add_volume_weighted_smoothed_field("PartType0", "Coordinates", "Masses",
                                                  "SmoothingLength", "Density","metalmass",
                                                  ds.field_info)
    ds.add_field(('metalsmoothedmasses'),function=_metalsmoothedmasses,units='code_metallicity',particle_type=True)



    ds.add_field(('starmasses'),function=_starmasses,units='g',particle_type=True)
    ds.add_field(('starcoordinates'),function=_starcoordinates,units='cm',particle_type=True)
    ds.add_field(('starformationtime'),function=_starformationtime,units='dimensionless',particle_type=True)

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
    
    if starages == True:
        ds.add_field(('stellarages'),function=_stellarages,units='Gyr',particle_type=True)

    
    return ds
