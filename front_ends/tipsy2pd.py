import numpy as np
import astropy.units as u
import yt
from yt import derived_field
import ipdb
import config as cfg
from astropy.cosmology import Planck13
import astropy.units as u
from yt.config import ytcfg

ytcfg["yt","skip_dataset_cache"] = "True"

#need - 
#1. PartType4_Metallicity (newstar metals)
#2. PartType4_coordinates
#3. Parttype4_StelarFormationTiome
#4. PartType0_Density
#5. PartType0_Metallicity
#6. Parttype0_Coordinates
#7. Parttype0_Smoothed_Density

def tipsy_field_add(fname,bounding_box = None ,ds=None,starages=False):

    def _starmetals(field,data):
        return data[('Stars', 'Metals')]
    
    def _starcoordinates(field,data):
        return data[('Stars', 'Coordinates')]

    def _starformationtime(field,data):
        return data[('Stars', 'FormationTime')]

    def _stellarages(field,data):
        ad = data.ds.all_data()
        simtime = data.ds.current_time.in_units('Gyr')
        simtime = simtime.value
        age = simtime - ad["starformationtime"].value
        age[np.where(age < 1.e-3)[0]] = 1.e-3
        return age

    def _starmasses(field,data):
        return data[('Stars', 'Mass')]
    
    def _gasdensity(fied,data):
        return data[('Gas', 'Density')]
        
    def _gasmetals(field,data):
        return data[('Gas', 'Metals')]
        
    def _gascoordinates(field,data):
        return data[('Gas','Coordinates')]

    def _gassmootheddensity(field,data):
        return data[('deposit', 'Gas_density')]
    
    def _gassmoothedmetals(field,data):
        return data[('deposit', 'Gas_smoothed_metallicity')]
        
    def _gassmoothedmasses(field,data):
        return data[('deposit', 'Gas_mass')]

    def _metaldens(field,data):
        return (data["Gas","Density"]*data["Gas","Metals"])


    if fname != None:
        ds = yt.load(fname,bounding_box=bounding_box,over_refine_factor=cfg.par.oref,n_ref=cfg.par.n_ref)
        ds.index
    
    ds.add_field(('starmetals'),function=_starmetals,units="code_metallicity",particle_type=True)
    ds.add_field(('starcoordinates'),function=_starcoordinates,units='cm',particle_type=True)
    ds.add_field(('starformationtime'),function=_starformationtime,units='Gyr',particle_type=True)
    ds.add_field(('stellarages'),function=_stellarages,units='Gyr',particle_type=True)
    ds.add_field(('starmasses'),function=_starmasses,units='g',particle_type=True)
    ds.add_field(('gasdensity'),function=_gasdensity,units='g/cm**3',particle_type=True)
    ds.add_field(('gasmetals'),function=_gasmetals,units="code_metallicity",particle_type=True)
    ds.add_field(('gascoordinates'),function=_gascoordinates,units='cm',particle_type=True)
    ds.add_field(('gassmootheddensity'),function=_gassmootheddensity,units='g/cm**3',particle_type=True)
    ds.add_field(('gassmoothedmetals'),function=_gassmoothedmetals,units='code_metallicity',particle_type=True)
    ds.add_field(('metaldens'),function=_metaldens,units="g/cm**3", particle_type=True)
    ds.add_field(('gassmoothedmasses'),function=_gassmoothedmasses,units='g',particle_type=True)


    ad = ds.all_data()
    
    return ds

