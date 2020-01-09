import numpy as np
import yt
import powderday.config as cfg
from yt.config import ytcfg
from yt.data_objects.particle_filters import add_particle_filter

ytcfg["yt","skip_dataset_cache"] = "True"


def tipsy_field_add(fname,bounding_box = None ,ds=None,add_smoothed_quantities=True):

    def _starmetals(field,data):
        return data[('newstars', 'Metals')]
    
    def _starcoordinates(field,data):
        return data[('newstars', 'Coordinates')]

    def _starformationtime(field,data):
        return data[('newstars', 'FormationTime')]

    def _stellarages(field,data):
        ad = data.ds.all_data()
        simtime = data.ds.current_time.in_units('Gyr')
        simtime = simtime.value
        age = simtime - ad["starformationtime"].in_units('Gyr').value
        age[np.where(age < 1.e-3)[0]] = 1.e-3
        age = data.ds.arr(age,'Gyr')
        return age

    def _starmasses(field,data):
        return data[('newstars', 'Mass')]


    def _diskstarcoordinates(field,data):
        return data[('diskstars','Coordinates')]
        
    def _diskstarmasses(field,data):
        return data[("diskstars","Mass")]
    
    
    def _gasdensity(field,data):
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

    def _gasmasses(field,data):
        return data[('Gas','Mass')]

    def _metaldens(field,data):
        return (data["Gas","Density"]*data["Gas","Metals"])

    def _gasfh2(field,data):
        try: return data[('PartType0','FractionH2')] # change this
        except: return (data.ds.arr(data[("Gas", "Mass")].value*0, 'dimensionless'))

    def _gassfr(field,data):
        try: return data[('PartType0','StarFormationRate')] # change this
        except: return (data.ds.arr(data[("Gas", "Mass")].value*0, 'code_mass/code_time'))


    if fname != None:
        ds = yt.load(fname,over_refine_factor=cfg.par.oref,n_ref=cfg.par.n_ref)
        ds.index
    

    #set up particle filters to figure out which stars might have been
    #initalized with the simulation.  we'll call stars that are formed
    #in the simulation 'new_stars' (no matter how old they are), and
    #stars that are initalized with the simulation as 'diskstars'.
    #Assume that no cosmology someone uses will give us a Hubble time
    #> 15 Gyr.

    def newstars(pfilter, data):
        age = data.ds.current_time - data[pfilter.filtered_type, "creation_time"]
        filter = np.logical_and(age.in_units('Gyr') <= 15, age.in_units('Gyr') >= 0)
        return filter

    def diskstars(pfilter, data):
        age = data.ds.current_time - data[pfilter.filtered_type, "creation_time"]
        filter = np.logical_or(age.in_units('Gyr') >= 15, age.in_units('Gyr') <= 0)
        return filter


    add_particle_filter("newstars", function=newstars, filtered_type='Stars', requires=["creation_time"])
    add_particle_filter("diskstars", function=diskstars, filtered_type='Stars', requires=["creation_time"])

    ds.add_particle_filter("newstars")
    ds.add_particle_filter("diskstars")

    ad = ds.all_data()

    #add the regular powderday fields
    ds.add_field(('starmetals'),function=_starmetals,units="code_metallicity",particle_type=True)
    ds.add_field(('starcoordinates'),function=_starcoordinates,units='cm',particle_type=True)
    ds.add_field(('starformationtime'),function=_starformationtime,units='Gyr',particle_type=True)
    ds.add_field(('stellarages'),function=_stellarages,units='Gyr',particle_type=True)
    ds.add_field(('starmasses'),function=_starmasses,units='g',particle_type=True)
    ds.add_field(('gasdensity'),function=_gasdensity,units='g/cm**3',particle_type=True)
    ds.add_field(('gasmetals'),function=_gasmetals,units="code_metallicity",particle_type=True)
    ds.add_field(('gascoordinates'),function=_gascoordinates,units='cm',particle_type=True)
    if add_smoothed_quantities == True:
        ds.add_field(('gassmootheddensity'),function=_gassmootheddensity,units='g/cm**3',particle_type=True)
        ds.add_field(('gassmoothedmetals'),function=_gassmoothedmetals,units='code_metallicity',particle_type=True)
        ds.add_field(('gassmoothedmasses'),function=_gassmoothedmasses,units='g',particle_type=True)
    ds.add_field(('gasmasses'),function=_gasmasses,units='g',particle_type=True)
    ds.add_field(('gasfh2'),function=_gasfh2,units='dimensionless',particle_type=True)
    ds.add_field(('gassfr'),function=_gassfr,units='g/s',particle_type=True)
    ds.add_field(('metaldens'),function=_metaldens,units="g/cm**3", particle_type=True)
    #only add the disk star fields if there are any disk stars
    if len(ad["diskstars","Mass"]) > 0 :
        ds.add_field(('diskstarmasses'),function=_diskstarmasses,units='g',particle_type=True)
        ds.add_field(('diskstarcoordinates'),function=_diskstarcoordinates,units='cm',particle_type=True)



    ad = ds.all_data()
    return ds

