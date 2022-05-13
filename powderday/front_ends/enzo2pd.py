from __future__ import print_function
import numpy as np
import yt
import powderday.config as cfg

from yt.data_objects.particle_filters import add_particle_filter


def enzo_field_add(fname,ds = None, starages = False):

    def _starmetals(field,data):
        return data[('newstars','metallicity_fraction')]

    def _starcoordinates(field,data):
        #set units of cm then tack them back on because the column_stack loses them
        xpos = data[ ('newstars', 'particle_position_x')].in_units("cm")
        ypos = data[ ('newstars', 'particle_position_y')].in_units("cm")
        zpos = data[ ('newstars', 'particle_position_z')].in_units("cm")
        coordinates = np.column_stack((xpos,ypos,zpos))
        coordinates = data.ds.arr(coordinates,"cm")
        return coordinates

    def _stellarages(field,data):
        age = ds.current_time.in_units('Gyr')-data[('newstars', 'creation_time')].in_units('Gyr')
        age[np.where(age < 1.e-3)[0]] = 1.e-3
        return age

    def _starmasses(field,data):
        return data[('newstars', 'particle_mass')]

    def _gasdensity(field,data):
        return data[('gas', 'density')]
        
    def _gasmetals(field,data):
        return data[ ('gas', 'metallicity')]

    def _gasmasses(field,data):
        return data[('gas','cell_mass')]

    def _gasfh2(field, data):
        try: return data[('gas', 'FractionH2')]
        except: return data[('gas', 'metallicity')]*0. #just some dimensionless array
        

    #load the ds
    if fname != None:
        ds = yt.load(fname)
        ds.index

    #set up particle_filters to figure out which particles are stars.
    #we'll call particles that have ages > 0 stars.

    def newstars(pfilter,data):
        age = data[pfilter.filtered_type,"creation_time"]
        filter = age.in_units('Gyr') > 0
        return filter


        
    add_particle_filter("newstars",function=newstars,filtered_type='all',requires=["creation_time"])
    ds.add_particle_filter("newstars")
    ad = ds.all_data()




    ds.add_field(('star','metals'),function=_starmetals,units="code_metallicity",particle_type=True)
    ds.add_field(('star','coordinates'),function=_starcoordinates,units="cm",particle_type=True)
    ds.add_field(('stellar','ages'),function=_stellarages,units='Gyr',particle_type=True)
    ds.add_field(('star','masses'),function=_starmasses,units='g',particle_type=True)
    ds.add_field(('gas','density'),function=_gasdensity,units='g/cm**3')
    ds.add_field(('gas','metals'),function=_gasmetals,units="code_metallicity")
    ds.add_field(('gas','fh2'),function=_gasfh2,units='dimensionless')
    ds.add_field(('gas','masses'),function=_gasmasses,units='g')
    
    ad = ds.all_data()

    return ds
