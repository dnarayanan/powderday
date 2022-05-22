import numpy as np
import yt
import powderday.config as cfg
from yt.config import ytcfg
from yt.data_objects.particle_filters import add_particle_filter
from powderday.mlt.dgr_extrarandomtree_part import dgr_ert
from unyt import G

#ytcfg["yt","skip_dataset_cache"] = "True"


def tipsy_field_add(fname,bounding_box = None ,ds=None,add_smoothed_quantities=True):

    def _gassmoothinglength(field,data):
        return data[("Gas", "smoothing_length")].in_units('pc')

    def _starmetals(field,data):
        return data[("newstars", "Metals")]
    
    def _starcoordinates(field,data):
        return data[("newstars", "Coordinates")]

    def _starformationtime(field,data):
        return data[("newstars", "FormationTime")]

    def _stellarages(field,data):
        ad = data.ds.all_data()
        simtime = data.ds.current_time.in_units('Gyr')
        simtime = simtime.value
        age = simtime - data[("newstars","FormationTime")].in_units('Gyr').value
        age[np.where(age < 1.e-3)[0]] = 1.e-3
        age = data.ds.arr(age,'Gyr')
        return age

    def _starmasses(field,data):
        return data[("newstars", "Mass")]


    def _diskstarcoordinates(field,data):
        return data[("diskstars", "Coordinates")]
        
    def _diskstarmasses(field,data):
        return data[("diskstars", "Mass")]
    
    
    def _gasdensity(field,data):
        return data[("Gas", "Density")]

    def _gasmetals(field,data):
        return data[("Gas", "Metals")]
        
    def _gascoordinates(field,data):
        return data[("Gas", "Coordinates")]

    def _gassmootheddensity(field,data):
        if float(yt.__version__[0:3]) >= 4:
            return data.ds.parameters['octree'][("Gas", "Density")]
        else:
            return data[("deposit", "Gas_density")]
    
    def _gassmoothedmetals(field,data):
        if float(yt.__version__[0:3]) >= 4:
            return data.ds.parameters['octree'][("Gas", "Metals")]
        else:
            return data[("deposit", "Gas_smoothed_metallicity")]
            # Does this field exist?
        
    def _gassmoothedmasses(field,data):
        if float(yt.__version__[0:3]) >= 4:
            return data.ds.parameters['octree'][("Gas", "Mass")]
        else:
            return data[("deposit", "Gas_mass")]

    def _gasmasses(field,data):
        return data[("Gas","Mass")]

    def _dustmass_dtm(field,data):
        return (data[("Gas","Metals")]*data[("Gas","Mass")]*cfg.par.dusttometals_ratio)

    def _li_ml_dustmass(field,data):
        return (data.ds.arr(data.ds.parameters['li_ml_dustmass'].value,'code_mass'))


    def _dustmass_rr(field,data):
        #hard coded values from remy-ruyer table 1
        a = 2.21
        alpha = 2.02
        x_sun = 8.69

        x = 12.+np.log10(data["Gas","Metals"]/cfg.par.solar * 10.**(x_sun-12.) )
        y = a + alpha*(x_sun-np.asarray(x))
        gas_to_dust_ratio = 10.**(y)
        dust_to_gas_ratio = 1./gas_to_dust_ratio
        return dust_to_gas_ratio * data["Gas","Mass"]

    def _dustmass_li_bestfit(field,data):
        log_dust_to_gas_ratio = (2.445*np.log10(data["Gas","Metals"]/cfg.par.solar))-(2.029)
        dust_to_gas_ratio = 10.**(log_dust_to_gas_ratio)
        return dust_to_gas_ratio * data["Gas","Mass"]


    def _dustsmoothedmasses(field, data):
        if float(yt.__version__[0:3]) >= 4:
            return (data.ds.parameters['octree'][("dust","mass")])
            
        else:
            return (data.ds.arr(data[("Gas", "dust_mass")].value, 'code_mass'))


    def _li_ml_dustsmoothedmasses(field,data):
        return (data.ds.parameters['octree'][("Gas", "li_ml_dustmass")])


    def _metaldens(field,data):
        return (data["Gas","Density"]*data["Gas","Metals"])

    def _gasfh2(field,data):
        try: return data[("Gas","H2")]
        except: return (data.ds.arr(data[("Gas", "Mass")].value*0, 'dimensionless'))

    def _gassfr(field,data):
        instsfr = data.ds.arr(data[("Gas", "Mass")].value*0, 'Msun/yr')
        # temperature below which SF can happen (K)
        try: tempcut = ds.quan(ds.parameters["dTempMax"],'K')
        except: tempcut = ds.quan(10**3,'K') # default for g14 sims
        denscut = ds.quan(100*1.67*10**-24,'g/cm**3') # density above which SF can happen - default for g14 sims
        # star formation timescale (s)
        try: deltat = ds.quan(ds.parameters["dDeltaStarForm"]*3.15569*10**7,'s')
        except: deltat = ds.quan(10**6*3.15569*10**7,'s') # 1 Myr
        # Star formation efficiency
        try: cstarfac = ds.parameters["dCStar"]
        except: cstarfac = 0.15
        # Mass at which stars form (Msol)
        try: massform = ds.quan(ds.parameters["dInitStarMass"]*ds.parameters["dMsolUnit"],'Msun')
        except: massform = max(data[("newstars", "Mass")])
        SFposs = np.where((data[("Gas","Temperature")]<tempcut) & (data[("Gas","Density")]>denscut))
        try: cstar = cstarfac*data[("Gas","H2")][SFposs] # SF efficiency with molecular hydrogen
        except: cstar = cstarfac*ds.data.arr(np.ones(len(data[("Gas","Mass")][SFposs])),"dimensionless")
        dynt = 1.0/np.sqrt(data[("Gas","Density")][SFposs]*4.0*np.pi*G)
        SFeff = 1.0 - np.e**(-1.0*cstar*deltat/dynt)
        instsfr[SFposs] = (data[("Gas","Mass")][SFposs]/massform)*SFeff
        return instsfr

    if fname != None:
        ds = yt.load(fname)
        ds.index
    
    #if we're in the 4.x branch of yt, load up the octree for smoothing
    if float(yt.__version__[0:3]) >= 4:
        left = np.array([pos[0] for pos in bounding_box])
        right = np.array([pos[1] for pos in bounding_box])
        # Add option for scatter vs gather to parameters_master file?
        octree = ds.octree(left,right,n_ref=cfg.par.n_ref)
        ds.parameters['octree'] = octree
    

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
    ds.add_field(("star","metals"),function=_starmetals,units="code_metallicity",sampling_type='particle')
    ds.add_field(("star","coordinates"),function=_starcoordinates,units='cm',sampling_type='particle')
    ds.add_field(("star","formationtime"),function=_starformationtime,units='Gyr',sampling_type='particle')
    ds.add_field(("stellar","ages"),function=_stellarages,units='Gyr',sampling_type='particle')
    ds.add_field(("star","masses"),function=_starmasses,units='g',sampling_type='particle')
    ds.add_field(("gas","density"),function=_gasdensity,units='g/cm**3',sampling_type='particle')
    ds.add_field(("gas","metals"),function=_gasmetals,units="code_metallicity",sampling_type='particle')
    ds.add_field(("gas","coordinates"),function=_gascoordinates,units='cm',sampling_type='particle')
    if add_smoothed_quantities == True:
        ds.add_field(("gas","smootheddensity"),function=_gassmootheddensity,units='g/cm**3',sampling_type='particle')
        ds.add_field(("gas","smoothedmetals"),function=_gassmoothedmetals,units='code_metallicity',sampling_type='particle')
        ds.add_field(("gas","smoothedmasses"),function=_gassmoothedmasses,units='g',sampling_type='particle')
    ds.add_field(("gas","masses"),function=_gasmasses,units='g',sampling_type='particle')
    ds.add_field(("gas","fh2"),function=_gasfh2,units='dimensionless',sampling_type='particle')
    ds.add_field(("gas","sfr"),function=_gassfr,units='g/s',sampling_type='particle')
    ds.add_field(("gas","smoothinglength"),function=_gassmoothinglength,units='pc',sampling_type='particle')
    ds.add_field(("metal","dens"),function=_metaldens,units="g/cm**3", sampling_type='particle')
    #only add the disk star fields if there are any disk stars
    if len(ad["diskstars","Mass"]) > 0 :
        ds.add_field(("diskstar","masses"),function=_diskstarmasses,units='g',sampling_type='particle')
        ds.add_field(("diskstar","coordinates"),function=_diskstarcoordinates,units='cm',sampling_type='particle')

    #add the dust mass, but only if we're using the DTM dust mass
    if cfg.par.dust_grid_type == 'dtm':
        ds.add_field(("dust","mass"), function=_dustmass_dtm,units='code_mass',sampling_type='particle')

    if cfg.par.dust_grid_type == 'rr':
        ds.add_field(("dust","mass"),function=_dustmass_rr,  units='code_mass',sampling_type='particle')



    if cfg.par.dust_grid_type == 'li_ml':
        if float(yt.__version__[0:3]) < 4:
            raise KeyError("li_ml is not implemented for tipsy front ends with yt<4. Please set another dust grid type in the parameters or upgrade to yt4")
        #the issue is that we can't smooth the dust onto a grid that
        #is necessary in dust_grid_gen/li_ml_oct because yt3.x
        #requires a field: ('PartType0', 'particle_position').
        else:
        #get the dust to gas ratio
            ad = ds.all_data()
            li_ml_dgr = dgr_ert(ad["Gas","Metals"],ad["gas","sfr"],ad["Gas","Mass"])
            li_ml_dustmass = ((10.**li_ml_dgr)*ad["Gas","Mass"]).in_units('code_mass')
            #this is an icky way to pass this to the function for ds.add_field in the next line. but such is life.
            ds.parameters['li_ml_dustmass'] = li_ml_dustmass
            ds.add_field(("Gas","li_ml_dustmass"),function=_li_ml_dustmass,units='code_mass',sampling_type='particle')
            #just adding a new field that is called 'dustmass' so that we
            #can save it later in analytics.
            ds.add_field(("dust","mass"),function=_li_ml_dustmass,units='code_mass',sampling_type='particle')
            ds.add_deposited_particle_field(("Gas","li_ml_dustmass"),"sum")
            if add_smoothed_quantities == True:
                ds.add_field(("Gas","li_ml_dustsmoothedmasses"), function=_li_ml_dustsmoothedmasses,units='code_mass',sampling_type='particle')


    if cfg.par.dust_grid_type == 'li_bestfit':
        ds.add_field(("dust","mass"),function=_dustmass_li_bestfit,  units='code_mass',sampling_type='particle')

    ad = ds.all_data()
    return ds
