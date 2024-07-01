# -*- coding: utf-8 -*-


from __future__ import print_function
import numpy as np
import powderday.config as cfg

import pdb,sys

def manual_oct(reg,refined):
    wTrue = np.where(np.array(refined) == True)[0]
    wFalse = np.where(np.array(refined) == False)[0]
    
    
    density_smoothed = reg["gas","smootheddensity"]
    metallicity_smoothed = reg["gas","smoothedmetals"]
    masses_smoothed = reg["gas","smoothedmasses"]

    
    try:
        smoothed_dust_masses = reg[('dust','smoothedmasses')]
    except:
        raise KeyError('Dust mass information not present in this snapshot. Please set another dust grid type in the parameters.')
    dust_to_gas_ratio = smoothed_dust_masses.in_units('g')/masses_smoothed
    #masses_smoothed can be 0 at some places; this will make dtg nan
    #out even though it would eventually get multiplied to 0 when we
    #multiply by density smoothed.  So, to be stable we nan_to_num it.
    
    #print('[dust_grid_gen/manual_oct:] gas smoothed total :',np.sum(masses_smoothed.to('Msun')))
    
    vol_fact = reg.parameters['octree'][('index','dx')][~refined].to('kpc')*reg.parameters['octree'][('index','dy')][~refined].to('kpc')*reg.parameters['octree'][('index','dz')][~refined].to('kpc')

    # correct for dtg issues
    # yt does not necessarily have the smoothed dust and gas masses
    # offset from the true values by the same amount, i.e.
    # total smoothed gas mass/total particle gas mass not =
    # total smoothed dust mass/total particle dust mass
    cor_fact = (np.sum(masses_smoothed.to('Msun'))/np.sum(reg['gas','mass'].to('Msun')))/(np.sum(smoothed_dust_masses.to('Msun'))/np.sum(reg['dust','mass'].to('Msun')))
    # correct for total gas normalization
    cor_fact = cor_fact * np.sum(reg['gas','mass'].to('Msun'))/np.sum(density_smoothed.to('Msun/kpc**3')*vol_fact)

    dust_to_gas_ratio = np.nan_to_num(dust_to_gas_ratio)
    dust_smoothed = np.zeros(len(refined))
    dust_smoothed[wFalse]  = dust_to_gas_ratio * density_smoothed * cor_fact
    return dust_smoothed

def dtm_grid_oct(reg,refined):
    wTrue = np.where(np.array(refined) == True)[0]
    wFalse = np.where(np.array(refined) == False)[0]


    density_smoothed = reg["gas","smootheddensity"]
    metallicity_smoothed = reg["gas","smoothedmetals"]
    masses_smoothed = reg["gas","smoothedmasses"]

    gas_mass_particle = reg['gas','mass'].to('Msun')
    # correct for total gas normalization
    vol_fact = reg.parameters['octree'][('index','dx')][~refined].to('kpc')*reg.parameters['octree'][('index','dy')][~refined].to('kpc')*reg.parameters['octree'][('index','dz')][~refined].to('kpc')
    cor_fact = np.sum(gas_mass_particle)/np.sum(density_smoothed.to('Msun/kpc**3')*vol_fact)

    dust_smoothed = np.zeros(len(refined))
    
    print ('[grid_construction/dust_grid_gen/dtm_grid: ] len(wFalse) = ',len(wFalse))
    print ('[grid_construction/dust_grid_gen/dtm_grid: ] len(metallicity_smoothed) = ',len(metallicity_smoothed))

    dust_smoothed[wFalse] = metallicity_smoothed * density_smoothed * cfg.par.dusttometals_ratio * cor_fact
    return dust_smoothed



def remy_ruyer_oct(reg,refined):
    #remy ruyer 2014 A&A 563, A31 -- here, we take the Xco_Z
    #power-law, slope free parameterization to define the dust-to-gas
    #ratio, and hece the dust density

    #y = log10(G/D) = a+alpha(x_solar - x) 
    
    #x_solar = 12 + log10(O/H) = 8.69
    
    #hard coded values from remy-ruyer table 1
    a = 2.21
    alpha = 2.02
    x_sun = 8.69


    wTrue = np.where(np.array(refined) == True)[0]
    wFalse = np.where(np.array(refined) == False)[0]



    density_smoothed = reg["gas","smootheddensity"]
    metallicity_smoothed = reg["gas","smoothedmetals"]
    masses_smoothed = reg["gas","smoothedmasses"]
    gas_mass_particle = reg['gas','mass'].to('Msun')
    # correct for total gas normalization
    vol_fact = reg.parameters['octree'][('index','dx')][~refined].to('kpc')*reg.parameters['octree'][('index','dy')][~refined].to('kpc')*reg.parameters['octree'][('index','dz')][~refined].to('kpc')
    cor_fact = np.sum(gas_mass_particle)/np.sum(density_smoothed.to('Msun/kpc**3')*vol_fact)

    #anywhere the smoothing finds a cell with zero metallicity, set
    #this to some very low value
    metallicity_smoothed[metallicity_smoothed == 0] = 1.e-10

    x = 12.+np.log10(metallicity_smoothed/cfg.par.solar * 10.**(x_sun-12.) )
    
    y = a + alpha*(x_sun-np.asarray(x))

    gas_to_dust_ratio = 10.**(y)
    dust_to_gas_ratio = 1./gas_to_dust_ratio

    dust_smoothed = np.zeros(len(refined))
    dust_smoothed[wFalse]  = dust_to_gas_ratio * density_smoothed * cor_fact
    
    return dust_smoothed


def li_bestfit_oct(reg,refined):
    #li, narayanan & dave, 2019, arXiv/1906.09277.  here, we take the
    #results of their equation 12 which relates the dust to gas ratio
    #as a function of metallicity from the simba cosmological
    #simulations that include an on-the-fly prescription for dust
    #formation, growth, and destruction.  has a slightly steeper DGR
    #vs Z slope than the remy-ruyer relation.  note, this best fit
    #does not include passive galaxies which lie significantly off of
    #this (and the remy-ruyer) relation.

    
    wTrue = np.where(np.array(refined) == True)[0]
    wFalse = np.where(np.array(refined) == False)[0]


    density_smoothed = reg["gas","smootheddensity"]
    metallicity_smoothed = reg["gas","smoothedmetals"]
    masses_smoothed = reg["gas","smoothedmasses"]

    gas_mass_particle = reg['gas','mass'].to('Msun')
    # correct for total gas normalization
    vol_fact = reg.parameters['octree'][('index','dx')][~refined].to('kpc')*reg.parameters['octree'][('index','dy')][~refined].to('kpc')*reg.parameters['octree'][('index','dz')][~refined].to('kpc')
    cor_fact = np.sum(gas_mass_particle)/np.sum(density_smoothed.to('Msun/kpc**3')*vol_fact)

    #anywhere the smoothing finds a cell with zero metallicity, set
    #this to some very low value
    metallicity_smoothed[metallicity_smoothed == 0] = 1.e-10

    log_dust_to_gas_ratio = (2.445*np.log10(metallicity_smoothed/cfg.par.solar))-(2.029)
    dust_to_gas_ratio = 10.**(log_dust_to_gas_ratio)

    dust_smoothed = np.zeros(len(refined))
    dust_smoothed[wFalse]  = dust_to_gas_ratio * density_smoothed * cor_fact

    return dust_smoothed


'''
def li_ml(ds,refined):

    wTrue = np.where(np.array(refined) == True)[0]
    wFalse = np.where(np.array(refined) == False)[0]

    ad = ds.all_data()
    density_smoothed = ad["gassmootheddensity"]
    metallicity_smoothed = ad["gassmoothedmetals"]
    masses_smoothed = ad["gassmoothedmasses"]
    star_masses_smoothed = ad["starsmoothedmasses"]
    

    dgr = dgr_ert(metallicity_smoothed,star_masses_smoothed,masses_smoothed)
    dust_to_gas_ratio = 10.**(dgr)

    dust_smoothed = np.zeros(len(refined))
    dust_smoothed[wFalse] = dust_to_gas_ratio * density_smoothed

    return dust_smoothed
'''

def li_ml_oct(reg,refined):
    
    wTrue = np.where(np.array(refined) == True)[0]
    wFalse = np.where(np.array(refined) == False)[0]

    
    density_smoothed = reg["gas","smootheddensity"]
    metallicity_smoothed = reg["gas","smoothedmetals"]
    masses_smoothed = reg["gas","smoothedmasses"]


    try:
        smoothed_dust_masses = reg[('li_ml_dustsmoothedmasses')]
    except:
        raise KeyError('Li, Narayanan & Dave mass information not present in this snapshot. Please set another dust grid type in the parameters.')
    dust_to_gas_ratio = smoothed_dust_masses.in_units('g')/masses_smoothed
    #masses_smoothed can be 0 at some places; this will make dtg nan
    #out even though it would eventually get multiplied to 0 when we
    #multiply by density smoothed.  So, to be stable we nan_to_num it.
    
    
    gas_mass_particle = reg['gas','mass'].to('Msun')
    # correct for total gas normalization
    vol_fact = reg.parameters['octree'][('index','dx')][~refined].to('kpc')*reg.parameters['octree'][('index','dy')][~refined].to('kpc')*reg.parameters['octree'][('index','dz')][~refined].to('kpc')
    cor_fact = np.sum(gas_mass_particle)/np.sum(density_smoothed.to('Msun/kpc**3')*vol_fact)
    # correct for dtg issues                                                                                                                                    # yt does not necessarily have the smoothed dust and gas masses                                                                                             # offset from the true values by the same amount, i.e.                                                                                                      # total smoothed gas mass/total particle gas mass not =
    # total smoothed dust mass/total particle dust mass  
    cor_fact = cor_fact*(np.sum(masses_smoothed.to('Msun'))/np.sum(reg['gas','mass'].to('Msun')))/(np.sum(smoothed_dust_masses.to('Msun'))/np.sum(reg['dust','mass'].to('Msun')))
    
    dust_to_gas_ratio = np.nan_to_num(dust_to_gas_ratio)
    dust_smoothed = np.zeros(len(refined))
    dust_smoothed[wFalse]  = dust_to_gas_ratio * density_smoothed * cor_fact
    return dust_smoothed



def dtm_particle_mesh(reg):
    #calculates the dust based on the DTM ratio for either particles
    #directly (i.e. arepo quantities) or a mesh (i.e AMR simulations)

    metaldens = reg["metal","dens"]
    dustdens = (metaldens*cfg.par.dusttometals_ratio).to('g/cm**3').value

    return dustdens


def manual_particle_mesh(reg):

    #calculates the dust based on the DTM ratio for either particles
    #directly (i.e. arepo quantities) or a mesh (i.e AMR simulations)
    if cfg.par.otf_extinction == False:
        if ('PartType0','DustDensity') in reg.ds.derived_field_list:
            dustdens = reg.ds.arr(reg["PartType0","DustDensity"].value,'code_mass/code_length**3')
            dustdens = dustdens.to('g/cm**3').value
        else:
            raise ValueError('It looks like we cant find PartType0,DustDensity in your Arepo simulations. Please try another choice amongst [dtm, rr, li_bestfit, li_ml].  Alternatively, edit [dust_grid_gen/manual_particle_mesh] to change the value of the field assigned to dustdens')
    else:
        if ('PartType3','Dust_DustDensity') in reg.ds.derived_field_list:
            dustdens = reg.ds.arr(reg["PartType3","Dust_DustDensity"].value,'code_mass/code_length**3')
            dustdens = dustdens.to('g/cm**3').value
        else:
            raise ValueError('It looks like we cant find PartType3,Dust_DustDensity in your Arepo simulations. Please try another choice amongst [dtm, rr, li_bestfit, li_ml].  Alternatively, edit [dust_grid_gen/manual_particle_mesh] to change the value of the field assigned to dustdens')


    return dustdens




def remy_ruyer_particle_mesh(reg):
    #remy ruyer 2014 A&A 563, A31 -- here, we take the Xco_Z
    #power-law, slope free parameterization to define the dust-to-gas
    #ratio, and hece the dust density

    #y = log10(G/D) = a+alpha(x_solar - x)

    #x_solar = 12 + log10(O/H) = 8.69

    #hard coded values from remy-ruyer table 1
    a = 2.21
    alpha = 2.02
    x_sun = 8.69


    density = reg["gas","density"]
    metallicity=reg["gas","metals"]
    masses = reg["gas","masses"]

    #anywhere the smoothing finds a cell with zero metallicity, set
    #this to some very low value
    metallicity[metallicity == 0] = 1.e-10

    x = 12.+np.log10(metallicity/cfg.par.solar * 10.**(x_sun-12.) )

    y = a + alpha*(x_sun-np.asarray(x))

    gas_to_dust_ratio = 10.**(y)
    dust_to_gas_ratio = 1./gas_to_dust_ratio

    dustdens = dust_to_gas_ratio * density

    return dustdens



def li_bestfit_particle_mesh(reg):
    #li, narayanan & dave, 2019, arXiv/1906.09277.  here, we take the
    #results of their equation 12 which relates the dust to gas ratio
    #as a function of metallicity from the simba cosmological
    #simulations that include an on-the-fly prescription for dust
    #formation, growth, and destruction.  has a slightly steeper DGR
    #vs Z slope than the remy-ruyer relation.  note, this best fit
    #does not include passive galaxies which lie significantly off of
    #this (and the remy-ruyer) relation.

    density = reg["gas","density"]
    metallicity=reg["gas","metals"]
    masses = reg["gas","masses"]

    #anywhere the smoothing finds a cell with zero metallicity, set
    #this to some very low value
    metallicity[metallicity == 0] = 1.e-10
    
    log_dust_to_gas_ratio = (2.445*np.log10(metallicity/cfg.par.solar))-(2.029)
    dust_to_gas_ratio = 10.**(log_dust_to_gas_ratio)
    
    dustdens = dust_to_gas_ratio * density

    return dustdens


def li_ml_particle_mesh(reg):

    density = reg["gas","density"]
    metallicity=reg["gas","metals"]
    masses = reg["gas","masses"]


    try:
        dust_masses = reg[('li_ml_dustmass')]
    except:
        raise KeyError('Li, Narayanan & Dave mass information not present in this snapshot. Please set another dust grid type in the parameters.')
    dust_to_gas_ratio = dust_masses.in_units('g')/masses
    #masses_ can be 0 at some places; this will make dtg nan
    #out even though it would eventually get multiplied to 0 when we
    #multiply by density smoothed.  So, to be stable we nan_to_num it.
    dust_to_gas_ratio = np.nan_to_num(dust_to_gas_ratio)
    dustdens = dust_to_gas_ratio * density

    return dustdens



#For all the AMR dust definitions, we call functions slightly
#differently.  We set up a new field for ds1 (the zoomed dataset) in
#which we define the dust via whatever physics we are aassuming
#(i.e. dust to metals; Li ML etc.).  then this field is what gets
#added to the AMRGrid via the .from_yt method in the enzo_tributary

def dtm_amr(ds):

    
    def _dust_density_dtm_amr(field, data):
        return data[('gas', 'metal_density')].in_units("g/cm**3")*cfg.par.dusttometals_ratio
    ds.add_field(('gas', 'dust_density'), function=_dust_density_dtm_amr, units = 'g/cm**3', sampling_type='cell')


def remy_ruyer_amr(ds):

    #remy ruyer 2014 A&A 563, A31 -- here, we take the Xco_Z
    #power-law, slope free parameterization to define the dust-to-gas
    #ratio, and hece the dust density

    #y = log10(G/D) = a+alpha(x_solar - x)

    #x_solar = 12 + log10(O/H) = 8.69

    def _dust_density_rr_amr(field,data):

        #hard coded values from remy-ruyer table 1
        a = 2.21
        alpha = 2.02
        x_sun = 8.69


        density = data["gas","density"]
        metallicity = data["gas","metals"]
        masses = data["gasmasses"]

        #anywhere the smoothing finds a cell with zero metallicity, set
        #this to some very low value
        metallicity[metallicity==0] = 1.e-10

        x = 12.+np.log10(metallicity/cfg.par.solar * 10.**(x_sun-12.) )
        
        y = a + alpha*(x_sun-np.asarray(x))
        
        gas_to_dust_ratio = 10.**(y)
        dust_to_gas_ratio = 1./gas_to_dust_ratio
        
        dust_density = (dust_to_gas_ratio*density).in_units('g/cm**3')

        return dust_density

    ds.add_field(('gas', 'dust_density'), function=_dust_density_rr_amr, units = 'g/cm**3', sampling_type='cell')
    
    

def li_bestfit_amr(ds):
    
    def _dust_density_li_bestfit_amr(field,data):
        
        density = data["gas","density"]
        metallicity=data["gas","metals"]
        masses = data["gas","masses"]

        metallicity[metallicity == 0] = 1.e-10

        log_dust_to_gas_ratio = (2.445*np.log10(metallicity/cfg.par.solar))-(2.029)
        dust_to_gas_ratio = 10.**(log_dust_to_gas_ratio)

        dust_density = (dust_to_gas_ratio*density).in_units('g/cm**3')
    
        return dust_density

    ds.add_field(('gas', 'dust_density'), function=_dust_density_li_bestfit_amr, units = 'g/cm**3', sampling_type='cell')

def li_ml_amr(ds):

    raise KeyError('The Li, Narayanan & Dave 2019 Dust model is not currently implemented for Enzo simulations.  Please email desika.narayanan@gmail.com and bug him to put this in.')


