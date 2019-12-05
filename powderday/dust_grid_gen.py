# -*- coding: utf-8 -*-


from __future__ import print_function
import numpy as np
import powderday.config as cfg
from powderday.mlt.dgr_extrarandomtree_part import dgr_ert


def manual(pf,refined):
    wTrue = np.where(np.array(refined) == True)[0]
    wFalse = np.where(np.array(refined) == False)[0]
    
    ad = pf.all_data()
    density_smoothed = ad["gassmootheddensity"]
    metallicity_smoothed = ad["gassmoothedmetals"]
    masses_smoothed = ad["gassmoothedmasses"]

    
    try:
        smoothed_dust_masses = ad[('dustsmoothedmasses')]
    except:
        raise KeyError('Dust mass information not present in this snapshot. Please set another dust grid type in the parameters.')
    dust_to_gas_ratio = smoothed_dust_masses.in_units('g')/masses_smoothed
    #masses_smoothed can be 0 at some places; this will make dtg nan
    #out even though it would eventually get multiplied to 0 when we
    #multiply by density smoothed.  So, to be stable we nan_to_num it.
    dust_to_gas_ratio = np.nan_to_num(dust_to_gas_ratio)
    dust_smoothed = np.zeros(len(refined))
    dust_smoothed[wFalse]  = dust_to_gas_ratio * density_smoothed
    return dust_smoothed

def dtm_grid(pf,refined):
    wTrue = np.where(np.array(refined) == True)[0]
    wFalse = np.where(np.array(refined) == False)[0]


    ad = pf.all_data()
    density_smoothed = ad["gassmootheddensity"]
    metallicity_smoothed = ad["gassmoothedmetals"]
    masses_smoothed = ad["gassmoothedmasses"]
    
    dust_smoothed = np.zeros(len(refined))
    
    print ('[grid_construction/dust_grid_gen/dtm_grid: ] len(wFalse) = ',len(wFalse))
    print ('[grid_construction/dust_grid_gen/dtm_grid: ] len(metallicity_smoothed) = ',len(metallicity_smoothed))

    dust_smoothed[wFalse] = metallicity_smoothed * density_smoothed * cfg.par.dusttometals_ratio
    return dust_smoothed



def remy_ruyer(pf,refined):
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


    ad = pf.all_data()
    density_smoothed = ad["gassmootheddensity"]
    metallicity_smoothed = ad["gassmoothedmetals"]
    masses_smoothed = ad["gassmoothedmasses"]
 
    #anywhere the smoothing finds a cell with zero metallicity, set
    #this to some very low value
    metallicity_smoothed[metallicity_smoothed == 0] = 1.e-10

    x = 12.+np.log10(metallicity_smoothed/cfg.par.solar * 10.**(x_sun-12.) )
    
    y = a + alpha*(x_sun-np.asarray(x))

    gas_to_dust_ratio = 10.**(y)
    dust_to_gas_ratio = 1./gas_to_dust_ratio

    dust_smoothed = np.zeros(len(refined))
    dust_smoothed[wFalse]  = dust_to_gas_ratio * density_smoothed
    
    return dust_smoothed


def li_bestfit(pf,refined):
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

    ad = pf.all_data()
    density_smoothed = ad["gassmootheddensity"]
    metallicity_smoothed = ad["gassmoothedmetals"]
    masses_smoothed = ad["gassmoothedmasses"]

    #anywhere the smoothing finds a cell with zero metallicity, set
    #this to some very low value
    metallicity_smoothed[metallicity_smoothed == 0] = 1.e-10

    log_dust_to_gas_ratio = (2.445*np.log10(metallicity_smoothed/cfg.par.solar))-(2.029)
    dust_to_gas_ratio = 10.**(log_dust_to_gas_ratio)

    dust_smoothed = np.zeros(len(refined))
    dust_smoothed[wFalse]  = dust_to_gas_ratio * density_smoothed

    return dust_smoothed


'''
def li_ml(pf,refined):

    wTrue = np.where(np.array(refined) == True)[0]
    wFalse = np.where(np.array(refined) == False)[0]

    ad = pf.all_data()
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

def li_ml(pf,refined):
    
    wTrue = np.where(np.array(refined) == True)[0]
    wFalse = np.where(np.array(refined) == False)[0]

    ad = pf.all_data()
    density_smoothed = ad["gassmootheddensity"]
    metallicity_smoothed = ad["gassmoothedmetals"]
    masses_smoothed = ad["gassmoothedmasses"]


    try:
        smoothed_dust_masses = ad[('li_ml_dustsmoothedmasses')]
    except:
        raise KeyError('Li, Narayanan & Dave mass information not present in this snapshot. Please set another dust grid type in the parameters.')
    dust_to_gas_ratio = smoothed_dust_masses.in_units('g')/masses_smoothed
    #masses_smoothed can be 0 at some places; this will make dtg nan
    #out even though it would eventually get multiplied to 0 when we
    #multiply by density smoothed.  So, to be stable we nan_to_num it.
    dust_to_gas_ratio = np.nan_to_num(dust_to_gas_ratio)
    dust_smoothed = np.zeros(len(refined))
    dust_smoothed[wFalse]  = dust_to_gas_ratio * density_smoothed
    return dust_smoothed
