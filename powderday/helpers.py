import numpy as np
import powderday.config as cfg
from astropy import units as u
# from astropy.modeling.blackbody import blackbody_lambda,blackbody_nu
from astropy.modeling.models import BlackBody
import h5py

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    
    return idx

def get_J_CMB():
    #returns the mean intensity for the CMB integrated over min_lam to
    #max_lam (i.e. returns erg/s/cm**2; the same thing as doing 4*sigma
    #T^4)
    min_lam = (1.*u.angstrom).to(u.micron)
    max_lam = (1*u.cm).to(u.micron)

    wavelengths = np.linspace(min_lam,max_lam,1.e5)

    # flux = blackbody_lambda(wavelengths,cfg.model.TCMB) # blackbody_lambda deprecated in astropy v4.3
    ## equivalent function provided here https://docs.astropy.org/en/v4.1/modeling/blackbody_deprecated.html#blackbody-functions-deprecated
    bb_lam = BlackBody(cfg.model.TCMB * u.K, scale=1.0 * u.erg / (u.cm ** 2 * u.AA * u.s * u.sr))
    flux = bb_lam(wavelengths)

    J = np.trapz(flux,wavelengths).to(u.erg/u.s/u.cm**2/u.sr)
    solid_angle = 4.*np.pi*u.sr
    J = J*solid_angle
    return J




def energy_density_absorbed_by_CMB():
    extinction_file = cfg.par.dustdir+cfg.par.dustfile

    mw_df = h5py.File(extinction_file,'r')
    mw_o = mw_df['optical_properties']
    mw_df_nu = mw_o['nu']*u.Hz
    mw_df_chi = mw_o['chi']*u.cm**2/u.g
 

    # b_nu = blackbody_nu(mw_df_nu,cfg.model.TCMB) # blackbody_nu deprecated in astropy v4.3
    ## equivalent function provided here https://docs.astropy.org/en/v4.1/modeling/blackbody_deprecated.html#blackbody-functions-deprecated
    bb_nu = BlackBody(cfg.model.TCMB * u.K)
    b_nu = bb_nu(mw_df_nu)
    
    #energy_density_absorbed = 4pi int b_nu * kappa_nu d_nu since b_nu
    #has units erg/s/cm^2/Hz/str and kappa_nu has units cm^2/g.  this results in units erg/s/g
    steradians = 4*np.pi*u.sr
    energy_density_absorbed = (steradians*np.trapz( (b_nu*mw_df_chi),mw_df_nu)).to(u.erg/u.s/u.g)
    return energy_density_absorbed
