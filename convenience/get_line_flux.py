#example code to get nebular line flux from a pd simulation

#originally written by prerak garg (university of florida), updated a
#bit by desika narayanan

from hyperion.model import ModelOutput
import numpy as np
from scipy.optimize import curve_fit
from scipy import integrate
import astropy.units as u

def func(x, *params):
    return (params[1]*np.exp(-(((x - params[0])/abs(params[2]))**2)*(2.772588724)))+abs(params[3])

def get_flux(lam_cent,left_edge,right_edge,lam, spec, fit_level):
    # Calculates line luminosity by fitting a Gaussian. 
    # Needs wavelength in microns and SED in ergs/s/microns
    # fit_level is the minimum Rsquare value needed. It detemines the goodness of fit (the higher it is the better the fit) 
    # Returns -1.0 if there is an error or if Rsuqare of the fit is less than fit_level
    
    left = left_edge
    right = right_edge
    lam_peak = lam_cent


    #left = 0.6556
    #right = 0.6576
    #lam_peak = 0.65645131
    
    index_left = np.abs(lam - left).argmin()
    index_right = np.abs(lam - right).argmin()

    spec_clip = spec[int(index_left):int(index_right) + 1]
    lam_clip = lam[int(index_left):int(index_right) + 1]
   
    peak_index = np.abs(lam_clip - lam_peak).argmin()

    cont = (spec_clip[0] + spec_clip[-1]) / 2.
    spec_clip = spec_clip/cont
    guess = [lam_peak, spec_clip[peak_index], abs((lam_clip[-1] - lam_clip[0])), 1.]

    try:
        x = np.linspace(lam_clip[0], lam_clip[-1], 1000)
        popt, pcov = curve_fit(func, lam_clip, spec_clip, p0=guess)
        fit = func(x, popt[0], popt[1], popt[2], popt[3])
        fit = fit - abs(popt[3])
        line_flux = integrate.simps(fit,x)*cont
        y = func(lam_clip, popt[0], popt[1], popt[2], popt[3])

        ss_res = np.sum((spec_clip - y) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r2 = 1 - (ss_res / ss_tot)
        if (r2 <= fit_level):
            print ("**BAD FIT**")
            lune_flux = -1.0
        
        return line_flux

    except Exception as error:
        print (str(error))
        return -1.0



# Main Program Code

#get the SED and get the units right
fname = '/blue/narayanan/prerakgarg/pd_runs/m25n512_jp/snap305_mono/snap305.galaxy100.rtout.sed'
m = ModelOutput(fname)
wav,nufnu = m.get_sed(inclination='all',aperture=-1)
nufnu *= u.erg/u.s
wav *= u.micron
flam = (nufnu/wav).to(u.erg/u.s/u.micron)
flam = flam[0]

lam_cent = 0.6564 #micron of central wavelength
left_edge = lam_cent*0.999 #wavelength of left edge of line (guessed)
right_edge = lam_cent*1.001 #wavelength of right edge of line (guessed)


line_flux = get_flux(lam_cent,left_edge,right_edge,wav.value, flam.value, 0.98)
print ("line flux :"+ str(line_flux) +" ergs/s")

