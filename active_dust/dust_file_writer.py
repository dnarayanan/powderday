#Based on routines originally written by Qi Li, University of Florida

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

import astropy.units as u
import astropy.constants as constants
from hyperion.dust import IsotropicDust
import os


def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx

def Q_unit_transform(Q):
	# wavelength [micron] -> 1/wavelength [1/micron]
	# linear -> log10 
	Q[:,0] = np.log10(1/Q[:,0])
	Q[:,1] = np.log10(Q[:,1])
	return Q

def Q_interp_wl(wlen,Q0):
	Q = interpolate.interp1d(Q0[:,0],Q0[:,1],fill_value='extrapolate')
	Q_wl = Q(wlen)
	return Q_wl

def Qext_tab_load(fin_xtab = 'grain_size.txt',fin_gra='Gra_Optical/Gra_LD93_',fin_sil='Sil_Optical/Sil_LD93_'):
	# Return
	#  xtab, tuple Qtab = (Qabs_C, Qabs_Si, Qsca_C, Qsca_Si)
	xtab = np.log10(np.loadtxt(fin_xtab))
	ntab = len(xtab)
	Qabs_C = [[] for i in range(ntab)]
	Qsca_C = [[] for i in range(ntab)]
	Qabs_Si = [[] for i in range(ntab)]
	Qsca_Si = [[] for i in range(ntab)]

	# read Qext tables into arrays
	for i in range(ntab):
		num = str(i).zfill(2)
		Qabs_C[i] = np.loadtxt(fin_gra + num,usecols=(0,1))
		Qabs_C[i] = Q_unit_transform(Qabs_C[i])
		Qabs_Si[i] = np.loadtxt(fin_sil + num,usecols=(0,1))
		Qabs_Si[i] = Q_unit_transform(Qabs_Si[i])
		Qsca_C[i] = np.loadtxt(fin_gra + num,usecols=(0,2))
		Qsca_C[i] = Q_unit_transform(Qsca_C[i])
		Qsca_Si[i] = np.loadtxt(fin_sil + num,usecols=(0,2))
		Qsca_Si[i] = Q_unit_transform(Qsca_Si[i])

	return xtab,(Qabs_C,Qabs_Si,Qsca_C,Qsca_Si)


def Qext_get(x,wlen,cfrac,xtab,Qtab):
	# Input
	#  x: array of log10(grain_radii) of the simulation
	#  wlen: wavelengths of desired extinction curve
	#  cfrac: mass fraction of graphite
	#  xtab: array of log10(grain_radii) of the optical property table
    #  Qtab: tabulated extinction coefficient
	# Return
	#  tuple (Q_absorption, Q_scatter)

	ntab = len(xtab)
	dim = (len(x),len(wlen)) # dim: wavelength x grain size
	wlen = np.log10(1. / wlen)
	print(wlen)

 	# original tables
	Qabs_C,Qabs_Si,Qsca_C,Qsca_Si = Qtab[0],Qtab[1],Qtab[2],Qtab[3]

	# extended in wlen
	Qabs_C_wl = np.zeros((ntab,dim[1]))
	Qsca_C_wl = np.zeros((ntab,dim[1]))
	Qabs_Si_wl = np.zeros((ntab,dim[1]))
	Qsca_Si_wl = np.zeros((ntab,dim[1]))

	# fully extended table (additionally extended in grain size)
	Qabs_C_2D = np.zeros(dim)
	Qsca_C_2D = np.zeros(dim)
	Qabs_Si_2D = np.zeros(dim)
	Qsca_Si_2D = np.zeros(dim)


	# align arrays across grain sizes by linear interpolation
	for i in range(ntab):
		Qabs_C_wl[i,:] = Q_interp_wl(wlen,Qabs_C[i])
		Qabs_Si_wl[i,:] = Q_interp_wl(wlen,Qabs_Si[i])
		Qsca_C_wl[i,:] = Q_interp_wl(wlen,Qsca_C[i])
		Qsca_Si_wl[i,:] = Q_interp_wl(wlen,Qsca_Si[i])

	# bilinear interpolate to obtain the extended table
	for i in range(dim[1]):
		Q = interpolate.interp1d(xtab,Qabs_C_wl[:,i],fill_value='extrapolate')
		Qabs_C_2D[:,i] = Q(x)
		Q = interpolate.interp1d(xtab,Qabs_Si_wl[:,i],fill_value='extrapolate')
		Qabs_Si_2D[:,i] = Q(x)
		Q = interpolate.interp1d(xtab,Qsca_C_wl[:,i],fill_value='extrapolate')
		Qsca_C_2D[:,i] = Q(x)
		Q = interpolate.interp1d(xtab,Qsca_Si_wl[:,i],fill_value='extrapolate')
		Qsca_Si_2D[:,i] = Q(x)


	return (cfrac*10**Qabs_C_2D + (1.0-cfrac)*10**Qabs_Si_2D, cfrac*10**Qsca_C_2D + (1.0-cfrac)*10**Qsca_Si_2D)



def extinction_law(x,dsf,wlen,cfrac,t_Qext,t_Qext_V):
	# Input
	#  x: array of log10(grain_radii) of the simulation
	#  dsf: counts (not mass!) of dust grains within bins with centers x
	#  wlen: wavelengths of desired extinction curve
	#  cfrac: mass fraction of graphite
    #  t_Qext: tuple (Q_absorption, Q_scatter)
	# return
	#   A/Av, R (albedo)

	a = 10**x
	Qabs, Qsca = t_Qext[0], t_Qext[1]
	Qabs_V, Qsca_V = t_Qext_V[0].transpose()[0], t_Qext_V[0].transpose()[0]
	Qext = Qabs + Qsca
	Qext_V = Qabs_V + Qsca_V
		
	A = np.zeros(len(wlen))
	Asca = np.zeros(len(wlen))
	Av = 2.5 * np.log10(np.e) * np.pi * np.sum(Qext_V*a**2*dsf)
	for i in range(len(wlen)):
		Asca[i] = 2.5 * np.log10(np.e) * np.pi * np.sum(Qsca[:,i]*a**2*dsf)
		A[i] = 2.5 * np.log10(np.e) * np.pi * np.sum(Qext[:,i]*a**2*dsf)
	R = Asca / A
	A #/= Av
	print(Av)
	return A, R, Qext



if __name__ == "__main__":

    #grain size range that we're modeling.  we set it up as a linspace
    #array so that it can have two purposes: (a) to serve as an input
    #array for a test run where we compare against a known size
    #distribution function read in from DNSF_example.txt, and then (b)
    #to set the limits for the grain size bins that we're going to
    #create for powderday
    x=np.linspace(-4,0,41)
        
    #avelength array that we're modeling: 0.1-1000 micron
    wlen = 1. / np.logspace(-4,3,201)*u.micron
    ASSUMED_DENSITY_OF_DUST = 2.4*u.g/u.cm**3

    nu = (constants.c/wlen).to(u.Hz)
    
    #load an example dust size function for testing against
    dsf = np.loadtxt('DNSF_example.txt')

    #assumed graphite fraction
    cfrac = 0.54

    #compute quantities for exmaple dsf above
    xtab, Qtab = Qext_tab_load()
    t_Qext = Qext_get(x,wlen.value,cfrac,xtab,Qtab)
    t_Qext_V = Qext_get(x,np.array([0.551]),cfrac,xtab,Qtab)
    Aext, R, dumQext = extinction_law(x,dsf,wlen.value,cfrac,t_Qext,t_Qext_V)
        
    #nowbreak up the size distribution into bins and then
    #scale wiht the loaded up DSF to see if their co-added
    #extinction laws look reasonable or not.
    
    nbins = 25
    
    #array that holds the left edge, right edge, and edges of bins in between.  we use this to set the left and right edge arrays
    edges = np.linspace(np.min(x),np.max(x),nbins+1) 
    grain_size_left_edge_array = edges[0:-1]
    grain_size_right_edge_array = edges[1::]

    #grain_size_left_edge_array = np.linspace(np.min(x),np.max(x),nbins)
    #grain_size_right_edge_array = []
    outfile_filenames = []

    Aext_array = np.zeros([len(wlen.value),nbins])

    #set up an array for storing the weighted fractional contribution
    #each size bin contributes to the final extinction curve. note -
    #this is only done for testing here by comparing these co-added
    #extinction laws back to the original one computed from dsf
    frac = np.zeros(nbins)
    
    #loop through the bins (but not the right most one, hence the -1)
    for counter,i in enumerate(range(nbins)):
        
        grain_sizes_this_bin = np.linspace(grain_size_left_edge_array[i],grain_size_right_edge_array[i],41)#this 41 is an arbitrary choice
        
        #save the right edge of the bin
        #grain_size_right_edge_array.append(grain_size_left_edge_array[i+1])
        

        #used to figure out what bin size of the sample dsf the current size bin we're on corresponds to
        idx = find_nearest(x,grain_size_left_edge_array[i])
        
        #assuming a flat distribution of sizes within each bin
        temp_dsf = np.repeat(1.e59,len(grain_sizes_this_bin))
        
        #set the weighted fraction to scale the binned final dust size distribution by  for testing
        frac[i] = dsf[idx]

        #calculate Qext for the bins
        t_Qext = Qext_get(grain_sizes_this_bin,wlen.value,cfrac,xtab,Qtab)
        t_Qext_V = Qext_get(grain_sizes_this_bin,np.array([0.551]),cfrac,xtab,Qtab)
        
        #compute the stuff!
        temp_Aext,albedo,temp_Qext = extinction_law(grain_sizes_this_bin,temp_dsf,wlen.value,cfrac,t_Qext,t_Qext_V)
        
        #save the extinction laws to an array, just for testing
        Aext_array[:,i] = temp_Aext
        
        
        #calculate kappa
        #kappa (a) = 3 * Qext (a) / ( 4 * a * rho_grain)

        kappa_a_lambda = np.zeros(temp_Qext.shape)
        for i in range(kappa_a_lambda.shape[1]):
            kappa_a_lambda[:,i] = 3.*temp_Qext[:,i]/(4.* ((10.**grain_sizes_this_bin)*u.micron).to(u.cm)*ASSUMED_DENSITY_OF_DUST)
        #kappa_lambda = np.sum(kappa_a_lambda,axis=0)

        #to go from kappa_a_lambda-->kappa_lambda, we need to
        #integrate over a size distribution (as we want the
        #distribution-weighted average over the size range we're
        #modeling here).  this represents a minor inconsistency in
        #these methods as it {\it a priori} assumes a size
        #distribution.  we will assume MRN (dn/da = a^-3.5, so dn = a^-2.5)

        dN = 1400.*(10.**grain_sizes_this_bin)**(-2.5)
        kappa_lambda = np.absolute(np.trapz(kappa_a_lambda,dN,axis=0)/np.sum(dN))

        lam_micron = (constants.c/nu).to(u.micron)

        #----------------------------------
        #create the HDF5 file for powderday
        #----------------------------------


        d = IsotropicDust(nu.value,albedo,kappa_lambda)

        if not os.path.exists('dust_files/'):
            os.makedirs('dust_files/')
        filename = 'dust_files/binned_dust_sizes.'+str(counter)+'.hdf5'
        outfile_filenames.append(filename)

        d.write(filename)




    #----------------------------------
    #save the metadata
    #----------------------------------
    grain_size_right_edge_array = np.asarray(grain_size_right_edge_array)

    x = grain_size_left_edge_array[0:-1]
    y = grain_size_right_edge_array
    z = np.asarray(outfile_filenames)
    #np.savetxt('dust_files/binned_dust_sizes.key',np.transpose([grain_size_left_edge_array[0:-1],grain_size_right_edge_array,np.asarray(outfile_filenames)]))

    np.savez('dust_files/binned_dust_sizes.npz',grain_size_left_edge_array = grain_size_left_edge_array,grain_size_right_edge_array = grain_size_right_edge_array,outfile_filenames = outfile_filenames)



    #----------------------------------
    #Making some plots just for testing
    #----------------------------------
    final_Aext = np.average(Aext_array,axis=1,weights=frac)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.loglog(1/wlen,final_Aext/np.max(final_Aext),label='coadded')
    #ax.loglog(1/wlen,Aext/np.max(Aext),label='original')
    ax.loglog(wlen,final_Aext,label='coadded')
    ax.loglog(wlen,Aext,label='original')

    plt.legend(loc=2)

    fig.savefig('extinction.png',dpi=300)
        
