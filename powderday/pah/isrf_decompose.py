import numpy as np
import os,h5py,pdb
from astropy import units as u
from astropy import constants as const
#DEBUG - CHANGE TO THIS WHEN WE USE IN ACTUAL POWDERDAY CODE from powderday.pah.pah_file_read import read_draine_file
from pah_file_read import read_draine_file
from scipy.interpolate import interp1d,interp2d


import pdb

#DEBUG eventually change to be part of cfg.master
draine_data_dir = '/blue/narayanan/desika.narayanan/powderday_files/PAHs/dataverse_files/' 

#get a GSD 
#DEBUG DEBUG DEBUG - we'll want this to be an input from a
#function call eventually.  we'll want this function to read in both
#the sizes of the GSD and the dN.  this is just for testing for now.
data = np.load('extinction_m12n256_MW.npz')
gsd = data['a3N_loga'][0,:]
simulation_sizes = 10.**(data['loga'])*u.micron

#get the wavelengths of the simulation ISRF DEBUG DEBUG DEBUG this
#will eventually, like the above GSD, get read in via a function call.
#this is just for testing.
f = h5py.File('/blue/narayanan/desika.narayanan/pd_runs/powderday_testing/tests/SKIRT/arepo_idealized_extinction_lowres/pd_skirt_comparison.027.otf_dtm.rtout.sed')
dset = f['iteration_00007']
simulation_isrf_nu = dset['ISRF_frequency_bins'][:] * u.Hz
simulation_isrf_lam = (const.c/simulation_isrf_nu).to(u.micron)
simulation_specific_energy_sum = dset['specific_energy_sum_nu']*u.erg/u.s/u.g #is [n_nu, n_dust, n_cells] big  


#get the simulation_isrf in units of erg/s DEBUG DEBUG DEBUG
#eventually we'll want to feed in the dust mass per cell from
#powderday itself, instead of this janky reading in npz files
data = np.load('/blue/narayanan/desika.narayanan/pd_runs/powderday_testing/tests/SKIRT/arepo_idealized_extinction_lowres/grid_physical_properties.027_galaxy0.npz')
grid_dust_masses = data['grid_dustmass']*u.Msun
simulation_specific_energy_sum *= grid_dust_masses.cgs #now in erg/s
#simulation_specific_energy_sum_lam = np.flip(simulation_specific_energy_sum,axis=0) #just ordering in terms of the wavelength grid

#first list all the directories we're going to go into to build the table

draine_directories = []
print('powderday/pah/isrf_decompose]: building a ISRF table from the following directories: ')
for it in os.scandir(draine_data_dir):
    if it.is_dir():
        print(it.path)
        draine_directories.append(it.path)

isrf_files = []
for directory in draine_directories:
    for file in os.listdir(directory):
        if file.startswith("isrf"):
            isrf_files.append(file)

'''#note - this bit isn't formally needed, and even still, it reads
in 2x iout files (compared to the isrf files) since there's a Cabs for
the cation state of PAHs, and one for the neutral


iout_U0_files = [] #just saving the U=0 files since we just need them to grab C_abs
for directory in draine_directories:
    for file in os.listdir(directory):
        if file.startswith("iout_graD") and file.endswith("_0.00"):
            print(directory,file)
            iout_U0_files.append(file)

'''
        


#get the length of a basis ISRF vector
data = np.loadtxt(draine_directories[0]+'/'+isrf_files[0],skiprows=7,usecols=(0,1))
nlam = len(data[:,0])
basis_isrf_vectors = np.zeros([len(draine_directories),nlam])

for counter, (directory,file) in enumerate(zip(draine_directories,isrf_files)):
    data = np.loadtxt(directory+'/'+file,skiprows=7,usecols=(0,1))
    basis_isrf_vectors[counter] = data[:,1]

#add the units (as listed in the files)
basis_isrf_vectors *= u.erg/u.cm**3


#now follow equation 2 of Draine et al. 2021, ApJ, 917, D to get the
#basis ISRF vectors in units of erg/s
#h_ref = int(dnu * basis_isrF * c * C_abs(nu) 
#to do this - we need to read in the C_abs(nu) for every grain size 

#get the Cabs vectors - these are going to be
#(nwavelengths,n_grain_sizes).  These are the same for every radiation
#field (and every logU) so we only need to read in a single file.

for file in os.listdir(draine_directories[0]):
    if file.startswith("iout_graD") and file.endswith("_0.00") and "ib" in file: Cabsfile_cation = draine_directories[0]+'/'+file
    if file.startswith("iout_graD") and file.endswith("_0.00") and "nb" in file: Cabsfile_neutral = draine_directories[0]+'/'+file
PAH_list_cation = read_draine_file(Cabsfile_cation)    
PAH_list_neutral = read_draine_file(Cabsfile_neutral)

#i'm sure there's some awesome pythonic cool kids way to code this.
#but i'm an IDL coder at heart.
Cabs_cation = np.zeros([len(PAH_list_cation[0].size_list),len(PAH_list_cation[0].lam)])
for i in range(len(PAH_list_cation)):
    #this works bc PAH_list is an n_sizes long list of PAH objects
    Cabs_cation[i,:] = PAH_list_cation[i].cabs

Cabs_neutral = np.zeros([len(PAH_list_neutral[0].size_list),len(PAH_list_neutral[0].lam)])
for i in range(len(PAH_list_neutral)):
    #this works bc PAH_list is an n_sizes long list of PAH objects
    Cabs_neutral[i,:] = PAH_list_neutral[i].cabs



#We now have Cabs (for neutrals and cations) in terms of (size,nu),
#and we want to get these in terms of just (nu). To do this, we need
#to regrid (i.e. downsample) the first dimension so that the GSD is at
#the same sizes as the simulation that was run.  Then we can just
#multiply Cabs(size_regridded,nu) by dN, where dN is the number of
#grains at a given size.

# we can just take the very first size of the the cation PAH_list - they should all be the same.
draine_sizes = PAH_list_cation[0].size_list*u.cm 
simulation_sizes = simulation_sizes.to(u.cm)
draine_lam = PAH_list_cation[0].lam*u.micron

'''
#Cabs is in dimensions of draine_sizes,draine_lam. 
#resampling Cabs to a shape of simulation_sizes,draine_lam -- this
#will allow us to deal with the fact that Cabs is given at individual
#grain sizes (by multiplying by simulation GSD), while also keeping the draine_lam axis available for integration later.
'''

f_2d_interp_cation = interp2d(draine_sizes.value,draine_lam.value,Cabs_cation.T,kind = 'cubic') 
f_2d_interp_neutral = interp2d(draine_sizes.value,draine_lam.value,Cabs_neutral.T,kind = 'cubic')
#Cabs_cation_regrid = f_2d_interp_cation(simulation_sizes.value,simulation_isrf_lam.value).T
#Cabs_neutral_regrid = f_2d_interp_neutral(simulation_sizes.value,simulation_isrf_lam.value).T
Cabs_cation_regrid = f_2d_interp_cation(simulation_sizes.value,draine_lam.value).T
Cabs_neutral_regrid = f_2d_interp_neutral(simulation_sizes.value,draine_lam.value).T

#we now have to get Cabs(nu) -- right now we have Cabs(size,nu).
#To do this, we should ensure that the GSD is at the same sizes as the
#draine files (and if not, regrid so that it is), and then multiply
#Cabs(size,nu) * dN, where dN is the number of grains at a given size.

Cabs_cation_regrid = np.dot(gsd,Cabs_cation_regrid)*u.cm**2 # now in terms of just wavelength (at simulation_isrf_lam wavelengths)
Cabs_neutral_regrid = np.dot(gsd,Cabs_neutral_regrid)*u.cm**2

#now to get the ISRF from erg/cm**3-->erg/s just based on dimensional analysis 
basis_isrf_vectors*=Cabs_neutral_regrid*const.c
basis_isrf_vectors = basis_isrf_vectors.to(u.erg/u.s)

#4 now resample the hyperion ISRF to the wavelengths of the Draine
#basis ISRFs so that we can NNLS

f_1d_interp_lam = interp1d(simulation_isrf_lam.to(u.micron).value,simulation_specific_energy_sum.T.value,kind='cubic')
simulation_specific_energy_sum_regrid = f_1d_interp_lam(draine_lam.to(u.micron).value).T

#the regridding can turn some wavelengths where there was 0 emission
#(from too low photon count simulations) to negative, so we zero these
#back out.
simulation_specific_energy_sum_regrid[simulation_specific_energy_sum_regrid < 0] = 0

import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111)
#ax.loglog(draine_lam,basis_isrf_vectors[0,:])
#ax.loglog(draine_lam,simulation_specific_energy_sum_regrid[:,0,0])

kernel_size = 20
kernel = np.ones(kernel_size)/kernel_size
data_convolved= np.zeros(simulation_specific_energy_sum_regrid.shape)
for i in range(simulation_specific_energy_sum.shape[1]): #ndust
    for j in range(simulation_specific_energy_sum.shape[2]): #ncells
        data_convolved[:,i,j] = np.convolve(simulation_specific_energy_sum_regrid[:,i,j],kernel,mode='same')

ax.loglog(draine_lam,data_convolved[:,0,1000],label='cell 0, smoothed')
ax.loglog(draine_lam,np.sum(data_convolved,axis=2)[:,0],label='total for grid, smoothed')
ax.loglog(draine_lam,basis_isrf_vectors[3,:],label='draine #10')

#ax.loglog(draine_lam,simulation_specific_energy_sum_regrid[:,0,0])
    #ax.loglog(simulation_isrf_lam,simulation_specific_energy_sum[:,0,i])
ax.loglog(simulation_isrf_lam,np.sum(simulation_specific_energy_sum,axis=2)[:,0])
#ax.set_xlim([0.1,1.e3])
plt.legend()
fig.savefig('isrf_tests.png',dpi=300)





#5. then nnls!  with nnls, we can then sum the PAH components for each
#cell accordingly.  note - because the ISRF computed from hyperion has
#the infrared component saved, we need to cut off our ISRF for both
#the hyperion model and draine basis functions at some wavelength
#before thermal IR emission gets big, like maybe 10 micron.

