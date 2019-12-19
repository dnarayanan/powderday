# A tool to predict dust-to-gas mass ratio (DGR), using extreme randomized trees trained by Simba simulation data
# input:(1) parameter_list determining which input parameters are active, for now please follow the template param.txt
#       (2) table_input whose columns represent metallicity Z (optional but dangerous if not included, stellar mass Ms (optional), SFR (optional), gas mass Mgas (optional)
# example: python dgr_extrarandomtree.py parameter_list table_input [DGR_output]

#written by Qi Li, University of Florida


from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import numpy as np
#import matplotlib.pyplot as plt
from sklearn.ensemble import ExtraTreesRegressor as ETR
from sklearn.cluster import DBSCAN
import h5py
import powderday.config as cfg
import sys

def filt_gal_prop(X,ipar,labels):
    ### filter out unused galaxy properties and corresponding labels
	mask = []
	print("Disabled input parameters:")
	if ipar[0] == 0: # deactive Z -> Z
		print("-- Gas-phase metallicity")
		mask.append(0)
		labels = filter(lambda a: a != r'$Z$', labels)
	if True: # legacy: always deactivate Ms -> Ms, fgas
		print("-- Stellar mass")
		mask.append(1)
		mask.append(2)
		labels = filter(lambda a: a != r'$M_*$', labels)
		labels = filter(lambda a: a != r'$f_g$', labels)
	if ipar[2] == 0: # deactivate sfr -> tau_sf
		print("-- Star formation rate")
		mask.append(3)
		labels = filter(lambda a: a != r'$\tau_{\rm depletion}$', labels)
	if ipar[3] == 0: # deactive Mg -> fgas, tau_sf
		print("-- gas mass")
		mask.append(2)
		mask.append(3)
		labels = filter(lambda a: a != r'$f_g$', labels)
		labels = filter(lambda a: a != r'$\tau_{\rm depletion}$', labels)
	mask = np.unique(mask)

	return np.delete(X,mask,axis=1),labels


def input_prepare(ipar,tab):
	### calculate Z, Ms, fgas, tau_sf based on input parameters
	if len(tab.shape) == 1: # in case only one row
		NN = 1
		KK = tab.shape[0] # the number of input galaxy properties
		tab = tab.reshape((1,KK))
	else:
		NN = tab.shape[0]
		KK = tab.shape[1]
	if len(ipar[ipar!=0]) != KK:
		print("ERROR: the number of input parameters mismatch the number of active parameters") # it sounds confusing
		return np.zeros(0)
	
	tab_tmp = np.zeros((NN,len(ipar)))
	X = np.zeros((NN,4))
	tab_tmp[:,:] = 1 # avoid RuntimeWarning
	j = 0
	for i in range(len(ipar)):
		if ipar[i] == 1 and i!=1:  #legacy: no Ms for particle-based data
			tab_tmp[:,i] = tab[:,j]
			j += 1
	Z = np.log10(tab_tmp[:,0]/0.0134)
	sfr = tab_tmp[:,2]
	Mg = tab_tmp[:,3]
	tau_sf = Mg/sfr # depletion time 
	
	X[:,0] = Z #legacy: no Ms for particle-based data
	X[:,3] = np.log10(tau_sf)
	
	X,_ = filt_gal_prop(X,ipar,[])

	return X
	
	
def matrix_prepare(ipar,sim='m12n256.hdf5'): # load simulated dataset and prepare the matrix used for regression
	### load simulated dataset
	gal = h5py.File(sim,'r')
	gal = gal['PartType0']

	Z  = np.log10(gal['Metallicity'][:]/cfg.par.solar) # Zsolar # X0  #XC0
	Md = gal['Dust_Masses'][:] # Msolar
	Mg = gal['Masses'][:] - Md
	sfr = gal['StarFormationRate'][:] # Msolar/yr                
	tau_sf = Mg/sfr     # depletion time            # X3
	dgr = Md/Mg                                     # Y   #XC2
	dgr[np.where(dgr<0)] = 0.

	NN = len(Md)
	KK = 4
	labels = [r'$Z$',r'$M_*$',r'$f_g$',r'$\tau_{\rm depletion}$']

		
	### load matrix for regression
	X = np.zeros((NN,KK+2))
	X[:,0] = Z
	X[:,1] = 1.
	X[:,2] = 1.
	X[:,3] = np.log10(tau_sf)	
	X[:,KK] = np.log10(dgr)
	X[:,KK+1] = 1. #sSFR, used for clustering

	### clean up data and prepare vectors needed for fitting
	X = X[~np.isnan(X).any(axis=1)]
	X = X[~np.isinf(X).any(axis=1)]
	X = X[np.lexsort(np.rot90(X))] # sort the data in lexicographical order
	# ND points for clustering (seperate "main sequence" from quenched galaxies)
	NN = len(X[:,0])
	XC = np.zeros((NN,2))
	XC[:,0] = X[:,0]
	XC[:,1] = X[:,KK]
	# (X,Y) for regression
	Y = X[:,KK]
	X = X[:,:KK]

	### clustering, identify "main sequence"
	clt = DBSCAN(eps=0.1,min_samples=30)
	clt.fit(XC)

	# clustering plot
	#plt.figure()
	#plt.plot(X[:,0],Y,'co',alpha=0.3)
	for i in np.unique(clt.labels_):
		print("cluster label = {} {}".format(i,len(Y[clt.labels_==i])*1.0/NN))
	#	plt.plot(X[clt.labels_==i,0],Y[clt.labels_==i],'o',alpha=0.3,label=str(i))
	#plt.legend()
	#plt.xlabel(r'$\log (Z/Z_\odot)$')
	#plt.ylabel(r'$\log {\rm DGR}$')
	#plt.savefig('clustering.png',dpi=300)

	# only fit "main sequence"
	X = X[clt.labels_==0]
	Y = Y[clt.labels_==0]

	### load parameter list and filter out unused galaxy properties
	X,labels = filt_gal_prop(X,ipar,labels)

	return X,Y,labels

def split_train_cv(X,Y):
# split training and c.v. sets
	NN = len(Y)
	ID = np.unique((np.random.random(int(NN*7/10))*(NN-1)).astype(int))
	ID_cv = np.arange(NN)
	ID_cv = np.delete(ID_cv,ID)

	return X[ID],Y[ID],X[ID_cv],Y[ID_cv]


def mse(X,X_pred):
	return np.mean((X-X_pred)**2)

def regression(X,Y,X_cv,Y_cv,n_estimators=200,criterion='mse',max_depth=20):
####### ERT regressor trained by simulated data
	KK = X.shape[1]
	mse_train = []
	mse_cv = []
	importance = np.zeros(KK)

	regr = ETR(n_estimators=n_estimators,criterion=criterion,max_depth=max_depth)
	regr.fit(X,Y)
	Y_pred = regr.predict(X)
	mse_train.append(np.sqrt(mse(Y,Y_pred)))
	Y_pred = regr.predict(X_cv)
	mse_cv.append(np.sqrt(mse(Y_cv,Y_pred)))
	importance[:] = regr.feature_importances_

	return  regr, mse_train, mse_cv

def predict(regr,x):
#### predict DGR using trained ERT rgr and input parameters x
# question: how to propagate the uncertainty?
	return regr.predict(x)

def dgr_ert(metallicity,sfr,gas_mass):

        #groom the input arrays 
        #stellar_mass = stellar_mass.in_units('Msun').value
        gas_mass = gas_mass.in_units('Msun').value
        metallicity = metallicity.value
        sfr = sfr.in_units('Msun/yr').value

        #set values to lowest nonzero value so that the ERT doesn't choke in empty cells
        metallicity[np.where(metallicity == 0)[0]] = np.min(metallicity[np.nonzero(metallicity)])
        #stellar_mass[np.where(stellar_mass == 0)[0]] = np.min(stellar_mass[np.nonzero(stellar_mass)])
        sfr[np.where(sfr == 0)[0]] = np.min(sfr[np.nonzero(sfr)])
        gas_mass[np.where(gas_mass == 0)[0]] = np.min(gas_mass[np.nonzero(gas_mass)])

        param = cfg.par.pd_source_dir+'powderday/mlt/mlt_param.list'
        ipar = np.loadtxt(param)
        #tab  = np.array([metallicity,stellar_mass,gas_mass])
        tab = np.array([metallicity,sfr,gas_mass])
        #tab = np.array([metallicity,sfr])
        tab = tab.T


        X_in = input_prepare(ipar,tab)
        X,Y,labels = matrix_prepare(ipar,sim=cfg.par.pd_source_dir+'powderday/mlt/m12n256.hdf5')
        X,Y,X_cv,Y_cv = split_train_cv(X,Y)
        regr,mse_train,mse_cv = regression(X,Y,X_cv,Y_cv,n_estimators=400,criterion='mse',max_depth=20)
        print("mse_train, mse_c.v. = {} {}".format(mse_train,mse_cv))
	
        Y_out = predict(regr,X_in)

        '''
        if len(sys.argv) > 3:
                out_dgr = sys.argv[3]
                np.savetxt(out_dgr,Y_out)
        else:
                print(Y_out)
        '''
        return Y_out
