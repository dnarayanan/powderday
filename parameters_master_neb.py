# coding=utf-8

#-----------------------------------------------------------
# Parameters that can tuned for running CLOUDY models individually for each star particle. 
# There is an option to include nebular emission for Young stars, post-AGB stars, Diffuse Ionized Gas and AGNs.
# These only work when add_nebular_emission is set to true and use_cloudy_tables is set to false in parameters_master.py
#-----------------------------------------------------------

#===============================================
# Sub Resolution Modeling 
#===============================================
use_cmdf = False                            # If True, star particles that have mass greater than cmdf_mas_mass (defined below) are broken down using a 
                                            # cluster mass distribution function (cmdf) of the form dN/dM goes as M^(beta). This works irrespecitve of whether
                                            # nebular emission is turned on or not.  The cmdf is set by the following parameters defined below: 
                                            # cmdf_min_mass, cmdf_max_mass, cmdf_bins and cmdf_beta.

cmdf_min_mass = 3.5                         # Minimum mass of the star clusters in units of log(Msun). Note: Results might be inconsistent if
                                            # set lower than 3.5. (See Chandar et al.2014 for more info) (Default = 3.5)

cmdf_max_mass = 5.0                         # Maximum mass of the star clusters in units of log(Msun). (Default = 5.0). Note: Only star particles that
                                            # have a mass greater than this parameter are broken down. 

cmdf_bins = 6                               # The number of bins used for calulating the cluster mass distribution function (Default = 6.0)

cmdf_beta = -2.0                            # Beta (power law exponent) for calculating CMDF (dN/dM goes as M^(beta)) 

use_age_distribution = False                # Setting this to True, divides the star particles with ages between age_dist_min and age_dist_max (next parameters) into 
                                            # an ensemble of particles all of whom have the same properties except their age which is picked from a power law age 
                                            # distribution of the form dN/dt is proportional to t^-0.65 (Imp: This can only be used if use_cmdf is also set to True). 
                                            # Note: The function has a bunch of tunable parameters that can be changed though we feel that their default values
                                            # should be good enough for most cases. The function is located in cloudy_tools.py file under powderday/nebular_emission. 

age_dist_min = 3e-3                         # Star particle above this age are sub-divided into an age distribution if use_age_distribution is set to True
                                            # (Units: Gyr, Default = 3.e-3)

age_dist_max = 1e-2                         # Star particles below this age are sub-divided into an age distribution if use_age_distribution is set to True
                                            # (Units: Gyr, Default = 1.e-2)
                                            
#***********************
# COMMON PARAMETERS
#***********************
# NOTE: These parmeters take either three or four values as an input. 
# They correspond to the value of the pararmeter for young_stars, PAGB stars, AGN and DIG respectively.

FORCE_gas_logu = [False, False, False] 	    # If set, then we force the ionization parameter (gas_logu) to be 
                            			    # gas_logu (next parameter) else, it is taken to be variable and dependent on ionizing 
                            			    # radiation from star particles. (Default: [False, False, False])

gas_logu = [-2.0, -2.0, -2.0]        		# Gas ionization parameter. This is only relevant 
                            			    # if add_neb_emission is set to True and FORCE_gas_logu is set to True (Default: [-2.0, -2.0, -2.0])

gas_logu_init = [0.0, 0.0, 0.0]        	    # Force the ionization parameter to increase/decrease by this value (Scale: log). 
                            		    	# Useful if you want to run tests (Default: [0.0, 0.0, 0.0])

FORCE_gas_logz = [False, False, False]      # If set, then we force the metallicity (gas_logz) to be gas_logz (next parameter)
                            	            # else, it is taken to be the star particles metallicity. (Default: [False, False, False])

gas_logz = [0.0, 0.0, 0.0]  			    # Metallicity of the HII region in units of log(Z/Z_sun)
                            			    # only relevant if add_neb_emission = True and FORCE_gas_logz = True (Default: [0.0, 0.0, 0.0])

FORCE_logq = [False, False, False]      	# If set, then we force the number of ionizing photons to be source_logq (next parameter)
                                            # else, it is taken to be variable and dependent on ionizing radiation of the source. (Default: [False, False, False])

source_logq = [1.e47, 1.e47,1.e47]          # The number of ionizing photons emitted by the source in units of s^-1. Only relevant if add_neb_emission = True, 
    										# use_cloudy_tables = True and  FORCE_gas_logq = True (Default: [1.e47,1.e47,1.e47])  
                                            
FORCE_inner_radius = [False, False, True]   # If set, then we force the inner radius of the cloud to be inner_radius (next parameter). 
											# IMP Note: This works only for young stars and Post-AGB stars. 
    										# For AGN we keep the inner radius fixed at whatever is set by inner_radius (next parameter) 
    										# irrespective of what this parameter is set to. (Default: [False,False,True])

inner_radius = [1.e19, 1.e19, 2.777e+20]   	# This sets the inner radius of the cloud in cm. This is used only when add_neb_emission = True,
                            		    	# use_cloudy_tables = False and FORCE_inner_radius = True (Default: [1.e19, 1.e19, 2.777e+20], Units = cm)

FORCE_N_O_Pilyugin = [False, False, False, False]  # If set to True, Nitrogen abundances are set according to the N/O vs O/H relation from Pilyugin et al. 2012
                                                   # If FORCE_N_O ratio (next parameter) is set to True then this parameter is ignored.(Default: [False,False,False, False])

FORCE_N_O_ratio = [False, False, False, False]     # If set, then we force the Nitrogen abundance such that the log of N/O ratio is N_O_ratio (next parameter). 
                            			           # This can be used as a template fix adundance ratio of other elements (Default:  [False, False, False])

N_O_ratio = [-0.85, -0.85, -0.85, -0.85]           # This sets the log of N/O ratio. This is used only when add_neb_emission = True,
                            			           # use_cloudy_tables = False, FORCE_N/O ratio = True and neb_abund = "direct" (Default: = [-0.85, -0.85, -0.85])

neb_abund = ["dopita", "dopita", "dopita", "dopita"]  # This sets the HII region elemental abundances for generating CLOUDY models. 
                            			              # Available abundaces are.
                            			              #    dopita:    Abundances from Dopita (2001) with old solar abundances = 0.019 and ISM grains.
                            			              #    newdopita: Abundances from Dopita (2013). Solar Abundances from Grevasse 2010 - z= 0.013
                            			              #               includes smooth polynomial for N/O, C/O relationship functional form for He(z),
                            			              #               new depletion and factors in ISM grains.
                            			              #    gutkin:    Abundabces from Gutkin (2016) and PARSEC metallicity (Bressan+2012) based on Grevesse+Sauvel (1998) 
                            			              #               and Caffau+2011 
                            			              #    direct:    Abundances are taken directly from the simulation if possible. Defaults to using "dopita" if there is 
                            			              #               an error. (Note: Works only for AGNs and star particles that are added directly without binning.
                            			              #               Make sure to set FORCE_BINNED to False)
                            			              # This is used only when add_neb_emission = True and use_cloudy_tables = False. (Default: ["dopita", "dopita", "dopita"])
#***************************
# YOUNG STARS (HII regions)
#***************************

add_young_stars = True     			        # If set, the young stars are included when calculating nebular emission (Default: True)


HII_Rinner_per_Rs = 0.01        		    # Rinner for cloudy calculations is set to this value times the Stromgen Radius. 
                            			    # For example, if set to 0.01 Rinner is taken to be 1 % of Stromgren Radius. 
                            		    	# If FORCE_inner_radius (next parameter) is set to True then this is overridden
                            			    # and the value set by the inner_radius is used. This parameter is used 
                            			    # only when add_neb_emission = True and use_cloudy_tables = False (Default: 0.01)
   
HII_nh = 1.e2               			    # Gas hydrogen density for calcualting nebular emission in units if cm^-3. 
                            			    # This is used only when add_neb_emission = True and use_cloudy_tables = False (Default = 1.e2)

HII_min_age = 1.e-3                         # Sets the minimum age limit for calculating nebular emission in units of Gyr. 
                                            # This is used only when add_neb_emission = True and use_cloudy_tables = False (Default = 1.e-3)

HII_max_age = 1.e-2         			    # Sets the maximum age limit for calculating nebular emission in units of Gyr. 
                            			    # This is used only when add_neb_emission = True and use_cloudy_tables = False (Default = 1.e-2)

HII_escape_fraction = 0.0   			    # Fraction of H-ionizaing photons that escape the HII region. 
                            			    # This is used only when add_neb_emission = True and use_cloudy_tables = False (Default = 0.0)

HII_dust = False                            # If set, then dust grains are included in the CLOUDY model. We use grains orion command to add
                                            # dust grains which specifies graphitic and silicate grains with a size distribution and abundance
                                            #appropriate for those along the line of sight to the Trapezium stars in Orion (see CLOUDY documentation
                                            # Hazy 1 for more info). (Default: False)
#****************
# Post-AGB STARS
#****************

add_pagb_stars = False      			    # If set, the Post-AGB stars are included when calculating nebular emission
                                            # # This works only when add_neb_emission = True and use_cloudy_tables = False (Default: False)

PAGB_N_enhancement = 0.4    			    # Enhances the Nitrogen abundance Post-AGB stars by increasing the log(N/O) by this value. 
                            			    # This used only when add_neb_emission = True, use_cloudy_tables = False and add_pagb_stars = True (Default = 0.4)  

PAGB_C_enhancement = 0.4    			    # Enhances the Carbon abundance Post-AGB stars by increasing the log(C/O) by this value.
                            			    # This used only when add_neb_emission = True, use_cloudy_tables = False and add_pagb_stars = True (Default = 0.4)

PAGB_Rinner_per_Rs = 0.01        		    # Rinner for cloudy calculations is set to this value times the Stromgen Radius. 
                            			    # For example, if set to 0.01 Rinner is taken to be 1 % of Stromgren Radius. 
                            			    # If FORCE_inner_radius (next parameter) is set to True then this is overridden
                            			    # and the value set by the inner_radius is used. This parameter is used 
                            			    # only when add_neb_emission = True and use_cloudy_tables = False (Default: 0.01)

PAGB_nh = 1.e2               			    # Gas hydrogen density for calcualting nebular emission in units if cm^-3. 
                            			    # This is used only when add_neb_emission = True and use_cloudy_tables = False (Default = 1.e2)

PAGB_min_age = 0.1          		    	# Sets the minimum age limit for calculating nebular emission from post-AGB stars, in units of Gyr.
                            			    # This is used only when add_neb_emission = True, use_cloudy_tables = False and add_pagb_stars = True (Default = 0.1)

PAGB_max_age = 10           			    # Sets the maximum age limit for calculating nebular emission from post-AGB stars, in units of Gyr.
                            			    # This is used only when add_neb_emission = True, use_cloudy_tables = False and add_pagb_stars = True (Default = 10)

PAGB_escape_fraction = 0.0   			    # Fraction of H-ionizaing photons that escape the HII region. 
                            			    # This is used only when add_neb_emission = True and use_cloudy_tables = False (Default = 0.0)

#**************
# AGN
#**************

add_AGN_neb = False				            # If set, AGNs are included when calculating nebular emission.
                                            # # This works only when add_neb_emission = True and use_cloudy_tables = False (Default: False)

AGN_nh = 1.e3					            # Gas hydrogen density for calcualting nebular emission in units if cm^-3. 
                            			    # This is used only when add_neb_emission = True and use_cloudy_tables = False (Default = 1.e2)

AGN_num_gas = 32							# For CLOUDY calculations we use the distance weighted average metallicity of gas particles around the AGN. 
			            					# The number of gas particles used for doing so is set by this parameter. (Default: 32)

#**********************
# DIffused Ionized Gas (DIG)
#**********************

add_DIG_neb = False                         # If set, Contribution from DIG is included when calculating nebular emission 
                                            # This works only when add_neb_emission = True and use_cloudy_tables = False (Default: False)

DIG_nh = 1.e1                               # Gas hydrogen density for calcualting nebular emission in units of cm^-3. (Default: 10)
                                            

DIG_min_logU = -6.0                         # Only gas cells with ionization parameter greater than this are considered for DIG calculation. 
                                            # This is done so as to speed up the calculation by ignoring the cells that do not have enough energy 
                                            # to produce any substantial emission. (Defualt: -6.0)

use_black_sed = False                       # If set, Black et al.(1987) ISRF is used as the input SED shape for DIG CLOUDY calculations 
                                            # else, the input SED shape is calulated by by taking a distance weighted average of the CLOUDY 
                                            # output spectrum of nearby stars. The normalization of the SED is set by the total energy 
                                            # above the lyman limit dumped in each cell. (Default: False)

stars_max_dist = 1                          # Only stars within this distance are considered for getting the input spectrum shape. (Units: Kpc)
                                            # This is used only when use_black_sed = False (Default = 1)

max_stars_num  = 20                         # This sets the upper limit on the number of stars that are used for calculating the input spectrum shape.
                                            # This is used only when use_black_sed = False (Default = 20)
#*************************
# DEBUGGING AND CLEAN UP
#*************************

dump_emlines = False                        # If True, The emission lines are saved in a file before going through the dust radiative transfer. 
                                            # These are the cloudy computed emission line strengths, and are calculated for all lines
                                            # cloudy calculates (i.e. not just those undergoing radiative transfer). The format for the output 
                                            # is a wavelength array,followed by a (nlam+2) list for each nebular emission bearing particle. 
                                            # The +2 in the (nlam+2) list are the O/H ratio and the id of that particle. With id = 0 , 1, 2 and 3 
                                            # corresponds to young stars, PAGB stars, AGN and DIG respectively. There is a convenience package in 
                                            # /convenience to help read in this file. This can be used as a fast way getting emission lines for the 
                                            # purpose of debugging the code. Naming convention: emlines.galaxy*.txt where * is the galaxy number. 
                                            # This works only when add_neb_emission = True (Default: False) 

cloudy_cleanup = True                       # If set to True, all the CLOUDY files will be deleted after the source addition is complete. 
                                            # Only relevant if add_neb_emission = True and use_cloudy_tables = False (Default: True)

#===============================================
#DEBUGGING -THE PERFORMANCE OF THE CODE USING THESE PARAMETERS IS NOT GUARANTEED
#===============================================
NEB_DEBUG = False # Dumps parameters related to nebular line emission in a file for debugging.
                  # The file includes the ionization parameter, number of ionizing photons, 
                  # metallicity, inner radius, stellar mass and age for each particle.
                  # Naming convention: nebular_properties_galaxy*.txt where * is the galaxy number

SAVE_NEB_SEDS = False # If set, the CLOUDY output SEDs are saved in a file 
