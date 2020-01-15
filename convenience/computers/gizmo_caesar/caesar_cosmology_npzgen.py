#caesar_cosmology_npzgen.py

#purpose: to read in a series of caesar files for a cosmology run and
#write an npz file that has dictionaries for the top N galaxies and
#their positions
from __future__ import print_function
import numpy as np
import yt
import caesar
import ipdb
from glob2 import glob

#directory = '/ufrc/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/Groups//'
directory = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512_filtered/output/Groups/'
NGALAXIES_MAX = 10000
TESTING = False
outfile = '/ufrc/narayanan/desika.narayanan/pd_runs/simba/m25n512/simba_m25n512.galaxies_pos_for_pd.npz'
DRIVER_FORMAT = False


MEMBERS = np.sort(glob('%s/caesar*.hdf5' % (directory)))

pos = {}
ngalaxies = {}
for nh in range(NGALAXIES_MAX):
    pos['galaxy'+str(nh)] = {}


if TESTING:MEMBERS=[MEMBERS[-1]]


for file in MEMBERS:

    #this gets snaps in 3 digit format
    if DRIVER_FORMAT: #if the caesar files were written with driver
        snapnum = file[file.find('caesar_')+8:file.find('_z')] 
    else: #else they were written with the CLI
        snapnum = file[file.find('caesar_')+16:file.find('.hdf5')]
    #this try/except is in case the caesar was run on snaps so early there are no galaxies
    try:
        obj = caesar.load(file)
    
    
        if obj.ngalaxies < NGALAXIES_MAX:
            NGALAXIES = obj.ngalaxies
        else:
            NGALAXIES = NGALAXIES_MAX

        #this gets snaps in 3 digit format
        
        ngalaxies['snap'+snapnum] = NGALAXIES
    
   
        for nh in range(NGALAXIES):
            pos['galaxy'+str(nh)]['snap'+snapnum] = obj.galaxies[nh].pos.in_units('code_length').value

    except AttributeError as e:
        pos['galaxy0']['snap'+snapnum] = np.asarray([-1,-1,-1])
    
np.savez(outfile,ngalaxies=ngalaxies,pos=pos)
