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

directory = '/ufrc/narayanan/desika.narayanan/gizmo_runs/mufasa/simba/m50n512/Groups//'
NGALAXIES_MAX = 10000
TESTING = False
outfile = '/ufrc/narayanan/desika.narayanan/pd_runs/simba/m50n512/simba_m50n512.galaxies_pos_for_pd.npz'



MEMBERS = np.sort(glob('%s/caesar*.hdf5' % (directory)))

pos = {}
ngalaxies = {}
for nh in range(NGALAXIES_MAX):
    pos['galaxy'+str(nh)] = {}


if TESTING:MEMBERS=[MEMBERS[-1]]


for file in MEMBERS:
    obj = caesar.load(file)

    if obj.ngalaxies < NGALAXIES_MAX:
        NGALAXIES = obj.ngalaxies
    else:
        NGALAXIES = NGALAXIES_MAX


    snapnum = file[file.find('.hdf5')-3:file.find('.hdf5')]

    ngalaxies['snap'+snapnum] = NGALAXIES
    
   
    for nh in range(NGALAXIES):
        pos['galaxy'+str(nh)]['snap'+snapnum] = obj.galaxies[nh].pos.in_units('code_length').value

    
    
np.savez(outfile,ngalaxies=ngalaxies,pos=pos)
