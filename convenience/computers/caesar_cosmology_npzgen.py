#caesar_cosmology_npzgen.py

#purpose: to read in a series of caesar files for a cosmology run and
#write an npz file that has dictionaries for the top N galaxies and
#their positions

import numpy as np
import yt
import caesar
import ipdb
from glob2 import glob

directory = '/Volumes/pegasus2/gizmo_runs/mufasa/m25n512/fh_qr/'
NHALOS = 1000
TESTING = False
outfile = '/Volumes/pegasus/pd_runs/mufasa/m25n512/fh_qr/attenuation/mufasa_m25n512.halos_pos_for_pd.npz'



MEMBERS = np.sort(glob('%s/output/Groups/caesar*.hdf5' % (directory)))

pos = {}

for nh in range(NHALOS):
    pos['halo'+str(nh)] = {}

if TESTING:MEMBERS=[MEMBERS[-1]]

for file in MEMBERS:
    obj = caesar.load(file)
    snapnum = file[file.find('.hdf5')-3:file.find('.hdf5')]

    for nh in range(NHALOS):
        pos['halo'+str(nh)]['snap'+snapnum] = obj.halos[nh].pos.in_units('code_length').value
    
    
    
np.savez(outfile,NHALOS=NHALOS,pos=pos)
