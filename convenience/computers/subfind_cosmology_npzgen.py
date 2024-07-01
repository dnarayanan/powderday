#subfind_cosmology_npzgen.py

#purpose: to read in a series of filtered subfind files from arepo, and 
#write an npz file that has dictionaries for the top N galaxies and
#their positions
import numpy as np
import yt
import h5py
from glob2 import glob
import pdb

filter_directory = '/blue/narayanan/desika.narayanan/arepo_runs/m25n512b_09c6391_nolimits_filter_ready/output/filtered_snaps/'
outfile = '/blue/narayanan/desika.narayanan/arepo_runs/m25n512b_09c6391_nolimits_filter_ready/m25n512.galaxies_pos.npz'
NGALAXIES_MAX = 10000


#first find the snapnums based on the subdirs
snapdirs = np.sort(glob('%s/snap*' %(filter_directory)))
snaplist = []
for snapdir in snapdirs:
    snaplist.append(snapdir[snapdir.find('/snap')+5:snapdir.find('/snap')+8])


#first find the snapnums based on the 0th galaxy having been run
#zeroth_galaxy_files =  np.sort(glob('%s/*.0.hdf5' % (filter_directory)))
#snaplist = []
#for file in zeroth_galaxy_files:
#    snapnum = file[file.find('snapshot_')+9:file.find('.filtered')]
#    snaplist.append(snapnum)

ngalaxies = {}
pos = {}
for nh in range(NGALAXIES_MAX):
    pos['galaxy'+str(nh)] = {}

for snapnum in snaplist:
    galaxy_filtered_files = np.sort(glob('%s/*.hdf5' % (filter_directory+'/snap'+str(snapnum))))
    num_galaxies = len(galaxy_filtered_files)
    
    for galaxy in galaxy_filtered_files:

        print(galaxy)
        galnum = int(galaxy[galaxy.find('galaxy_')+7:galaxy.find('.hdf5')])
        
        #set the number of galaxies
        if num_galaxies < NGALAXIES_MAX:
            NGALAXIES = num_galaxies
        else:
            NGALAXIES = NGALAXIES_MAX

        ngalaxies['snap'+snapnum] = NGALAXIES

        
        #now open the galaxy file and find its COM in code units
        file = h5py.File(galaxy,"r")
        particle_pos = file['PartType0']['Coordinates'][:]
        particle_mass = file['PartType0']['Masses'][:]
        xpos,ypos,zpos = particle_pos[:,0],particle_pos[:,1],particle_pos[:,2]

        try:
            xpos_mean = np.average(xpos,weights=particle_mass)
            ypos_mean = np.average(ypos,weights=particle_mass)
            zpos_mean = np.average(zpos,weights=particle_mass)
        except:
            print("no galaxy positions found")
            xpos_mean,ypos_mean,zpos_mean = -1,-1,-1
        
        pos['galaxy'+str(galnum)]['snap'+snapnum] = np.asarray([xpos_mean,ypos_mean,zpos_mean])

np.savez(outfile,ngalaxies=ngalaxies,pos=pos)
