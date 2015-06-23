import numpy as np
from yt.mods import *
import matplotlib.pyplot as plt
import config as cfg
import pdb

def proj_plots(pf):

    
    print '\n[plot_generate/proj_plots] Saving Diagnostic Projection Plots \n'
    p = ProjectionPlot(pf,"z",("deposit","PartType0_smoothed_density"))
    p.save(cfg.par.PD_output_dir+'/proj_plot_x.png')
    p = ProjectionPlot(pf,"x",("deposit","PartType0_smoothed_density"))
    p.save(cfg.par.PD_output_dir+'/proj_plot_y.png')
    p = ProjectionPlot(pf,"y",("deposit","PartType0_smoothed_density"))
    p.save(cfg.par.PD_output_dir+'/proj_plot_z.png')
 
    return None



def mass_weighted_distribution(array,masses,fileout='junk.png',nbins=100):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.hist(array,bins=nbins,weights=masses,log=True,normed=True)
    
    fig.savefig(fileout)

