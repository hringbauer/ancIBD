import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import socket as socket
import os as os
import sys as sys
import multiprocessing as mp
import h5py
import allel

#sys.path.append("/n/groups/reich/hringbauer/git/hapBLOCK/python3/")  # hack to get development package first in path
from ancIBD.IO.h5_load import get_coverage,get_genos_pairs,get_idx_iid,get_idx_iid_exact, get_markers_good, get_opp_homos_f

def opp_homos(g1, g2):
    """Return opposing homozygotes"""
    o1 = (g1 == 0) & (g2 == 2)
    o2 = (g1 == 2) & (g2 == 0)
    return (o1 | o2)

def plot_hets(map_het=[], het=[], ax=None, het_c="slateblue", c_roh="seagreen", 
              figsize=(14,2), cm_lim=[], ylim=[-0.1, 1.1], fs = 12,
              alpha=0.3, ms=1, lw = 12, title="", min_cm=0, xtickl=True,
              xlabel="centimorgan", ylabel = f"Opp. Homozygote (y/n)", 
              plot=True, savepath=""):
    """Plot Heterozygotes against genenetic map,
    plus ROH calls.
    lw: Linewidth of ROH""" 
    if ax==None:
        plt.figure(figsize=figsize)
        ax = plt.gca()
    ax.plot(map_het*100, (het * 1.1 - 0.05), "o", ms=ms, alpha=alpha, zorder=0, color=het_c)
    ax.set_xlabel(xlabel, fontsize=fs)
    if not xtickl:
        ax.set_xticklabels([])
    ax.set_ylim(ylim)
    
    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(np.array([1,0]) * 1.1 - 0.05)
    ax2.set_yticklabels([])
    ax.set_yticklabels([])
    ax2.set_ylabel(ylabel, fontsize=fs*0.8, color=het_c, rotation=270, labelpad=fs)
    
    if len(cm_lim)==2:
        ax.set_xlim(cm_lim)
        
    if len(title)>0:
        ax.set_title(title, fontsize=fs)
        
    if len(savepath)>0:
        plt.savefig(savepath, bbox_inches = 'tight', pad_inches = 0, dpi=300)
        print(f"Saved figure to: {savepath}")
        #plt.savefig(folder + "posterior_cm.png", bbox_inches = 'tight', pad_inches = 0, dpi=300)
    
    if plot:
        plt.show()  
        
def plot_opp_homos(path_h5 = f"/n/groups/reich/hringbauer/hapsburg_runs/data/data_eirini/h5/all_ch", ch = 16,
                  iids = ["MYG001.A0101", "MYG006.A0101"], cutoff=0.99, output=True, exact=True,
                  cm_lim=[], ms=4, savepath=""):
    """Plot pairwise opp. homoyzgotes. Loads f5 data, and calls and plots opp. homos"""
    path_h5 = f"{path_h5}{ch}.h5"
    f = h5py.File(path_h5, "r") # Load for Sanity Check. See below!
    g1, g2, m = get_genos_pairs(f, sample1=iids[0], sample2=iids[1], 
                                cutoff=cutoff, output=output, exact=exact)
    o_homos = opp_homos(g1, g2)
    plot_hets(m, o_homos, title=f"{iids[0]} - {iids[1]}, Chromosome {ch}, Opposing Homozygotes", ms=ms,
              cm_lim=cm_lim, savepath=savepath) # ./figs/principle_proof/siblings.png
    
def plot_identical_gts(path_h5 = f"/n/groups/reich/hringbauer/hapsburg_runs/data/data_eirini/h5/all_ch", ch = 16,
                  iids = ["MYG001.A0101", "MYG006.A0101"], cutoff=0.99, output=True, exact=True,
                  cm_lim=[], ms=4, savepath=""):
    """Plot pairwise opp. homoyzgotes. Loads f5 data, and calls and plots opp. homos"""
    path_h5 = f"{path_h5}{ch}.h5"
    f = h5py.File(path_h5, "r") # Load for Sanity Check. See below!
    g1, g2, m = get_genos_pairs(f, sample1=iids[0], sample2=iids[1], 
                                cutoff=cutoff, output=output, exact=exact)
    diff = (g1!=g2)
    plot_hets(m, diff, title=f"{iids[0]} - {iids[1]}, Chromosome {ch}, Different Diploid Genotyeps", ms=ms,
              cm_lim=cm_lim, ylabel = f"Diff. Genotypes (y/n)", savepath=savepath) # ./figs/principle_proof/siblings.png