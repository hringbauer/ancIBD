import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import h5py as h5py
from matplotlib.patches import Rectangle


from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'   # Set the default
rcParams['font.sans-serif'] = ['Arial']  # Make sure to have the font installed (it is on cluster for Harald)

import sys
sys.path.append("/n/groups/reich/hringbauer/git/hapBLOCK/python3/")  # hack to get development package first in path
from IO.h5_load import get_coverage,get_genos_pairs,get_idx_iid,get_idx_iid_exact, get_markers_good, get_opp_homos_f

def plot_posterior(ax=0, morgan=[], post=[], het=[], het_m=[], 
                   df_ibd=[], df_truth=[], state=0, figsize=(12,3), 
                   xlim=[], ylim=[-0.08,1.27], ylabel="Posterior", xlabel="Position",
                   c="maroon", het_c="gray", c_truth="green", ms=1,
                   lw=3, lw_ibd=10, c_ibd="slateblue", y_ibd=1.2, dpi=400, 
                   fs_l=12, show=True, min_cm=4, title="", savepath=""):
    """Plot Posterior [k,l] array. If morgan given, plot in centimorgan.
    Can then also plot hapROH formatted IBD blocks (df_ibd).
    And plot ground truth hapROH formatted IBD blocks (df_truth).
    If het is given [array boolean], plot het using het_m coordinates"""
    if ax==0:
        plt.figure(figsize=figsize)
        ax=plt.gca()
    if len(morgan)==0:
        morgan = np.arange(np.shape(post)[1])
    ax.plot(morgan*100, post[state,:], color=c, lw=lw)
    ax.set_yticks([0., 0.2, 0.4, 0.6, 0.8, 1.0])
    ### Do optional plotting
    # Hets
    if len(xlabel)>0:
        ax.set_xlabel(xlabel, fontsize=fs_l)
    if len(ylabel)>0:
        ax.set_ylabel(ylabel, fontsize=fs_l)
    if len(xlim)>0:
        ax.set_xlim(xlim)
    if len(ylim)>0:
        ax.set_ylim(ylim)
    if len(het)>0:
        plot_hets(ax, het_m, het, ms=ms, het_c=het_c, fs_l=fs_l, ylim=ylim)
    if len(title)>0:
        ax.set_title(title, fontsize=fs_l)
    if len(df_ibd)>0:
        df_ibd = df_ibd[df_ibd["lengthM"] > min_cm/100]  # Filter out long enough calls
        ax.hlines(y=[y_ibd]*len(df_ibd), xmin=100 * df_ibd["StartM"], xmax= 100 * df_ibd["EndM"], 
                        colors=c_ibd, linewidth=lw_ibd)
    if len(df_truth)>0:  # Plot Ground truth dataframe
        ### Plot them
        plt.hlines(y=[1.12]*len(df_truth), xmin=100 * df_truth["IBD_Begin"], 
                   xmax=100 * df_truth["IBD_End"], colors=c_truth, linewidth=lw_ibd)
        
    if len(savepath)>0:
        plt.savefig(savepath, bbox_inches ='tight', pad_inches = 0, dpi=dpi)
        print(f"Saved to {savepath}")
    if show:
        plt.show()
    else: 
        return ax
    
def plot_hets(ax, het_m, het, alpha=0.3, ms=1, het_c="slateblue",
              ylabel = "Opp. Homozygotes (no/yes)", fs_l=12, ylim=[]):
    """Plot Heterozygote Markers onto Axis"""
    ax2 = ax.twinx()
    ax2.plot(het_m*100, (het * 1.1 - 0.05), "o", ms=ms, alpha=alpha, zorder=0, color=het_c)
    ax2.set_yticks([-0.05, 1.05])
    ax2.set_yticklabels([])
    ax2.set_ylabel(ylabel, fontsize=fs_l, color=het_c)
    if len(ylim)>0:
        ax2.set_ylim(ylim)
    
########################################################
    
def plot_posterior_panel(post=[], figsize=(12,9), c="maroon", c_hw="g",
                         hspace=0.12, wspace=0.15, xlim=[], ylim=[-0.05,1.05], 
                         lw=3, fs_l=12, ch=0, savepath=""):
    """Plot Posterior [k,l] array.
    ch: If bigger 0, load 1240K map postion"""
    
    if ch>0:
        m = get_map(ch=ch)
        assert(len(m)==np.shape(post)[1]) # Sanity Check
    else:
        m=[]
        
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(ncols=2, nrows=3, figure=fig, 
                           wspace=wspace, hspace=hspace)
    ax = fig.add_subplot(gs[0, 0])

    pos = [[1,0], [1,1], [2,0], [2,1]]
    axs= [fig.add_subplot(gs[i[0], i[1]]) for i in pos]
    labels = ["0-0", "1-0", "0-1", "1-1"]
    labels = ["Posterior:" + l + " IBD" for l in labels]
    
    plot_posterior(ax, morgan=m, post=post, state=0, c=c_hw, 
                       xlim=xlim, ylim=ylim, xlabel="",
                       ylabel="Posterior: Non-IBD", lw=lw, show=False)
    
    for i, ax in enumerate(axs):
        plot_posterior(ax, morgan=m, post=post, state=i+1, c=c, 
                       xlim=xlim, ylim=ylim, xlabel="",
                       ylabel=labels[i], lw=lw, show=False)
    if len(savepath)>0:
        plt.savefig(savepath, bbox_inches ='tight', pad_inches = 0, dpi=400)
        print(f"Saved to {savepath}")
    plt.show()
    
def get_map(path_h5= "./data/hdf5/1240k_v43/ch", ch=3, 
            col_map="variants/MAP"):
    """Get Map position in Morgan"""
    path_load = f"{path_h5}{ch}.h5"
    f = h5py.File(path_load, "r") # Load for Sanity Check. See below!
    m = f[col_map][:]
    return m

def plot_posterior_5States(basepath, start=-1, end=-1):
    # I assume the two files, map.npy and posterior.npy, reside in basepath
    r_map = np.load(f'{basepath}/map.npy', 'r')
    r_map = 100*r_map # cM is more intuitive for me
    post = np.load(f'{basepath}/posterior.npy', 'r')
    _, l = post.shape
    assert(len(r_map) == l) # sanity check

    i = np.searchsorted(r_map, start) if start != -1 else 0
    j = np.searchsorted(r_map, end) if end != -1 else -1

    for k in range(1,5):
        plt.plot(r_map[i:j], post[k, i:j], label=f'state {k}', linewidth=0.25, alpha=0.75)
    plt.plot(r_map[i:j], 1-post[0, i:j], label='sum of IBD states', color='black', linewidth=0.6)
    plt.xlabel('Genomic Position')
    plt.ylabel('Posterior')
    plt.legend(loc='upper right', fontsize='xx-small')
    plt.savefig(f'{basepath}/posterior.png', dpi=300)
    plt.clf()
    
def plot_posterior_7States(basepath, start=-1, end=-1, prefix="", simplify=True):
    # I assume the two files, map.npy and posterior.npy, reside in basepath
    r_map = np.load(f'{basepath}/map.npy', 'r')
    r_map = 100*r_map # cM is more intuitive for me
    post = np.load(f'{basepath}/posterior.npy', 'r')
    _, l = post.shape
    assert(len(r_map) == l) # sanity check

    i = np.searchsorted(r_map, start) if start != -1 else 0
    j = np.searchsorted(r_map, end) if end != -1 else -1

    if not simplify:
        for k in range(1,5):
            plt.plot(r_map[i:j], post[k, i:j], label=f'state {k}', linewidth=0.25, alpha=0.75)
    plt.plot(r_map[i:j], np.sum(post[1:5, i:j], axis=0), label='sum of IBD1 states', color='black', linewidth=0.75)
    plt.plot(r_map[i:j], np.sum(post[5:7, i:j], axis=0), label='sum of IBD2 states', color='grey', linewidth=0.75)


    plt.xlabel('Genomic Position')
    plt.ylabel('Posterior')
    plt.legend(loc='upper right', fontsize='xx-small')

    if len(prefix) == 0:
        plt.savefig(f'{basepath}/posterior.png', dpi=300)
    else:
        plt.savefig(f'{basepath}/posterior.{prefix}.png', dpi=300)
    plt.clf()
    
def plot_posterior_7States_plusGeno(basepath, start=-1, end=-1, iids=[], truth=None, path2hdf5="", ch=1, prefix="", outFolder=""):
    # I assume the two files, map.npy and posterior.npy, reside in basepath
    r_map = np.load(f'{basepath}/map.npy', 'r')
    r_map = 100*r_map # cM is more intuitive for me
    post = np.load(f'{basepath}/posterior.npy', 'r')
    _, l = post.shape
    assert(len(r_map) == l) # sanity check

    i = np.searchsorted(r_map, start) if start != -1 else 0
    j = np.searchsorted(r_map, end) if end != -1 else -1
    plt.plot(r_map[i:j], np.sum(post[1:5, i:j], axis=0), label='sum of IBD1 states', color='black', linewidth=0.75)
    plt.plot(r_map[i:j], np.sum(post[5:7, i:j], axis=0), label='sum of IBD2 states', color='grey', linewidth=0.75)

    ########################### grabbing genotypes #########################
    f = h5py.File(f'{path2hdf5}{ch}.h5', 'r')
    gt1, gt2, map = get_genos_pairs(f, sample1=iids[0], sample2=iids[1], cutoff=0.99, output=True, phased=False, exact=False)
    map = 100*map
    i = np.searchsorted(map, start) if start != -1 else 0
    j = np.searchsorted(map, end) if end != -1 else -1
    plt.scatter(map[i:j], [-0.05]*(j-i), marker='o', color='indigo', s=0.5, alpha=0.3) # to get a sense of the density of markers
    gt1, gt2, map = gt1[i:j], gt2[i:j], map[i:j]
    diff_gts = (gt1 != gt2)
    plt.scatter(map[diff_gts], [1.05]*np.sum(diff_gts), marker='o', color='blue', label='Diff. Genotypes', s=0.5, alpha=0.3)

    # plot oppo homozygotes
    index = np.logical_or(np.logical_and(gt1 == 0, gt2 == 2), np.logical_and(gt1 == 2, gt2 == 0))
    plt.scatter(map[index], [1.1]*np.sum(index), marker='o', color='red', label='oppo homo', s=0.5, alpha=0.3) 

    # plot truth IBD region if provided
    if truth:
        print('plot IBD new version')
        plt.hlines(y=1.15, xmin=truth[0], xmax=truth[1], 
                        colors="blue", linewidth=6)



    plt.xlabel('Genomic Position')
    plt.ylabel('Posterior')
    plt.legend(loc='upper right', fontsize='xx-small')
    if len(outFolder) > 0:
        basepath = outFolder
    if len(prefix) == 0:
        plt.savefig(f'{basepath}/posterior.png', dpi=300)
        plt.savefig(f'{basepath}/posterior.pdf', dpi=300)
    else:
        plt.savefig(f'{basepath}/posterior.{prefix}.png', dpi=300)
        plt.savefig(f'{basepath}/posterior.{prefix}.pdf', dpi=300)
    plt.clf()