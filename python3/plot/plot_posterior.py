import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import h5py as h5py

from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'   # Set the default
rcParams['font.sans-serif'] = ['Arial']  # Make sure to have the font installed (it is on cluster for Harald)

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