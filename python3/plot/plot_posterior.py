import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np

from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'   # Set the default
rcParams['font.sans-serif'] = ['Arial']  # Make sure to have the font installed (it is on cluster for Harald)

def plot_posterior(ax=0, post=[], state=0, figsize=(12,3), c="maroon", 
                   xlim=[], ylim=[], ylabel="Posterior", xlabel="Position",
                   lw=3, fs_l=12, show=True, title="", savepath=""):
    """Plot Posterior [k,l] array"""
    if ax==0:
        plt.figure(figsize=figsize)
        ax=plt.gca()
    ax.plot(np.exp(post[state,:]), color=c, lw=lw)
    
    ### Do optional plotting
    if len(xlabel)>0:
        ax.set_xlabel(xlabel, fontsize=fs_l)
    if len(ylabel)>0:
        ax.set_ylabel(ylabel, fontsize=fs_l)
    if len(xlim)>0:
        ax.set_xlim(xlim)
    if len(ylim)>0:
        ax.set_ylim(ylim)
    if len(title)>0:
        ax.set_title(title, fontsize=fs_l)
    if len(savepath)>0:
        plt.savefig(savepath, bbox_inches ='tight', pad_inches = 0, dpi=400)
        print(f"Saved to {savepath}")
    if show:
        plt.show()
    else: 
        return ax
    
def plot_posterior_panel(post=[], figsize=(12,9), c="maroon", c_hw="g",
                         hspace=0.12, wspace=0.15, xlim=[], ylim=[-0.05,1.05], lw=3, fs_l=12,
                         savepath=""):
    """Plot Posterior [k,l] array"""
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(ncols=2, nrows=3, figure=fig, 
                           wspace=wspace, hspace=hspace)
    ax = fig.add_subplot(gs[0, 0])

    pos = [[1,0], [1,1], [2,0], [2,1]]
    axs= [fig.add_subplot(gs[i[0], i[1]]) for i in pos]
    labels = ["0-0", "1-0", "1-0", "1-1"]
    labels = ["Posterior:" + l + " IBD" for l in labels]
    
    plot_posterior(ax, post=post, state=0, c=c_hw, 
                       xlim=xlim, ylim=ylim, xlabel="",
                       ylabel="Posterior: Non-IBD", lw=lw, show=False)
    
    for i, ax in enumerate(axs):
        plot_posterior(ax, post=post, state=i+1, c=c, 
                       xlim=xlim, ylim=ylim, xlabel="",
                       ylabel=labels[i], lw=lw, show=False)
    if len(savepath)>0:
        plt.savefig(savepath, bbox_inches ='tight', pad_inches = 0, dpi=400)
        print(f"Saved to {savepath}")
    plt.show()