import matplotlib.pyplot as plt
import numpy as np
import h5py as h5py
import pandas as pd

from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'   # Set the default
rcParams['font.sans-serif'] = ['Arial']  # Make sure to have the font installed (it is on cluster for Harald)



def plot_scatter_ibd(df_ibds, bins=[], min_cm=12, 
                     xlim=[-5,3600], ylim=[0,50], 
                     title="", savepath=""):
    """Plot Scatter Plot of IBD"""

    if len(bins)==0:
        bins = np.linspace(0,3600,73)

    plt.figure(figsize=(8,8))
    ax = plt.gca()
    ax.scatter(df_ibds[f"sum_IBD>{min_cm}"], df_ibds[f"n_IBD>{min_cm}"], s=40,
               ec="k", linewidth=0.5, color="deepskyblue")

    ax.set_xlabel(f"Sum IBD >{min_cm}cM [cM]", fontsize=14)
    ax.set_ylabel(f"n IBD >{min_cm}cM [cM]", fontsize=14)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    if len(savepath)>0:
        plt.savefig(savepath, bbox_inches ='tight', pad_inches = 0, dpi=400)
        print(f"Saved to {savepath}")
    if len(title)>0:
        ax.set_title(title)
    plt.show()