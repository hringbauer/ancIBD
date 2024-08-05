import matplotlib.pyplot as plt
import numpy as np
import h5py as h5py
import pandas as pd

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
    
    
def plot_sactter_ibd_target(df_ibds, df_target=[], bins=[], min_cm=12, s=40,
                            xlim=[-5,3600], ylim=[0,50], c_target = "maroon", c_back="gray",
                            title="", savepath=""):
    """Plot Scatter plot of sum and number of IBD segments (x & y axis).
    df_ibds: Background plot
    df_target: Foreground plot"""
    
    if len(bins)==0:
        bins = np.linspace(0,3600,73)

    plt.figure(figsize=(8,8))
    ax = plt.gca()
    ax.scatter(df_ibds[f"sum_IBD>{min_cm}"], df_ibds[f"n_IBD>{min_cm}"], s=s,
               ec="k", linewidth=0.3, color=c_back, zorder=0)
    
    ### Plot Foreground
    ax.scatter(df_target[f"sum_IBD>{min_cm}"], df_target[f"n_IBD>{min_cm}"], s=s*2,
               ec="k", linewidth=0.8, color=c_target, zorder=10)

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
    