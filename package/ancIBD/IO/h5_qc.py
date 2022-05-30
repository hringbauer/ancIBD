"""
Functions to QC imputed data (in hdf5 format)
@ Author: Harald Ringbauer, 2021
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os as os
import sys as sys
import h5py

def get_gp_df(path_h5 = "", chs = [1,], iids = ["iid1",], 
              cutoffs = [0.5,0.6,0.7,0.8,0.9]):
    """Check Genotype Probabilities for list of [iids] and list of [chs].
    Calculate max gp fractions <[cutoffs]. A value >0.01 for c=0.5 seems
    highly problematic, and should be flagged!
    Return summary dataframe for each chromosome and iid.
    path_h5: Path of the hdf5 up to the chr. number."""

    columns = ["iid", "ch", "n"] + cutoffs
    df = pd.DataFrame(columns=columns)

    for ch in chs:
        ### load the data
        print(f"Loading Chromosome {ch}...")
        with h5py.File(f"{path_h5}{ch}.h5", "r") as f: # Load for Sanity Check. See below!
                gt = f["calldata/GT"][:]
                gp = f["calldata/GP"][:]
                samples = f["samples"][:].astype("str")
                cm = f["variants/MAP"][:]

        ### iterate over iids
        for iid in iids:
            idx = (samples==iid)
            assert(np.sum(idx)==1)
            i = np.where(idx)[0][0]
            y = np.max(gp[:,i,:], axis=1)
            vs = [np.mean(y<=c) for c in cutoffs] # Values to append

            ### Append Values
            to_append = [iid, ch, len(y)] + vs
            a_series = pd.Series(to_append, index = df.columns)
            df = df.append(a_series, ignore_index=True)  
    
    return df

def plot_gp(ch = 1, iid = "Sz2",
            path_h5 = "/mnt/archgen/users/hringbauer/data/lango.2021.may/hdf5/ch",
            figsize = (9,2.5),
            color = "maroon", 
            savefolder=""):
    """Plot Max. Genotype Posteriors for
    [iid] on [ch]."""
    ### Load the data
    with h5py.File(f"{path_h5}{ch}.h5", "r") as f: # Load for Sanity Check. See below!
        print(list(f["variants"]))
        print(np.shape(f["calldata/GT"]))
        gt = f["calldata/GT"][:]
        gp = f["calldata/GP"][:]
        samples = f["samples"][:].astype("str")
        cm = f["variants/MAP"][:]
        #iids_found = f["samples"][:].astype("str")
    #print(f"Loaded {len(samples)} Individuals.")
    idx = (samples==iid)
    assert(np.sum(idx)==1)
    i = np.where(idx)[0][0]
    y = np.max(gp[:,i,:], axis=1)

    ### Do the plotting
    plt.figure(figsize=figsize)
    ax = plt.gca()
    ax.plot(cm, y, "o", c=color, ms=1, alpha=0.5)
    ax.set_ylabel("Max. GP")
    ax.set_xlabel("Map Position [centimorgan]")
    ax.set_title(f"{iid}, Chromosome {ch}")
    ax.set_ylim([0.42, 1.02])

    if len(savefolder)>0:
            savepath = os.path.join(savefolder, f"{iid}_ch{ch}.png")
            plt.savefig(savepath, bbox_inches ='tight', pad_inches = 0, dpi=400)
            print(f"Saved to {savepath}")

    plt.show()
