"""
Functions to QC imputed data (in hdf5 format) and select samples for IBD calling.
Mainly based on Leipzig big IBD run
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
    
##############
## Leipzig functions (moved here by Harald Jan 2025)

def get_h5_gp_stats(path_h5="", ch=3, het=False, gp_cutoff=0.99):
    """Return dataframe with coverage statistics of HDF5 file at 
    path_h5 and chromosome ch.
    het: Whether to calcualte het rate too."""
    with h5py.File(f"{path_h5}/ch{ch}.h5", "r") as f: # Load for Sanity Check. See below!
        print(list(f["variants"]))
        print(np.shape(f["calldata/GT"]))
        #gt = f["calldata/GT"][:]
        gp = f["calldata/GP"][:][:,:,:] # Fill in SNP number for TS
        samples = f["samples"][:].astype("str")

    ### Get Fraction of Genotypes imputed very well
    gp_max = np.max(gp, axis=2)
    gp = 0 # To free up memory
    gp_high = gp_max >= gp_cutoff
    f_gp_high = np.mean(gp_high, axis=0)

    ### Get missing data
    msg = np.isnan(gp_max)
    msg_avg = np.mean(msg, axis=0)
    
    ### Calculate het rate
    hets = 0
    if het:
        with h5py.File(f"{path_h5}/ch{ch}.h5", "r") as f:
            gt = f["calldata/GT"][:]
            gt2 = np.sum(gt, axis=2)
            hets = np.mean(gt2==1, axis=0)# Average Het Rate
            
    df = pd.DataFrame({"iid":samples, "frac_gp": f_gp_high, "frac_missing":msg_avg, "frac_het":hets})
    return df

def filter_df_suitable(df, frac_gp = 0.5, gp_cutoff=0.99,  
                       max_het=0.32, max_missing=0.1,
                       remove_iids=[]):
    """Filter Statistic df to samples suitable for IBD calling.
    Return new filtered df"""
    dft = df[df["frac_gp"]>frac_gp].copy().reset_index(drop=True)
    print(f"Filtering for maxGP {frac_gp}: {len(dft)}/{len(df)}")
    
    ### Filter for Heterozygosity
    idx = dft["frac_het"]<max_het
    dft = dft[idx].copy().reset_index(drop=True)
    print(f"Filtering for max_het {max_het}: {len(dft)}/{len(idx)}")
    
    ### Remove IIDs
    idx= dft["iid"].isin(remove_iids)
    dft = dft[~idx].copy().reset_index(drop=True)
    print(f"Filtering for IID: {len(dft)}/{len(idx)}")
    
    idx = dft["frac_missing"]>max_missing
    if np.sum(idx)>0:
        print(f"Warning: Samples with too much missing data!")
        dft = dft[~idx].copy().reset_index(drop=True)
    print(f"Suitable for IBD after full filtering: {len(dft)}/{len(df)}")
    return dft

def check_consistent_h5s(path_h5, chs=range(1,23)):
    """Check that h5 files are consistent"""

    sampless = []
    for ch in chs:
        with h5py.File(f"{path_h5}/ch{ch}.h5", "r") as f: # Load for Sanity Check. See below!
            samples = f["samples"][:].astype("str")
            sampless.append(samples)
            n = np.shape(f["calldata/GP"])[0]
            print(f"Chr. {ch}: {n} SNPs")
    there = np.all([np.array_equal(sampless[0], arr) for arr in sampless])
    
    if there:
        print(f"Green light! Found {len(sampless[0])} identical IIDs in all chromosome HDF5. ")
    else:
        print("Warning, Indivdiuals are not identical across chromosomes HDF")
        ls = [len(s) for s in sampless]
        print(ls)
        
def find_duplicates(path_h5, ch=3):
    """Find duplicate samples"""
    with h5py.File(f"{path_h5}/ch{ch}.h5", "r") as f: # Load for Sanity Check. See below!
        samples = f["samples"][:].astype("str")
        
    c = pd.value_counts(samples)
    dups = c[c>1].index.values
    
    if np.max(c)==1:
        print("Green Light! No duplicates")
    else:
        print(f"Warning, duplicate IIDs detected: \n{c[c>1]}")
    
    return dups, c
