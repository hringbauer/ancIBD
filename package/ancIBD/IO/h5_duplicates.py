"""
Functions to detect duplicate Samples in imputed HDF5 file
@ Author: Harald Ringbauer, 2023
"""

import numpy as np
import pandas as pd
import h5py
import itertools as it

def get_match_df(path_h5="./data/hdf5/1240k_v54.1/ch", ch=3,
                 iids=[], gp_cutoff=0.98):
    """Create diploid Match Rate for all pairs in list of iids"""

    iids1, iids2, ms, os = [],[],[],[]

    with h5py.File(f"{path_h5}{ch}.h5", "r") as f:
        for iid1,iid2 in it.combinations(iids, r=2): 

            m,o = get_fraction_identical(f, sample1=iid1, sample2=iid2, 
                                   gp_cutoff=gp_cutoff, output=False)# Load for Sanity Check. See below!
            iids1.append(iid1)
            iids2.append(iid2)
            ms.append(m)
            os.append(o)

    dft = pd.DataFrame({"iid1":iids1,"iid2":iids2, "frac_match":ms, "frac_snps":os})
    return dft

def get_fraction_identical(f, sample1="SUC006", sample2="R26.SG", 
                           gp_cutoff=0.98, output=False):
    """Get Fraction of Identical Genotype Configurations.
    Return Fraction same IID, and fraction SNPs both IIDs tested"""
    if output:
        print(f"Running {sample1}-{sample2}")
    j1 = get_exact_idx_iid(f, sample1)[0]
    j2 = get_exact_idx_iid(f, sample2)[0]
    
    idx1 = get_markers_good(f, j1, cutoff=gp_cutoff, output=output)
    idx2 = get_markers_good(f, j2, cutoff=gp_cutoff, output=output)
    idx = (idx1 & idx2)
    snp_frac = np.mean(idx)
    if output:
        print(f"Filtering to common GP variants: {snp_frac:.3f}x")
    
    gt1 = f["calldata/GT"][:, j1, :][idx,:]
    gt2 = f["calldata/GT"][:, j2, :][idx,:]
    g1, g2 = np.sum(gt1, axis=1), np.sum(gt2, axis=1)
    frac_same = np.mean(g1 == g2)
    return frac_same,snp_frac

def get_exact_idx_iid(f, sample, unique=True):
    """Return Index of sample samples in hdf5 f"""
    samples = pd.Series(f["samples"][:].astype("str"))
    idx = samples==sample
    idx = np.where(idx)[0]
    assert(len(idx)==1)
    return idx

def get_markers_good(f, j, output=True, cutoff=0.99):
    """Get markers"""
    m = np.max(f["calldata/GP"][:,j,:], axis=1)
    idx = (m>cutoff)
    if output:
        c1 = np.mean(m>cutoff)
        print(f"Filtering to {cutoff} GP variants: {c1:.3f}x")
    return idx