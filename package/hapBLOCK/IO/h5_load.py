import numpy as np
import pandas as pd
import socket as socket
import os as os
import sys as sys
import h5py


def get_idx_iid(f,sample, unique=True):
    """Return Index of sample samples in hdf5 f"""
    samples = pd.Series(f["samples"][:])
    idx = samples.str.contains(sample)
    idx = np.where(idx)[0]
    if unique:
        assert(len(idx)==1)
        return idx[0]
    else:
        return idx
    
def get_idx_iid_exact(f,sample, unique=True):
    """Return Index of sample samples in hdf5 f"""
    samples = pd.Series(f["samples"][:])
    idx = samples==sample
    idx = np.where(idx)[0]
    if unique:
        assert(len(idx)==1)
        return idx[0]
    else:
        return idx
    
def get_coverage(f, j):
    """Get Coverage of sample j in hdf5 f"""
    ads =  f["calldata/AD"][:,j,:2]
    ads[ads<0]=0
    cov = np.mean(ads)*2
    return cov

def get_markers_good(f, j, output=True, cutoff=0.99):
    """Get markers"""
    m = np.max(f["calldata/GP"][:,j,:], axis=1)
    idx = (m>cutoff)
    if output:
        c1 = np.mean(m>cutoff)
        print(f"Filtering to {cutoff} GP variants: {c1:.3f}x")
    return idx

def get_genos_pairs(f, sample1="SUC006", sample2="R26.SG", 
                    cutoff=0.98, output=True, phased=False, exact=False):
    """Return Genotypes and Map of pairs at intersection with GP>cutoff.
    phased: Whether to return [lx2] phased vector or [l] vetor of #derived.
    exact: Whether IID has to be an exact match"""
    if exact:
        j1 = get_idx_iid_exact(f, sample1)
        j2 = get_idx_iid_exact(f, sample2)
        
    else:
        j1 = get_idx_iid(f, sample1)
        j2 = get_idx_iid(f, sample2)
    
    idx1 = get_markers_good(f, j1, cutoff=cutoff, output=output)
    idx2 = get_markers_good(f, j2, cutoff=cutoff, output=output)
    idx = (idx1 & idx2)
    snp_frac = np.mean(idx)
    if output:
        print(f"Filtering to common GP variants: {snp_frac:.3f}x")
    
    m = f["variants/MAP"][:][idx]
    
    gt1 = f["calldata/GT"][:, j1, :][idx,:]
    gt2 = f["calldata/GT"][:, j2, :][idx,:]
    if not phased:
        gt1, gt2 = np.sum(gt1, axis=1), np.sum(gt2, axis=1)
    return gt1, gt2, m

def opp_homos(g1, g2):
    """Return opposing homozygotes"""
    o1 = (g1 == 0) & (g2 == 2)
    o2 = (g1 == 2) & (g2 == 0)
    return (o1 | o2)

def get_opp_homos_f(f_path="/n/groups/reich/hringbauer/git/hapBLOCK/data//hdf5/1240k_v43/ch",
                    iid1="SUC006", iid2="R26.SG", ch=3,
                    cutoff=0.99, output=True, exact=False):
    """Return opposing homozygotes boolean array and map array at intersection with
    GP>cutoff."""
    load_path = f_path + str(ch) + ".h5"
    with h5py.File(load_path, "r") as f:
        g1, g2, m = get_genos_pairs(f, sample1=iid1, sample2=iid2, 
                                cutoff=cutoff, output=output, exact=exact)
        o_homos = opp_homos(g1, g2)
    return o_homos, m