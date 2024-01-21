"""
Various Functions to prepare parameters for a batched run on a cluster
Splits up individuals into batches, run all pw. batches with ancIBD, and then collect results.
Functions here splits up input into batches, in a standardized way
@ Author: Harald Ringbauer, 2023
"""

import numpy as np
import pandas as pd
import os as os
import itertools as it

from ancIBD.IO.ind_ibd import create_ind_ibd_df, filter_ibd_df # Post-process IBD list to individual

#########################################################

def get_iids(path_meta="/n/groups/reich/hringbauer/git/yamnaya/data/meta_v2.tsv", min_snps = 600000):
    """Return list of iids to run"""
    df1 = pd.read_csv(path_meta, sep="\t")
    df1 = df1[df1["include"]==1].copy().reset_index(drop=True) # Only include unique,unrelated samples
    if min_snps>0:
        df2 = df1[df1["n_cov_snp"]>min_snps].copy()
        print(f"Loaded {len(df2)}/{len(df1)} samples with min snps: >{min_snps}")
    else:
        df2 = df1.copy()
        print(f"Loaded {len(df2)}/{len(df1)} samples with include==1")
    
    iids = df2["iid"].values
    return iids

def get_batch_idcs(i, batch_size):
    """Return the Index of batch and the within batch index"""
    iw = i % batch_size       # Within Index
    ib = int(i / batch_size)  # Batch index
    return iw, ib

def get_batch_pair_idx(i, batch_nr):
    """Return the Index of the two batches two run,
    using only triangular comparisons"""
    l1, l2 = np.tril_indices(n=batch_nr, k=0)
    return (l1[i], l2[i])

def create_savepath(folder_base="", ch=1, b1=1, b2=2, output=True):
    """Create the savepath in standardized output format"""
    savepath = os.path.join(folder_base, f"batch{b1}_{b2}", "ch" + str(ch) + ".tsv") # 
    if output:
        print(f"Created Save Path for IBD table: {savepath}")
    return savepath

def clean_double(run_idc=[], output=True):
    """Removes self comparisons as well as double comparisons.
    Input: list of pair-wise IIDs to run
    Return cleaned list."""
    idx_include = (run_idc[:,0] < run_idc[:,1])
    run_idc = run_idc[idx_include]
    if output:
        print(f"Filtered to {np.sum(idx_include)}/{len(idx_include)} pairs.")
    return run_idc

def get_idx_batch(b, batch_size):
    """Return all the indices of Indivdiuals of batch"""
    l,u = b*batch_size,(b+1)*batch_size
    return l,u

def get_unique_iid_pairs(iids, b1, b2, batch_size):
    """Get List of unique iid pairs to run.
    Return list of unique pairs, and list of all iids"""
    if b1>b2:
        l1,u1 = get_idx_batch(b1, batch_size)
        l2,u2 = get_idx_batch(b2, batch_size)
        iids1, iids2 = iids[l1:u1], iids[l2:u2]
        iids = np.concatenate((iids1, iids2)) # Append them to each other
        run_iids =  [[i1, i2] for i1 in iids1 for i2 in iids2] # All possible combinations. 
        
    elif b1==b2:
        l,u = get_idx_batch(b1, batch_size)
        iids = iids[l:u]
        run_iids = list(it.combinations(iids, 2))
    else:
        raise RuntimeError(f"Serious error, batch-index {b1}-{b2} not valid.")
    return iids, run_iids

def get_run_lists_batch(i = 67, k=3500, batch_size = 400, output=True):
    """Get IID list and Run lists for batches of samples.
    i: Run number
    k: number of total indivdiuals
    batch_size: number of individuals in one batch
    Return batch indices """
    batch_nr = int(k/batch_size)+1 # The number of batches to run
    pw_batch_nr = int(batch_nr * (batch_nr+1)/2)

    iw, ch = get_batch_idcs(i, pw_batch_nr) # Within batch number
    ch = ch+1 # For Autosome
    b1, b2 = get_batch_pair_idx(iw, batch_nr)

    if output:
        print(f"Run Nr: {i}")
        print(f"# of batches of individuals: {batch_nr}")
        print(f"# of pw. batches to run per chromosome: {pw_batch_nr}")
        print(f"Running Chromosome: {ch}")
        print(f"Running Batch Pair: {b1}-{b2}")
    return b1, b2, ch

def save_ibd_df(df_ibd, savepath, create=True):
    """Saves IBD Dataframe"""
    if create:
        folder = os.path.dirname(savepath)
        if not os.path.exists(folder):
            os.makedirs(folder)
    df_ibd.to_csv(savepath, sep="\t", index=False) # Save the Values
    print(f"Saved {len(df_ibd)} IBD blocks.")
    
def get_run_params_from_i(i, metapath="./data/iid_lists/iid_ibd_eurasia_v1.tsv", batch_size = 400, min_snps=0, 
                         output = True, folder_out = "/n/groups/reich/hringbauer/git/ibd_euro/output/ibd/v1/"):
    """Return the run parameters for run i
    min_snps: Minimum number of SNPs covered (for potential filtering on meta file)
    Returns iids, run_iids, and the output folder"""
    iids = get_iids(path_meta=metapath, min_snps=min_snps)
    b1, b2, ch = get_run_lists_batch(i = i, k=len(iids), batch_size = batch_size, output=output)
    iids, run_iids = get_unique_iid_pairs(iids, b1, b2, batch_size)
    path_ibd = create_savepath(folder_out, ch=ch, b1=b1, b2=b2, output=output)
    return iids, run_iids, ch, path_ibd

######################################################################################################
######################################################################################################
### Post-process batch runs

def join_chromosomes(base_path, chs=range(1,23), file_out="ch_all.tsv", output=True):
    """Join different Chromosomes together and save output.
    Return joined dataframe.
    file_out: If given, save the joint file with that name into the base_path folder"""
    df_ibds = []
    for ch in chs:
        path_save = os.path.join(base_path, f"ch{ch}.tsv")
        if not os.path.exists(path_save):
            print(f"Warning: {path_save} does not exist. Skipping!")
        else:
            df = pd.read_csv(path_save, sep="\t")
            if output:
                print(f"Chromosome {ch}; Loaded {len(df)} IBD")
            df_ibds.append(df)
    df_ibds = pd.concat(df_ibds) # Join together
    
    if len(file_out)>0:
        path_save = os.path.join(base_path, file_out)
        df_ibds.to_csv(path_save, sep="\t", index=False) # Save the Values
        if output:
            print(f"Saved {len(df_ibds)} IBD to {path_save}.")
    return df_ibds

def to_ind_df_batch(b1, b2, folder_out = "/n/groups/reich/hringbauer/git/ibd_euro/output/ibd/v1/",
                    chs=range(1,23), min_cms = [8, 12, 16, 20], snp_cm = 220, min_cm = 8, output=False):
    """Post-process a batch of individals.
    Returns individal IBD dataframe"""

    folder_batch = os.path.join(folder_out, f"batch{b1}_{b2}/")

    df_ibds = join_chromosomes(folder_batch, file_out="", chs=chs, output=output)

    df_res = create_ind_ibd_df(ibd_data=df_ibds,
                               min_cms = min_cms, snp_cm = snp_cm, min_cm = min_cm, sort_col = 0, 
                               output=output, savepath = "")
    return df_res

def to_ind_df_batches(batches=8, folder_out = "/n/groups/reich/hringbauer/git/ibd_euro/output/ibd/v1/",
                      chs=range(1,23), min_cms = [8, 12, 16, 20], snp_cm = 220, min_cm = 8, output=False,
                      savepath=""):
    """Runs multiple combinations of batches (wrapper of single batch function and then combines the processed batches. 
    Postprocess IBD to individal summary dataframe.
    Return merged IBD dataframe
    batches: If int: create all possible combinations. Otherwise needs to be array [n,2] of all
    pairs to run.
    savepath: If given, save IBD dataframe to there"""
    if isinstance(batches, int):
        batches = np.column_stack(np.tril_indices(n=batches, k=0))
    res = []
    for b1,b2 in batches:
        df_res = to_ind_df_batch(b1, b2, folder_out = folder_out, chs=chs,
                                 min_cms = min_cms, snp_cm = snp_cm, 
                                 min_cm = min_cm, output=output)
        res.append(df_res)
    df_res = pd.concat(res)
    return df_res 

def to_ibd_df_batches(batches=8, folder_out = "/n/groups/reich/hringbauer/git/ibd_euro/output/ibd/v1/",
                      chs=range(1,23), min_cms = [8, 12, 16, 20], snp_cm = 220, min_cm = 8, output=False,
                      savepath=""):
    """Runs multiple combinations of batches (wrapper of single batch function and then combines the processed batches. 
    Postprocess IBD to individal summary dataframe.
    Return merged IBD dataframe
    batches: If int: create all possible combinations. Otherwise needs to be array [n,2] of all
    pairs to run.
    savepath: If given, save IBD dataframe to there"""
    if isinstance(batches, int):
        batches = np.column_stack(np.tril_indices(n=batches, k=0))
    res = []
    
    for b1,b2 in batches:
        folder_batch = os.path.join(folder_out, f"batch{b1}_{b2}/")
        df_ibds = join_chromosomes(folder_batch, file_out="", 
                                   chs=chs, output=output)
        # Do Filtering here per batch to avoid giant file size
        df_ibds = filter_ibd_df(df_ibds, min_cm=min_cm, snp_cm=snp_cm, output=output) 
        res.append(df_ibds)
    df_res = pd.concat(res)
    return df_res 

def print_runid_missing(b = 1, folder_out = "", output=False):
    """Finds and prints indices of missing output (chXX.tsv) for batchwise runs.
    Return list of missing indices. Ideal for rerunning batch scripts.
    Uses C Indexing as would be used in submission script."""
    batches = np.column_stack(np.tril_indices(n=b, k=0))

    ls = []
    for ch in range(1,23):
        for b1,b2 in batches:
            folder_batch = os.path.join(folder_out, f"batch{b1}_{b2}/ch{ch}.tsv")
            exist = os.path.exists(folder_batch)
            if not exist:
                if output:
                    print(f"Missing! ch: {ch} batch: {b1}-{b2}")
                l = b * (b+1)/2 * (ch-1) + (b1+1) * b1/2 + b2 + 1 # The +1 is c indexing
                ls.append(str(int(l)))
    return ls

def find_output_missing(metapath="", folder_out="",
                        batch_size = 400, rge = [10, 20]):
    """Return List of all run nr.s that are missing.
    metapath: Path to .tsv of IIDs run for IBD screening [str]
    folder_out: Output folder [str]
    batch_size: How many individuals have been run per batch [int]"""
    
    df = pd.read_csv(metapath, sep="\t")
    k = len(df)

    run_nr_missing = []
    for i in range(rge[0], rge[1]):
        b1, b2, ch  = get_run_lists_batch(i, k=k, batch_size=batch_size, output=False)
        path = os.path.join(folder_out, f"batch{b1}_{b2}", f"ch{ch}.tsv")
        there = os.path.exists(path)

        if not there:
            run_nr_missing.append(i)

    print(f"Did not find output for: {len(run_nr_missing)} / {rge[1]-rge[0]} Runs")
    if len(run_nr_missing) == 0:
        print("Success! Everything found")
    return run_nr_missing

def get_batch_nr(n_iids, batchsize=400, n_chr=22):
    """Get the number of jobs to submit to cluster.
    Return n [int]"""
    n_batch = np.ceil(n_iids/batchsize)
    print(f"{n_iids} in batches of {batchsize}: {n_batch} Batches")
    n = int((n_batch * (n_batch+1) * n_chr)/2)
    print(f"Need {n} Submissions in total (0-{n-1})")