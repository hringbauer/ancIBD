"""
Function for running hapBLOCK on single chromsome, with all relevant keywords.
Function to run hapBLOCK on a full individual and full reference Dataset,
with all relevant keywords.
@Author: Harald Ringbauer, 2020
"""

from hypothesis import assume
import numpy as np
import multiprocessing as mp
import pandas as pd
import itertools as it
from time import time
import os
import sys
import h5py
import sys
import os
from main import HMM_Full  # To run the main plotting.
from ancIBD.plot.plot_posterior import plot_posterior # to plot the posterior.
from ancIBD.IO.h5_load import get_opp_homos_f

def hapBLOCK_chrom(folder_in="./data/hdf5/1240k_v43/ch", iids = ["", ""], 
                   ch=2, folder_out="", output=False, prefix_out="", logfile=False,
                   l_model="hdf5", IBD2=False, p_col="variants/AF_ALL", 
                   ibd_in=1, ibd_out=10, ibd_jump=500, min_cm1=6, min_cm2=2,
                   cutoff_post=0.99, max_gap=0.0075, snp_cm=0, maf=0.0, save=0):
    """Run IBD for pair of Individuals.
    folder_in: hdf5 path up to chromosome.
    iids: List of IIDs to compare [length 2]
    folder_out: Where to save the hapBLOCK output to
    min_cm: Minimal block length to call and save [cM]
    savepath: Where to save the IBD plot to.
    p_col: The dataset to use in hdf5 for der. AF. If default use p=0.5.
    If empyt use in sample AF.
    Return df_ibd, posterior, map, tot_ll"""
    iids = np.array(iids) # For better props when indexing
    assert(len(iids)==2) # Sanity Check of Input IIDs

    e_model = "haploid_gl2" if not IBD2 else "IBD2"
    h_model = "FiveStateScaled" if not IBD2 else "SevenStateScaled"
    t_model = "standard" if not IBD2 else "IBD2"
    p_model = "hapROH" if not IBD2 else "IBD2"

    h = HMM_Full(folder_in=folder_in, l_model=l_model, t_model=t_model, 
                     e_model=e_model, h_model = h_model, p_model=p_model,
                     output=output, load=True)
    h.t_obj.set_params(ibd_in = ibd_in, ibd_out = ibd_out, ibd_jump = ibd_jump)
    h.l_obj.set_params(iids=iids, ch=ch, p_col=p_col, maf=maf)

    if len(folder_out)>0:
        folder_out = h.prepare_path(folder_out, ch=ch, prefix_out=prefix_out, logfile=logfile)    

    h.p_obj.set_params(ch=ch, min_cm1=min_cm1, min_cm2=min_cm2, cutoff_post=cutoff_post, max_gap=max_gap, snp_cm=snp_cm, folder=folder_out, save=save)
    #post, r_vec, fwd, bwd, tot_ll = h.run_fwd_bwd()

    post, r_vec, bp_vec =  h.run_fwd_bwd(full=False)

    assert(len(r_vec) == len(bp_vec))

    df_ibd, _, _ = h.p_obj.call_roh(r_vec, bp_vec, post, iid1=iids[0], iid2=iids[1])

    
    if len(folder_out)>0:
        folder_out = h.prepare_path(folder_out, ch=ch, prefix_out=prefix_out, logfile=logfile)
        save_path = os.path.join(folder_out, f"ch{ch}.tsv")
        h.p_obj.save_ibd_df(df_ibd=df_ibd, save_path = save_path)
    return df_ibd, post, r_vec


def multi_run(fun, prms, processes = 4, output=False):
    """Implementation of running in Parallel.
    fun: Function
    prms: The Parameter Files
    processes: How many Processes to use"""
    if output:
        print(f"Running {len(prms)} total jobs; {processes} in parallel.")
    
    if len(prms)>1:
        if output:
            print("Starting Pool of multiple workers...")    
        with mp.Pool(processes = processes) as pool:
            results = pool.starmap(fun, prms)
    elif len(prms)==1:
        if output:
            print("Running single process...")
        results = fun(*prms[0])
    else:
        raise RuntimeWarning("Nothing to run! Please check input.")
    return results


def hapBLOCK_pair(folder_in="./data/hdf5/1240k_v43/ch", iids = ["", ""], 
                   chs=range(1,23), folder_out="", output=False, prefix_out="", logfile=False,
                   l_model="hdf5", IBD2=False, p_col="variants/AF_ALL", 
                   ibd_in=1, ibd_out=10, ibd_jump=500, min_cm1=6, min_cm2=2, 
                   cutoff_post=0.99, max_gap=0.0075, snp_cm=0, maf=0.0, save=0):
    # running time on each chromosome is minimal, so I don't think there is any need to parallelize across chromosomes.
    # it makes more sense to only parallize across pairs
    assert(len(iids) == 2)
    dfs = []
    for ch in chs:
        df, *_ = hapBLOCK_chrom(folder_in, iids, ch, folder_out, output, prefix_out, logfile, l_model, IBD2, p_col, ibd_in, ibd_out, ibd_jump, min_cm1, min_cm2, cutoff_post, max_gap, snp_cm, maf, save)
        dfs.append(df)
    path_ibd = os.path.join(folder_out, f"{min(iids[0], iids[1])}_{max(iids[0], iids[1])}_ibd.tsv")
    df = pd.concat(dfs, ignore_index=True)
    df.to_csv(path_ibd, sep="\t", index=False)
    if IBD2:
        IBD_total_length = 100*np.sum(df[df['segment_type']=='IBD1']['lengthM'])
        IBD2_length = 100*np.sum(df[df['segment_type']=='IBD2']['lengthM'])
        print(f'IBD1+IBD2 region total length: {IBD_total_length:.3f}')
        print(f'IBD2 region total length: {IBD2_length:.3f}')
        return IBD_total_length, IBD2_length
    else:
        IBD_total_length = 100*np.sum(df['lengthM'])
        print(f'total length of IBD: {IBD_total_length}')
        return IBD_total_length

def hapBLOCK_all(folder_in=None, iids=[], chs=range(1,23), folder_out="", output=False, prefix_out="", logfile=False,
                    l_model='hdf5', IBD2=True, p_col='variants/AF_ALL', ibd_in=1, ibd_out=10, ibd_jump=500, min_cm1=6, min_cm2=2,
                    cutoff_post=0.99, max_gap=0.0075, snp_cm=0, maf=0.0, save=3, nprocesses=4):
    """
    Run hapBLOCK on all pairs listed in $iids. If $iids is empty, then run hapBLOCK on all pairs in the hdf5 input.
    """
    if len(iids) == 0:
        with h5py.File(f'{folder_in}{chs[0]}.h5', 'r') as f:
            iids = f['samples'][:].astype('str')
    
    iids_dedup = set(iids)
    if len(iids_dedup) != len(iids): # Check whether duplicates
        raise RuntimeWarning("Duplicate IIDs detected!")
    print(f'Run hapBLOCK on {len(iids_dedup)} samples...')
    
    prms = [ [folder_in, [id1, id2], 
                   chs, os.path.join(folder_out, f'{min(id1, id2)}_{max(id1, id2)}'),
                    output, prefix_out, logfile,
                   l_model, IBD2, p_col, 
                   ibd_in, ibd_out, ibd_jump, min_cm1, min_cm2, 
                   cutoff_post, max_gap, snp_cm, maf, save] for id1, id2 in it.combinations(iids_dedup, 2)]
    
    # now running everything in parallel
    multi_run(hapBLOCK_pair, prms, processes=nprocesses)



def prep_param_list_chrom(folder_in, iids = [], ch=3,
                    folder_out="", output=True, logfile=False, prefix_out="default/",
                    l_model="hdf5", e_model="haploid_gl", h_model="FiveStateScaled", 
                    t_model="standard", p_model="hapROH", p_col="variants/AF_ALL", ibd_in=1, ibd_out=1, ibd_jump=500, min_cm=2,
                    cutoff_post=0.99, max_gap=0.0, save=0):
    """Prepare parameter lists for multirun of hapBLOCK_chrom. Ideal for multi-processing,
    as it gives a list of parameters - one for each iid pair."""
    params = [[folder_in, iid2, ch, folder_out, output, prefix_out, logfile, l_model, e_model,
              h_model, t_model, p_model, p_col, ibd_in, ibd_out, ibd_jump, min_cm, cutoff_post, max_gap, save] for iid2 in iids]
    assert(len(params[0])==20)
    return params

#################################################################################
#################################################################################


def hapBLOCK_chroms(folder_in="./data/hdf5/1240k_v43/ch", iids = [], run_iids=[],
                   ch=2, folder_out="", output=False, prefix_out="", logfile=False,
                   l_model="hdf5", e_model="haploid_gl", h_model="FiveStateScaled", 
                   t_model="standard", p_col="variants/AF_ALL", ibd_in=1, ibd_out=10, ibd_jump=400, min_cm=2,
                   cutoff_post=0.99, max_gap=0.0075):
    """Run IBD for list of Individuals, and saves their IBD csv into a single 
    output folder.
    folder_in: hdf5 path up to chromosome.
    iids: List of IIDs to load [k indivdiuals]
    run_iids: If given: list of IID pairs to run. If not run all pairs
    folder_out: Where to save the hapBLOCK output to
    min_cm: Minimal block length to call and save [cM]
    savepath: Where to save the IBD plot to.
    Return df_ibd, posterior, map, tot_ll"""
    ### Run all pairs if empty
    iids = np.array(iids) # For better props when indexing
    if not len(set(iids))==len(iids): # Check whether duplicates
        raise RuntimeWarning("Duplicate IIDs detected!")
    if len(run_iids)==0:
        run_iids = it.combinations(iids, 2)
        
    ### Load all the objects
    h = HMM_Full(folder_in=folder_in, l_model=l_model, t_model=t_model, 
                     e_model=e_model, h_model = h_model,
                     output=output, load=True)
    h.t_obj.set_params(ibd_in = ibd_in, ibd_out = ibd_out, ibd_jump = ibd_jump)
    h.l_obj.set_params(iids=iids, ch=ch, p_col=p_col)
    h.p_obj.set_params(ch=ch, min_cm=min_cm, cutoff_post=cutoff_post, max_gap=max_gap)
    
    ### Load all data
    htsl, p, r_vec, samples =  h.l_obj.load_all_data()
    
    ### Load transition matrix
    t_mat = h.t_obj.full_transition_matrix(r_vec, n=4, submat33 = h.submat33)
    
    ### loop over all Run Pair Individuals
    df_ibds = []
    for iid1,iid2 in run_iids:
        i1 = get_sample_index(samples, iid1)
        i2 = get_sample_index(samples, iid2) 
        idcs = [i1*2, i1*2+1, i2*2, i2*2+1] # Get the right indices
        e_mat =  h.e_obj.give_emission_matrix(htsl[idcs,:], p)
        post =  h.fwd_bwd(e_mat, t_mat, in_val =  h.in_val, 
                            full=False, output= h.output)

        df_ibd, _, _ = h.p_obj.call_roh(r_vec, post, iid1, iid2)
        df_ibds.append(df_ibd)
    
    df_ibds = pd.concat(df_ibds)
    
    if len(folder_out)>0:
        folder_out = h.prepare_path(folder_out, ch=ch, prefix_out=prefix_out, logfile=logfile)
        save_path = os.path.join(folder_out, f"ch{ch}.tsv")
        h.p_obj.save_ibd_df(df_ibd=df_ibds, save_path = save_path)

    return df_ibds

def get_sample_index(iids, sample):
    """Get Index of sample - check if really there"""
    idx = np.where(iids[:]==sample)[0]
    assert(len(idx)==1)
    return idx[0]

############################################################################
### Run and plot in one go - for one pair of iids.

def run_plot_pair(path_h5="/n/groups/reich/hringbauer/git/hapBLOCK/data/hdf5/1240k_v43/ch", 
                  iids = ["", ""], ch=2, xlim=[], folder_out="", 
                  plot=False, path_fig="", output=False, exact=True,
                  ibd_in=1, ibd_out=10, ibd_jump=400, min_cm=2, 
                  cutoff_post=0.99, max_gap=0.0075, 
                  l_model="hdf5", e_model="haploid_gl", h_model="FiveStateScaled",
                  p_col="variants/AF_ALL", t_model="standard",
                  title="", c="gray", c_hw="maroon", 
                  state=0, return_post=False, **kwargs):
    """Run and plot IBD for pair of Individuals.
    folder_out: Where to save the hapBLOCK output to
    iids: list of two iids [List of Length 2]
    path_fig: Where to save the IBD plot to [String]
    p_col: The dataset to use in hdf5 for der. AF. If default use p=0.5.
           If empty string use in sample AF.
    return_post: Whether to return posterior [Boolean]
    kwargs: Optional Keyword Arguments for Plotting (e.g. c_ibd)
    """
    assert(len(iids)==2) # Sanity Check of Input IIDs - as here it should be pairs
    
    df_ibd, post, r_vec = hapBLOCK_chrom(
           folder_in=path_h5, iids = iids,  ch=ch, folder_out=folder_out, 
           output=output, prefix_out="", logfile=False,
           l_model=l_model, e_model=e_model, h_model=h_model, 
           t_model=t_model, p_col=p_col, ibd_in=ibd_in, ibd_out=ibd_out, ibd_jump=ibd_jump, 
           min_cm=min_cm, cutoff_post=cutoff_post, max_gap=max_gap)
        
    if plot:
        if len(title)==0:
            title = f"ancIBD v0.5, {iids[0]} - {iids[1]}, Chr. {ch}"
            
        ### Load the data from the HDF5
        o_homos, m = get_opp_homos_f(iid1=iids[0], iid2=iids[1], 
                                     f_path=path_h5, ch=ch, exact=exact)
        print(f"Plotting {len(r_vec)} markers")
        plot_posterior(post=post, morgan=r_vec, df_ibd=df_ibd, 
                       het=o_homos, het_m=m, state=state,
                       min_cm=min_cm, title=title, xlim=xlim, show=True, 
                       savepath=path_fig, xlabel="Chromosome Position [cM]", **kwargs)
    if return_post:
        return post
