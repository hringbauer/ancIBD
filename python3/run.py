"""
Function for running hapBLOCK on single chromsome, with all relevant keywords.
Function to run hapBLOCK on a full individual and full reference Dataset,
with all relevant keywords.
@Author: Harald Ringbauer, 2020
"""

import numpy as np
import multiprocessing as mp
import pandas as pd
import sys as sys
sys.path.append("/n/groups/reich/hringbauer/git/hapBLOCK/python3/") 
from main import HMM_Full

def hapBLOCK_chrom(folder_in="./data/hdf5/1240k_v43/ch", iids = ["", ""], 
                   ch=2, folder_out="", output=False, prefix_out="", logfile=False,
                   l_model="hdf5", e_model="haploid_gl", h_model="FiveStateScaled", 
                   t_model="standard", ibd_in=1, ibd_out=10, ibd_jump=500, min_cm=2,
                   cutoff_post=0.99, max_gap=0.0075):
    """Run IBD for pair of Individuals.
    folder_in: hdf5 path up to chromosome.
    iids: List of IIDs to compare [length 2]
    folder_out: Where to save the hapBLOCK output to
    min_cm: Minimal block length to call and save [cM]
    savepath: Where to save the IBD plot to.
    Return df_ibd, posterior, map, tot_ll"""
    assert(len(iids)==2) # Sanity Check of Input IIDs
    h = HMM_Full(folder_in=folder_in, l_model=l_model, t_model=t_model, 
                     e_model=e_model, h_model = h_model,
                     output=output, load=True)
    h.t_obj.set_params(ibd_in = ibd_in, ibd_out = ibd_out, ibd_jump = ibd_jump)
    h.l_obj.set_params(iids=iids, ch=ch)
    h.p_obj.set_params(ch=ch, min_cm=min_cm, cutoff_post=cutoff_post, max_gap=max_gap)
    #post, r_vec, fwd, bwd, tot_ll = h.run_fwd_bwd()
    post, r_vec =  h.run_fwd_bwd(full=False)
    df_ibd, _, _ = h.p_obj.call_roh(r_vec, post)
    
    if len(folder_out)>0:
        folder_out = h.prepare_path(folder_out, iid=iids, ch=ch, prefix_out=prefix_out, logfile=logfile)
        h.p_obj.save_output(df=df_ibd, save_folder=folder_out) # r_map=[], post=[]
    return df_ibd, post, r_vec


def hapBLOCK_chroms(folder_in="./data/hdf5/1240k_v43/ch", iids = [], run_iids=[],
                   ch=2, folder_out="", output=False, prefix_out="", logfile=False,
                   l_model="hdf5", e_model="haploid_gl", h_model="FiveStateScaled", 
                   t_model="standard", ibd_in=1, ibd_out=10, ibd_jump=400, min_cm=2,
                   cutoff_post=0.99, max_gap=0.0075, processes=1):
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
    iids = np.array(iids) # For better seach props
    if len(run_iids)==0:
        run_iids = it.combinations(iids, 2)
        
    ### Load all the objects
    h = HMM_Full(folder_in=folder_in, l_model=l_model, t_model=t_model, 
                     e_model=e_model, h_model = h_model,
                     output=output, load=True)
    h.t_obj.set_params(ibd_in = ibd_in, ibd_out = ibd_out, ibd_jump = ibd_jump)
    h.l_obj.set_params(iids=iids, ch=ch)
    h.p_obj.set_params(ch=ch, min_cm=min_cm, cutoff_post=cutoff_post, max_gap=max_gap)
    
    ### Load all data
    t = time()
    htsl, p, r_vec, samples =  h.l_obj.load_all_data()
    
    e = time()
    print(f"Runtime Loading: {(e-t)} s")
    
    ### Load transition matrix
    t = time()
    t_mat = h.t_obj.full_transition_matrix(r_vec, n=4, submat33 = h.submat33)
    e = time()
    print(f"Runtime T Mat.: {(e-t)} s")
    
    ### loop over all Run Pair Individuals
    df_ibds = []
    for iid1,iid2 in run_iids:
        t = time()
        i1 = get_sample_index(samples, iid1)
        i2 = get_sample_index(samples, iid2) 
        idcs = [i1*2, i1*2+1, i2*2, i2*2+1] # Get the right indices
        e_mat =  h.e_obj.give_emission_matrix(htsl[idcs,:], p)
        e = time()
        print(f"Runtime Loading Emission Matrix: {(e-t)} s")
        
        t = time()
        post =  h.fwd_bwd(e_mat, t_mat, in_val =  h.in_val, 
                            full=False, output= h.output)
        e = time()
        print(f"Runtime FWD-BWD: {(e-t)} s")
        
        t = time()
        df_ibd, _, _ = h.p_obj.call_roh(r_vec, post, iid1, iid2)
        df_ibds.append(df_ibd)
        e = time()
        print(f"Runtime Postprocessing: {(e-t)} s")
    
    df_ibds = pd.concat(df_ibds)
    
    if len(folder_out)>0:
        folder_out = h.prepare_path(folder_out, iid=iids, ch=ch, prefix_out=prefix_out, logfile=logfile)
        h.p_obj.save_output(df=df_ibd, save_folder=folder_out) # r_map=[], post=[]

    return df_ibds

def get_sample_index(iids, sample):
    """Get Index of sample - check if really there"""
    idx = np.where(iids[:]==sample)[0]
    assert(len(idx)==1)
    return idx[0]