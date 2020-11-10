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
                   l_model="hdf5", e_model="haploid_gl", h_model="FiveStateFast", 
                   t_model="standard", ibd_in=1, ibd_out=1, ibd_jump=500, min_cm=2,
                   cutoff_post=0.99, max_gap=0.0):
    """Run IBD for pair of Individuals.
    folder_in: hdf5 path up to chromosome.
    iids: List of IIDs to compare [length 2]
    folder_out: Where to save the hapBLOCK output to
    min_cm: Minimal block length to call and save [cM]
    savepath: Where to save the IBD plot to"""
    assert(len(iids)==2) # Sanity Check of Input IIDs
    h = HMM_Full(folder_in=folder_in, l_model=l_model, t_model=t_model, 
                     e_model=e_model, h_model = h_model,
                     output=output, load=True)
    h.t_obj.set_params(ibd_in = ibd_in, ibd_out = ibd_out, ibd_jump = ibd_jump)
    h.l_obj.set_params(iids=iids, ch=ch)
    h.p_obj.set_params(ch=ch, min_cm=min_cm, cutoff_post=cutoff_post, max_gap=max_gap)
    post, r_vec, fwd, bwd, tot_ll = h.run_fwd_bwd()
    df_ibd, _, _ = h.p_obj.call_roh(r_vec, post)
    
    if len(folder_out)>0:
        folder_out = h.prepare_path(folder_out, iid=iids, ch=ch, prefix_out=prefix_out, logfile=logfile)
        h.p_obj.save_output(df=df_ibd, save_folder=folder_out) # r_map=[], post=[]