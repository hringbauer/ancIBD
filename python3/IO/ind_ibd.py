"""
Class for post-processing pw. IBD list into a single 
summary dataframe for each pair of individuals, for various 
IBD length classes
@ Author: Harald Ringbauer, 2020
"""

import numpy as np
import pandas as pd
import os as os


def filter_ibd_df(df, min_cm=4, snp_cm=60, output=True):
    """Post Process ROH Dataframe. Filter to rows that are okay.
    min_cm: Minimum Length in CentiMorgan
    snp_cm: How many SNPs per CentiMorgan"""

    # Filter for Length:
    length_okay = (df["lengthM"] * 100) > min_cm
    
    # Filter for SNP Density:
    densities = df["length"] / (df["lengthM"] * 100)
    densities_ok = (densities > snp_cm)
    df["SNP_Dens"] = densities
    
    df = df[densities_ok & length_okay].copy()
    
    if output==True:
        print(f"> {min_cm} cM: {np.sum(length_okay)}/{len(length_okay)}")
        print(f"Of these with suff. SNPs per cM> {snp_cm}: \
              {np.sum(densities_ok & length_okay)}/{np.sum(length_okay)}")   
    return df

def roh_statistic_df(df, min_cm=0, col_lengthM="lengthM"):
    """Gives out Summary statistic of ROH df"""
    if len(df)==0:   # If no ROH Block found
        max_roh, sum_roh, n_roh = 0, 0, 0
    
    else:
        idx = df["lengthM"]>min_cm/100.0 # convert to Morgan
        if np.sum(idx)==0:
            max_roh, sum_roh, n_roh = 0, 0, 0
            
        else:
            l = df["lengthM"][idx]
            max_roh = np.max(l)
            sum_roh = np.sum(l)
            n_roh = len(l) 
    return sum_roh, n_roh, max_roh

def roh_statistics_df(df, min_cms=[8, 12, 16, 20], col_lengthM="lengthM"):
    """Gives out IBD df row summary statistics.
    Return list of sum_roh, n_roh, max_roh for each of them [as list]
    min_cm: List of minimum IBD lengths [in cM]"""  
    res = [roh_statistic_df(df, c, col_lengthM) for c in min_cms]
    return res

def create_ind_ibd_df(path_ibd = "/n/groups/reich/hringbauer/git/yamnaya/output/ibd/v43/ch_all.tsv",
                      min_cms = [8, 12, 16, 20], snp_cm = 220, min_cm = 6, sort_col = -1,
                      savepath="", output=True):
    """Create dataframe with summary statistics for each individual.
    Return this novel dataframe in hapROH format [IBD in cM]
    path_ibd: What ibd file to load.
    min_cms: What IBD lengths to use as cutoff in analysis [cM].
    snp_cm: Minimum Density of SNP per cM of IBD block.
    sort_col: Which min_cms col to use for sort. If <0 no sort conducted."""
    df_ibd = pd.read_csv(path_ibd, sep="\t")
    df_ibd = filter_ibd_df(df_ibd, min_cm=min_cm, snp_cm=snp_cm, output=output) # Be aggressive here for good set for relatedness. Cutoff the slack
    
    if output:
        print(df_ibd["ch"].value_counts())

    #df_ibd = df_ibd.sort_values(by = ["iid1", "iid2", "ch"]) # Sort so it is easier to split up
    grouped = df_ibd.groupby(['iid1','iid2'])

    df_res = pd.DataFrame(grouped.groups.keys())
    df_res.columns = ["iid1", "iid2"]

    ### Create the actual statistics
    stats = [roh_statistics_df(df, min_cms=min_cms) for _, df in grouped]
    stats = np.array(stats)

    df_res["max_IBD"] = stats[:, 0, 2] * 100

    ### Add Values for varying cM cutoffs:
    for j, m_cm in enumerate(min_cms):
        df_res[f"sum_IBD>{m_cm}"] = stats[:,j, 0] * 100
        df_res[f"n_IBD>{m_cm}"] = stats[:,j,1]

    if sort_col>=0:
        df_res = df_res.sort_values(by=f"sum_IBD>{min_cms[sort_col]}", ascending=False)  # Sort output by minimal Cutoff  
        
    if len(savepath) > 0:
        df_res.to_csv(savepath, sep="\t", index=False)
        print(f"Saved {len(df_res)} individual IBD pairs to: {savepath}")

    return df_res