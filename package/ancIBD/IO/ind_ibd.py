"""
Class for post-processing pw. IBD list into a single 
summary dataframe for each pair of individuals, for various 
IBD length classes
@ Author: Harald Ringbauer, 2020
"""

import numpy as np
import pandas as pd
import os as os

#########################################################

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

def create_ind_ibd_df(ibd_data = "/n/groups/reich/hringbauer/git/yamnaya/output/ibd/v43/ch_all.tsv",
                      min_cms = [8, 12, 16, 20], snp_cm = 220, 
                      min_cm = 6, sort_col = -1,
                      savepath="", output=True):
    """Create dataframe with summary statistics for each individual.
    Return this novel dataframe in hapROH format [IBD in cM]
    ibd_data: If string, what ibd file to load. Or IBD dataframe.
    savepath: If given: Save post-processed IBD dataframe to there.
    min_cms: What IBD lengths to use as cutoff in analysis [cM].
    snp_cm: Minimum Density of SNP per cM of IBD block.
    sort_col: Which min_cms col to use for sort. If <0 no sort conducted."""
    if isinstance(ibd_data, str): 
        df_ibd = pd.read_csv(ibd_data, sep="\t")
    else:
        df_ibd = ibd_data
        
    ### Filter. Be aggressive here for good set for  relatedness. Cutoff the slack
    df_ibd = filter_ibd_df(df_ibd, min_cm=min_cm, snp_cm=snp_cm, output=output) 
    
    ### Fish out the no IBD data
    if len(df_ibd)==0:
        df_res = pd.DataFrame(columns=['iid1','iid2', "max_IBD"])
        for m_cm in min_cms:
            df_res[f"sum_IBD>{m_cm}"] = 0
            df_res[f"n_IBD>{m_cm}"] = 0
    
    else: # In case there are IBD
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
    
    ### Save if needed
    if len(savepath) > 0:
        df_res.to_csv(savepath, sep="\t", index=False)
        print(f"Saved {len(df_res)} individual IBD pairs to: {savepath}")

    return df_res

def ind_all_ibd_df(path_ibd = "/n/groups/reich/hringbauer/git/yamnaya/output/ibd/v43/ch_all.tsv",
                   col_lengthM="lengthM", snp_cm = 220, min_cm = 5,
                   output=True, sort=True, decimals=2, col_new="ibd", savepath=""):
    """Create dataframe with all IBD for each indivdiual pair
    Return this novel dataframe in hapROH format [IBD in cM]
    path_ibd: What ibd file to load.
    snp_cm: Minimum Density of SNP per cM of IBD block.
    sort: If True sort by longest IBD
    decimals: To how many decimals to round"""
    df_ibd = pd.read_csv(path_ibd, sep="\t")
    df_ibd = filter_ibd_df(df_ibd, min_cm=min_cm, snp_cm=snp_cm, output=output) # Cut the Slack
    
    if output:
        print(df_ibd["ch"].value_counts())

    grouped = df_ibd.groupby(['iid1','iid2']) # Group by Length of IBD

    df_res = pd.DataFrame(grouped.groups.keys())
    df_res.columns = ["iid1", "iid2"]

    ### Create the actual statistics
    stats = [ibd_lengths(df, col_lengthM=col_lengthM, string=True, sort=True, decimals=decimals) for _, df in grouped]
    df_res[col_new] = stats
    
    ### Sort by max IBD
    if sort:
        ibd_max = df_res[col_new].str.split(",").str[0].astype("float") # Get the longest IBD
        idx = np.argsort(-ibd_max) # Return sorted from highest to lowest
        df_res = df_res.iloc[idx,:]  # Sort output by minimal Cutoff  
        df_res=df_res.copy().reset_index(drop=True) # Reset Index
    
    ### Save if savepath given
    if len(savepath) > 0:
        df_res.to_csv(savepath, sep="\t", index=False)
        print(f"Saved {len(df_res)} individual IBD pairs to: {savepath}")

    return df_res

def ibd_lengths(df, col_lengthM="lengthM", string=True, sort=True, decimals=2, mpl=100):
    """Returns list of IBD lengths in IBD dataframe df. [in cM]
    string: If True - return comma seperated string.
    sort: Whether to sort IBD list"""  
    res = df[col_lengthM].values * mpl
    if sort:
        res = -np.sort(-res) # top to bottom
    if decimals>-np.inf:
        res = np.around(res, decimals=decimals)
    if string:
        res = ",".join(res.astype("str"))
    return res

def all_pairs_ibd(df_res, df_iid):
    """Return a new IBD dataframe with all possible IID pairs,
    set to 0 IBD if not in IBD dataframe
    df_iid: where to take iids from
    df_res: IBD dataframe (standard format)"""

    pairs = np.array(list(it.combinations(df_iid["iid"],2)))
    df_all = pd.DataFrame({'iid1':pairs[:,0], 'iid2':pairs[:,1]})
    df_all[df_res.columns[2:]]=0  ## Set new entries to 0
    df_all.set_index(["iid1", "iid2"], inplace=True)
    df_all.update(df_res.set_index(["iid1", "iid2"]))
    df_all.update(df_res.set_index(["iid2", "iid1"]))
    df_all.reset_index()
    return df_all


def combine_all_chroms(chs=[], folder_base="PATH/ch", 
                       path_save = "PATH/ch_all.tsv"):
    """Combine All Chromosomes.
    chs: Which Chromosomes to run [list]
    folder_base: Where to load from (path part up to including ch)
    path_save: Where to save the combined file to."""
    if len(chs)==0:
        chs = range(1,23)
                       
    df_ibds = []
    for ch in chs:
        path_load = f"{folder_base}{ch}.tsv"
        df = pd.read_csv(path_load, sep="\t")
        print(f"Chromosome {ch}; Loaded {len(df)} IBD")
        df_ibds.append(df)
    df_ibds = pd.concat(df_ibds)

    ### Save IBD
    df_ibds.to_csv(path_save, sep="\t", index=False) # Save the Values
    print(f"Saved {len(df_ibds)} IBD to {path_save}.")

