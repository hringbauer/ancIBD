"""
Function to estimate Ne from IBD results
@ Author: Harald Ringbauer, 2026
"""

import numpy as np
import pandas as pd
from TTNe.analytic import singlePop_2tp_given_Ne_negLoglik
from ancIBD.ibd_stats.ibd_sites_stats import get_ibd_stats_unrelated

def fit_Ne_ibd(df, n_pairs, bins = np.arange(8, 30, 0.25), Ne_list = np.arange(100, 20000, 100)):
    """Fit Ne for IBD data frame using ML approach.
    Return Ne (maxL approach), lower CI, and upper CI (CIs based on maxLL - 1.92)
    df: IBD Dataframe of all IBD segments
    n_pairs: Number of Pairs (of IIDs)
    Ne_list: List of Ne values to try
    """
    binMidpoint = (bins[1:] + bins[:-1])/2
    histo, _ = np.histogram(100*df['lengthM'], bins=bins)
    
    ### Hard Coded Human Genome Parameters
    ch_len_dict = {1:286.279, 2:268.840, 3:223.361, 4:214.688, 5:204.089, 6:192.040, 7:187.221, 8:168.003, 9:166.359, \
        10:181.144, 11:158.219, 12:174.679, 13:125.706, 14:120.203, 15:141.860, 16:134.038, 17:128.491, 18:117.709, \
        19:107.734, 20:108.267, 21:62.786, 22:74.110}
    G = np.array([val for k, val in ch_len_dict.items()]) # Total Genom length

    mat = np.zeros(len(Ne_list))

    for i, Ne in enumerate(Ne_list):
        mat[i] = -singlePop_2tp_given_Ne_negLoglik(Ne, histo, binMidpoint, G, 0, 0, 4*n_pairs, [(0,0), (0,0)])

    maxll = np.max(mat)
    i1 = np.where(mat == maxll)
    i1 = i1[0][0]

    CIregionNe = np.where(mat >= maxll - 1.92)
    Ne_CI = Ne_list[CIregionNe]
    
    print(f'MLE Ne: {Ne_list[i1]} (95% CI: {np.min(Ne_CI)}-{np.max(Ne_CI)})')
    return Ne_list[i1], np.min(Ne_CI), np.max(Ne_CI)

def fit_Ne_IBD_df(df_ibd_ind, df_ibd, df_meta, iids, min_cm=50, 
                  bins=np.arange(8, 30, 0.5), Ne_list=np.arange(1e4, 1e6, 1e3)):
    """Fit Ne (point estimated and 95% CI) from IBD result dataframes.
    Return Ne, Ci low, CI high.
    Input: df_ibd_ind: Standard IBD summary dataframe (one row one pair)
    df_ibd: Standard IBD segment dataframe (one row one IBD)
    iids: Which individuals to fit.
    min_cm: Sum of IBD >12cm long above which individual pairs get filtered out [in cM"""

    dft, n_pairs = get_ibd_stats_unrelated(df_ibd_ind, df_ibd, df_meta, iids1=iids, iids2=iids, min_cm=min_cm)
    maxL, ci1, ci2 = fit_Ne_ibd(dft, n_pairs, bins = bins, Ne_list = Ne_list)
    return maxL, ci1, ci2