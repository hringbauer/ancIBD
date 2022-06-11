"""
Function to do statistics on individual IBD dataframes
In particular functions to get population level IBD rates,
normalized per pair.
@ Author: Harald Ringbauer, 2020, All rights reserved
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import itertools as it
from scipy.special import gammaincinv

def new_columns(df, df_meta, col="New Clade", col_new="", match_col="iid"):
    """Maps Entries from meta dataframe onto the IBD dataframe.
    Return modified dataframe"""
    if isinstance(col, str):
        col = [col]
     
    for c in col:
        if len(col_new)==0:
            col_new1=c
    
        dct = pd.Series(df_meta[c].values, index=df_meta[match_col]).to_dict()
        for i in range(1,3):    
            df[col_new1 + str(i)] =df["iid" + str(i)].map(dct)
            df[col_new1 + str(i)] =df["iid" + str(i)].map(dct)
    return df

#################################################################################
### Helper Stats to find relatives

def find_relatives(df, iid="", min_ibd=20, cm=20):
    """ Identify all relatives of a given sample"""
    idx = (df["iid1"]==iid) | (df["iid2"]==iid)
    df1 = df[idx]
    df1 = df1[df1[f"sum_IBD>{cm}"]>=min_ibd]
    df1 = df1.sort_values(by=f"sum_IBD>{cm}", ascending=False)
    df1 = df1.copy().reset_index(drop=True)
    return df1

def rc_date(age):
    """Ascertain whether age is rc.
    age: Can be Array"""
    rc = np.mod(age, 10)
    rc = rc>0
    return rc

def plot_age_diff(df,figsize=(8,8), title="", xlim=[-2000,2000], ylim=[10,5000], 
                  cs=["red", "green"], fs=14, yscale="log", rcdate=False):
    """Plot the Age Difference between two samples"""
    df["delta_age"] = df["age1"] - df["age2"]
    if not rcdate:
        rc = rc_date(df["age1"]) &  rc_date(df["age2"])
    else:
        rc = (df["RC_age_min1"]!=0) & (df["RC_age_min2"]!=0)
    
    colors = [cs[1] if r else cs[0] for r in rc]
    
    plt.figure(figsize=figsize)
    ax = plt.gca()
    ax.scatter(df["delta_age"], df["sum_IBD>20"], ec="k", lw=0.5, 
               color=colors)

    ax.set_yscale(yscale)
    ax.set_ylabel("Sum ROH>20 cM [cM]", fontsize=fs)
    ax.set_xlabel("Age Difference of pair [years]", fontsize=fs)
    
    #
    p1 = mpatches.Patch(color=cs[0], label="Missing RC Date")
    p2 = mpatches.Patch(color=cs[1], label='Both RC Date')
    ax.legend(handles=[p1, p2], fontsize=fs)
    
    if len(title)>0:
        ax.set_title(title)
    if len(xlim)>0:
        ax.set_xlim(xlim)
    if len(ylim)>0:
        ax.set_ylim(ylim)
    plt.show()

    
#################################################################################

def give_sub_df(df, pop1="La Caleta", pop2="La Caleta", col="clst", 
                output=True, exact=False):
    """Return sub dataframe where pair across pop1 and pop2"""
    if exact:
        idx1 = (df[col+"1"]==pop1) & (df[col + "2"]==pop2)
        idx2 = (df[col+"1"]==pop2) & (df[col + "2"]==pop1)
    else:
        idx1 = (df[col+"1"].str.contains(pop1) & df[col + "2"].str.contains(pop2))
        idx2 = (df[col+"1"].str.contains(pop2) & df[col + "2"].str.contains(pop1))
        
    idx = idx1 | idx2
    if output:
        print(f"Found: {np.sum(idx)} Pairs fitting in dataframe.\n")
    return df[idx]

def give_stats_cm_bin(df, cms = [8,12,16,20], binary=True, output=True):
    """Return counts of IBD in bins. 
    df: IBD dataframe df of pw. individuals
    cms: Which bins to look into
    binary: Only count existing
    If upper bound 0: Take infinite bin"""
    counts = np.zeros(len(cms)-1) # Where the counts will be stored
    n = len(df) # The length of df (pw. comparisons)
    
    for i in range(len(cms)-1):
        cm1, cm2 = cms[i], cms[i+1]
        if cm2==0: # if no upper bound
            col1 = f"n_IBD>{cm1}"
            cts = df[col1].values
        else:
            col1, col2 = f"n_IBD>{cm1}", f"n_IBD>{cm2}"
            cts = df[col1].values - df[col2].values
        if binary:
            cts = np.clip(cts,0,1) # Only count existence
        tot_cts = np.sum(cts)
        counts[i] = tot_cts
        frac = np.mean(cts)
        if output:
            print(f"#IBD in {cm1}-{cm2} cM: {tot_cts}/{n}: {frac:.2f}")
    return counts

def get_IBD_stats(df, pop1="", pop2="", col="clade", exact=False,
                  cms=[4,6,8,10,12], binary=True, output=False, a=0.05):
    """Get IBD fraction statistics.
    a: Signficance level
    Return fractions, confidence intervalls as well as 
    number of pairsise comparisons"""
    if len(pop1)>0 or len(pop2)>0:
        df1 = give_sub_df(df, pop1=pop1, pop2=pop2, col=col, 
                          output=output, exact=exact)
    else:
        df1 = df
    
    if len(df1)==0: # If not enough pairs for comparison
        k = len(cms)-1
        fracs, cis, n = [0]*k, np.zeros((k,2)), 0 # Empty Comparison
        if len(df1)==0 and output:
            print(f"Empty Comparison: {pop1}-{pop2}")
    else:
        counts = give_stats_cm_bin(df1, output=output, cms=cms, binary=binary)
        n = len(df1) # The number of pw. comparisons
        fracs = counts/n
        cis = get_ci_counts(counts, n, a=a)
    #stds = np.sqrt(counts)/n
    return fracs, cis, n

def get_IBD_stats_pops(df, pops1=[], pops2=[], col="clade", 
                       cms=[4,6,8,10,12], output=False, 
                       binary=True, exact=False, a=0.05):
    """Get IBD fraction statistics for list of pop pairs.
    Return lists of fractions, confidence intervalls as well as 
    number of pairwise comparisons.
    a: Significance Level"""
    assert(len(pops1)==len(pops2)) # Sanity Check
    fracss, ciss, ns = [], [], []
    
    for i in range(len(pops1)):
        fracs, cis, n = get_IBD_stats(df, pop1=pops1[i], pop2=pops2[i], col=col,
                                      binary=binary, exact=exact, cms=cms, output=output, a=a)
        fracss.append(fracs)
        ciss.append(cis)
        ns.append(n)
    return fracss, ciss, ns

def get_ci_counts(counts, n, a=0.05, minc=1e-4):
    """Get Confidence Intervalls from counts and
    trials. Return list of CIS (lentght 2 each)
    counts: Array of Counts
    n: Total number of trials
    a: Signficance level"""
    cis = []
    for k in counts:
        if k==0:  # In case of no Count
            c0=minc/n
        else:
            c0 = gammaincinv(k, 0.5 * a)     # Lowe Bound
        c1 = gammaincinv(k + 1, 1 - 0.5 * a) # Upper bound
        cis.append([c0/n,c1/n])
    cis = np.array(cis) # Transform to numpy
    return cis

### Calulate pairwise IBD for all pops
def ibd_stats_pop_pairs(df, pop_pairs, cms = [8, 12, 16, 20, 0],
                        col = "label_region", a=0.32, binary=True,
                        output=False, exact=True):
    """For a list of population pairs, get ibd summary statistics for all pairs of
    individuals.
    Return list of ractions, list of CIs, list of #comparisons"""
    fracs, cis, ns = [], [], []
    
    for p1,p2 in pop_pairs:    
        frac, ci, n  = get_IBD_stats(df, pop1=p1, pop2=p2, col=col, binary=binary,
                                         cms=cms, output=output, a=a, exact=exact)
        fracs.append(frac)
        cis.append(ci)
        ns.append(n)
    return fracs, cis, ns

def create_ibd_pop_pair_df(pop_pairs, ns, fracs, cis, 
                           cms=[8, 12, 16, 20, 0], savepath="", ):
    """Save Population Pair dataframe if needed"""
    df_ibd = pd.DataFrame({"s1": pop_pairs[:,0], "s2": pop_pairs[:,1]})
    df_ibd["n_pairs"] = ns

    for i in range(len(cms)-1):
        cm1, cm2 = cms[i], cms[i+1]
        df_ibd[f"{cm1}-{cm2}cm"] = fracs[:,i]
        df_ibd[f"{cm1}-{cm2}cm_ci1"] = cis[:, i, 0]
        df_ibd[f"{cm1}-{cm2}cm_ci2"] = cis[:, i, 1]

    ### Save the F_ST-Matrix as well as the populaition:
    if len(savepath)>0:
        df_ibd.to_csv(savepath, sep="\t", index=False) # Save the Values
        print(f"Saved {len(df_ibd)} pw. pop comparisons.")
    return df_ibd