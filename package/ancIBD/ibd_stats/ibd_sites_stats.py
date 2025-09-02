# Various Functions to depict IBD sharing between sites
### Harald Ringbauer, Mar 2025

### Imports
import numpy as np
import pandas as pd
import itertools as it
import matplotlib.pyplot as plt

from TTNe.analytic import singlePop_2tp_given_Ne_negLoglik
from TTNe.analytic_multi import inferConstNe_singlePop_MultiTP
from hapsburg.PackagesSupport.roh_expectations import Expected_Roh
from scipy.sparse.linalg import expm_multiply

### Functions
def get_ibd_sites(df_ibd=[],  df_meta=[], site1="", site2="", output=True):
    """Return all IBD shared between site1 and site2.
    Return IBD dataframe as well as number of matching IIDs in the two sites in df_meta"""
    # Default of empty site2 is site1
    if len(site2)==0:
        site2 = site1
    
    iids1 = df_meta[df_meta['iid'].str.startswith(site1)]['iid'].to_list()
    iids2 = df_meta[df_meta['iid'].str.startswith(site2)]['iid'].to_list()
    n1, n2 = len(iids1), len(iids2)
    if output:
        print(f'Number of {site1} samples: {n1}')
        if site2!=site1: # Only Plot if second site different
            print(f'Number of {site2} samples: {n2}')
    
    dft = get_ibd_iids(df_ibd, iids1, iids2) # Get matching IBD
    
    return dft, n1, n2

def get_ibd_iid_lists(df_ibd=[], iids1=[], iids2=[], output=True):
    """Return all IBD shared between list of iids1 and iids2.
    Return IBD dataframe as well as number of matching IIDs in the two sites in df_meta"""
    # Default of empty site2 is site1
    if len(iids2)==0:
        iids2 = iids1

    n1, n2 = len(iids1), len(iids2)
    if output:
        print(f'Number of sample 1: {n1}')
        if set(iids1)!=set(iids2):
            print(f'Number of sample 2: {n2}')
    
    dft = get_ibd_iids(df_ibd, iids1, iids2) # Get matching IBD
    return dft, n1, n2

def get_ibd_iids(df_ibd, iids1=[], iids2=[]):
    """Get subset of all IBDs between iids1 (list) and iids2 (list).
    Return sub dataframe of df_ibd"""
    # select rows whose iid1 column is in the list of ETA_iid and iid2 column is in the list of PHA_iid
    # or vice versa
    subset1 = df_ibd[df_ibd['iid1'].isin(iids1) & df_ibd['iid2'].isin(iids2)]
    subset2 = df_ibd[df_ibd['iid1'].isin(iids2) & df_ibd['iid2'].isin(iids1)]
    dft = pd.concat([subset1, subset2])
    dft = dft.drop_duplicates(keep="first") # Drop Duplicate Entries. To deal with overlapping iids
    return dft

def find_relatives(df, col_rel="sum_IBD>12", min_cm=200):
    """Return all pairs of related IIDs in IBD dataframe df as well as nr of relatives"""
    nrelatives = np.sum(df[col_rel] > min_cm)
    print(f'Number of relatives: {nrelatives}')
    
    relatives = [[row['iid1'], row['iid2']] for index, row in df.iterrows() if row[col_rel] > min_cm]
    return nrelatives, relatives

def filter_ibd(df, relatives=[]):
    """Filter IBD df for all pairs in pairs (list of 2 Element lists)"""
    for relative in relatives:
        df = df[~((df['iid1'] == relative[0]) & (df['iid2'] == relative[1]))]
        df = df[~((df['iid1'] == relative[1]) & (df['iid2'] == relative[0]))]
    return df

def get_iids_in_meta(iids=[], df_meta=[]):
    """Subset to iids in Meta File.
    Return list of those iids and its length"""
    if len(df_meta)==0:
        iids_m = iids
    else:
        idx =  df_meta['iid'].isin(iids)
        iids_m = df_meta["iid"][idx].values
    return iids_m, len(iids_m)

def get_ibd_stats_unrelated(df_ibd_ind, df_ibd, df_meta, site1="", site2="", 
                            iids1=[], iids2=[], col_rel="sum_IBD>12", min_cm=200, output=True):
    """Get IBD dataframe and number tested between site1 and site2
    Return an IBD dataframe and pair number.
    If df_ibd_ind given, filter close relatives based on 
    its relative column using min_cm as cutoff.
    If iids1 given, use the iids1 and iids2 lists of individuals,
    otherwise the site1 and site2 values
    If iids2 or sits2 is empty, run a within-sample IBD comparison
    df_meta: Meta file of every individual ran for IBD
    df_ibd: File of all IBD segments
    df_ibd_ind: File of all IBD summary stats per pair"""
    
    ### Case 1 no iids given: Get IBD from matching Pandora Entries (from Meta)
    if len(iids1)==0:  
        dft, n1, n2 = get_ibd_sites(df_ibd, df_meta, site1=site1, site2=site2, output=output)
        dft_ind, _, _ = get_ibd_sites(df_ibd_ind, df_meta, site1=site1, site2=site2, output=False)
        nrelatives, relatives = find_relatives(dft_ind, min_cm=min_cm)

        ### Calculate the number of possible pairs
        if (site1==site2) or (len(site2)==0):
            n_pairs = int(n1 * (n1-1) / 2) # int just for printing below
        else:
            n_pairs = int(n1 * n2)

    ### Case 2: Get IBD for matching iids
    elif len(iids1)>0:
        if len(iids2)==0:
            iid2 = iids1
        ### Get iids and ibd segments matching meta file
        iids1, n1 = get_iids_in_meta(iids1, df_meta)
        iids2, n2 = get_iids_in_meta(iids2, df_meta)
        dft, _, _ = get_ibd_iid_lists(df_ibd, iids1=iids1, iids2=iids2)

        ### Find relatives if IBD summary stats given
        if len(df_ibd_ind)>0:
            dft_ind, _, _ = get_ibd_iid_lists(df_ibd_ind, iids1=iids1, iids2=iids2, output=False)
            nrelatives, relatives = find_relatives(dft_ind, min_cm=min_cm)
        else:
            nrelatives, relatives = 0, []

        ### Calculate the number of possible pairs
        if (set(iids1)==set(iids2)):
            n_pairs = int(n1 * (n1-1) / 2) # int just for printing below
        else:
            n_pairs = int(n1 * n2)

    #assert(n_pairs>0) # Sanity Check (for n=1 wrong)
    
    ### Remove Relatives
    dft=filter_ibd(dft, relatives) # Filter out the relative pairs
    n_pairs1 = n_pairs - nrelatives # Substract related indvidiuals
    print(f"Filtered to {n_pairs1}/{n_pairs} non-related iids ({col_rel}<{min_cm})")
    
    return dft, n_pairs1

def plot_ibd_2sites(df_ibd, df_ibd_ind, df_meta,
    site1="", site2="", iids1=[], iids2=[], min_cm=200,
    figsize = (6,6),
    bins = np.arange(8, 30, 1),
    Ne_plot = [4000, 35000], dts_plot=[], dts_a=1, 
    c_dts= ["#fde725", "#5ec962", "#21918c", "#3b528b", "#440154", "k"],
    c_plot = ["#5ec962", "#440154"], alpha=0.6, xlim_plot=[], 
    ylim_plot=[], yscale="log", 
    output_level=1, savepath=""):
    """Plot IBD across two sites together with predictions.
    If a site label is given, use it to pull its IIDs.
    If iids1 and iids2 list is given, use the iids.
    site1 and site2: First three letters of site (Pandora=specific)
    min_cm: Cutoff of summed IBD [in cm] of relative pairs to filter
    df_ibd: IBD dataframe, each row one IBD segment
    df_ibd_ind: IBD summary stats dataframe, each row one pair.
    If given, it is used to filter out relatives.
    Ne_plot: Diploid population size: Plot analytical predictions of those.
    Use hard-coded human chromosome lengths for calculations
    xlim_plot: [StartCM, EndCM] for plotting the IBD histogram [list].
    Default setting is tom in and max of bins
    dts_plot: Time differences of sampling to calculate and plot [list]
    Use IBD in group1 as dt=0 starting point and evolve IBD sharing forward.
    c_dts: Colors of dt calculations [list]. Matches dts_plot
    dts_a: Admixture Fraction to use when plotting dts.
    output_level: Detail of output to print. Default is 1.
    2 means all output is printed"""
    ### Pre-Process Input:
    if len(xlim_plot)==0:
        xlim_plot = [np.min(bins), np.max(bins)]
        
    ### Filter to the relevant pairs
    df_ibd1, npairs_1 = get_ibd_stats_unrelated(df_ibd_ind, df_ibd, df_meta, site1=site1, site2=site1, 
                                                iids1=iids1, iids2=iids1, min_cm=min_cm)
    df_ibd2, npairs_2 = get_ibd_stats_unrelated(df_ibd_ind, df_ibd, df_meta, site1=site2, site2=site2, 
                                                iids1=iids2, iids2=iids2, min_cm=min_cm)
    df_ibda, npairs_across = get_ibd_stats_unrelated(df_ibd_ind, df_ibd, df_meta, site1=site1, site2=site2, 
                                                     iids1=iids1, iids2=iids2, min_cm=min_cm, output=False)
    
    ### Plot the figure
    plt.figure(figsize=figsize)
    ax=plt.gca()

    ax.hist(100*df_ibd1['lengthM'], bins=bins, alpha=alpha, label=site1, edgecolor="gray", color="blue",
             weights=np.ones(len(df_ibd1))/npairs_1)
    ax.hist(100*df_ibd2['lengthM'], bins=bins, alpha=alpha, label=site2, edgecolor="gray", color="violet",
             weights=np.ones(len(df_ibd2))/npairs_2)
    ax.hist(100*df_ibda['lengthM'], bins=bins, alpha=alpha, label='Between', edgecolor="gray", color="lime",
                weights=np.ones(len(df_ibda))/npairs_across)

    ### Plot the Expected IBD
    e_roh = Expected_Roh()  
    binmidpoint = (bins[1:] + bins[:-1])/2
    binwidth = bins[1] - bins[0]

    for N, c in zip(Ne_plot, c_plot):
        y = e_roh.roh_pdf_allchr_N(binmidpoint/100, N=2*N )*4*binwidth/100
        ax.plot(binmidpoint, y, color=c, linestyle='dashed')
        ax.text(binmidpoint[0], y[0], f"Ne={N}", color='black', fontsize=12)

    if len(dts_plot)>0:
        if output_level==2:
            print(f"Plotting time differences: {dts_plot}")
        t0, te = dts_plot[0], dts_plot[-1]
        num_t = len(dts_plot)

        ### Use 0.1cm bins (more accurate for calculation)
        b1=np.linspace(0,100,1001)*0.01 # in 0.1cm bins from 0 to 100
        p0 = get_ibd_sharing_prob(df_ibd1, bins=b1, n_pairs=npairs_1)
        
        ### Calculate IBD sharing forward
        P = calc_IBD_decay_delta_t(b1, p0, t0=t0, te=te, num_t=num_t, a=dts_a)

        ### Add bins of 10 length bins (0.1 ->1cm bins here)
        P10 = P.reshape(P.shape[0],-1, 10).sum(axis=2)

        ### Plot the dt histograms
        x_plot = b1[:-1:10] * 100 # Convert to centimorgan
        t_plot = np.arange(num_t)

        for i in t_plot:
            ax.plot(x_plot+x_plot[1]/2, P10[i,:], color=c_dts[i], 
                    label=f"Predicted dt={dts_plot[i]}", ls="-", marker="", ms=5)

    ax.legend(loc='upper right', title="IBD Sharing")
    ax.set_xlabel('IBD segment length (cM)')
    ax.set_ylabel('Average # of IBD segments per pair')
    ax.set_xlim(xlim_plot)
    if len(ylim_plot)>0:
        ax.set_ylim(ylim_plot)
    ax.set_yscale(yscale)
    
    if len(savepath)>0:
        plt.savefig(savepath, bbox_inches ='tight', pad_inches = 0, dpi=400)
        print(f"Saved to {savepath}")
    
    plt.show()
    
def plot_ibd_single_site(df_ibd, df_ibd_ind, df_meta,
    site="", iids=[], min_cm=200, title="",  
    figsize = (6,6), npairs_force=0,
    bins = np.arange(8, 30, 1), ax=0,
    Ne_plot = [4000, 35000], # Diploid Pop Size
    c_plot = ["#5ec962", "#440154"], alpha=0.6, c="sandybrown",
    legend = True, yscale="log", savepath="", show=True):
    """Plot IBD Histogram of a single site.
    iids: If given, use iid list directly instead of site.
    npairs_force: If given do not calculate pair number but set this value"""
    if ax==0:
        plt.figure(figsize=figsize)
        ax=plt.gca()
        
    df_ibd, npairs = get_ibd_stats_unrelated(df_ibd_ind, df_ibd, df_meta, site1=site, site2=site, 
                                             iids1=iids, iids2=iids, min_cm=min_cm)
    if npairs_force>0:
        npairs=npairs_force
    
    ax.hist(100*df_ibd['lengthM'], bins=bins, alpha=alpha, label=site, edgecolor="gray", color=c,
             weights=np.ones(len(df_ibd))/npairs)
    ax.axhline(1/npairs, linestyle="-", label="One segment", zorder=0, color="lightgray") # Horizontal Line
    
    ### Plot the Expected IBD
    e_roh = Expected_Roh()  
    binmidpoint = (bins[1:] + bins[:-1])/2
    binwidth = bins[1] - bins[0]
    for N, c in zip(Ne_plot, c_plot):
        y = e_roh.roh_pdf_allchr_N(binmidpoint/100, N=2*N )*4*binwidth/100
        ax.plot(binmidpoint, y, color=c, linestyle='dashed')
        ax.text(binmidpoint[0], y[0], f"Ne={N}", color='black', fontsize=12)
    
    if legend:
        ax.legend(loc='upper right', title="IBD Sharing")
        
    if title:
        ax.set_title(title)
        
    ax.set_xlabel('IBD segment length (cM)')
    ax.set_ylabel('Average # of IBD segments per pair')
    ax.set_yscale(yscale)
    
    if len(savepath)>0:
        plt.savefig(savepath, bbox_inches ='tight', pad_inches = 0, dpi=400)
        print(f"Saved to {savepath}")
    
    if show:
        plt.show()   
        
        
def fit_Ne_ibd(df, n_pairs, bins = np.arange(8, 30, 0.25), Ne_list = np.arange(100, 20000, 100)):
    """Fit Ne and time gap. 
    df: IBD Dataframe of all IBD segments
    n_pairs: Number of Pairs (of IIDs)
    Ne_list: List of Ne values to try
    
    Return Ne and gap CI object"""
    binMidpoint = (bins[1:] + bins[:-1])/2
    histo, _ = np.histogram(100*df['lengthM'], bins=bins)
    
    ### Hard-Coded Human Genome Parameters
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
    return Ne_CI

###########################################
### Functions for Delta T sampling
def get_q_mat(x):
    """Get infinitesimal transition matrix.
    x: List of bins (edges). In Morgan."""
    n = len(x)-1
    w = x[1]-x[0] # Width of one Bin
    q = np.ones((n,n)) * 2 * w
    q = np.tril(q, k=-1)
    v = -np.sum(q, axis=1)/2 + w
    np.fill_diagonal(q,v)
    return q

def calc_IBD_decay_delta_t(x, p0, t0=0, te=50, num_t=6, a=1):
    """Calculate IBD sharing after dt steps.
    t step size is from t0 to te in num_t equal steps
    te is included
    x: Bins to use [array n+1, in morgan units]
    p0: Initial Probability or counts in IBD bins (at dt=0)"""
    assert(len(x) == len(p0)+1)

    q = get_q_mat(x) # Get Transition Rate Matrix
    # Exponentiate Rate Matrix and multiply with initial conditions
    P = expm_multiply(q.T, p0, start=t0, stop=te, num=num_t, endpoint=True)
    return P*a # Add admixture fraction

def get_ibd_sharing_prob(df_ibd, bins, n_pairs=1, col_L="lengthM"):
    """Get the probability of pw. IBD sharing in length bins.
    Return count array.
    df_ibd: Standard IBD dataframe (1 IBD segment per row) with lengthM column
    bins: Bins in Morgan
    n_pairs: Number of pairwise comparisons"""
    hist, _ = np.histogram(df_ibd[col_L], bins=bins)
    pr = hist/n_pairs
    return pr




