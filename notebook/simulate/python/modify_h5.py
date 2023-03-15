"""
Class to load, modify and save HDF5 (simulated or empirical data).
Has functionality to down-sample and add errors to the data.
@ Author: Harald Ringbauer, 2020
"""

import allel
import h5py  # Python Package to do the HDF5.
import numpy as np
import pandas as pd
import socket
import os
import random


class ModifyHDF5Genotypes(object):
    """Class for Modifying HDF5 genotypes and
    saving new HDF5s. Can downsample/throw down error/Create Readcound.
    Plan: Also do contamination"""

    f = 0    # The hdf5 object to modify
    original_path = "" # Where to find the original HDF5
    save_path = ""  # Where to save the modified HDF5 to
    output = True # Whether to print any output
    gt_new = []

    def __init__(self, original_path="", save_path="", output=True):
        """pop_path: Where to load a HDF5 from
           save_path: Where to save the new HDF5 to"""
        self.output = output
        self.save_path = save_path
        
        if output == True:
            print("Heyho back old friend. I started running")
        
        if len(original_path)>0:
            self.original_path = original_path
            self.load_data()
        else:
            print("No HDF5 Loaded! Alarm. Alarm. Alarm.")

    def load_data(self, path=""):
        """Load the HDF5 Data"""
        if len(path)==0:
            path = self.original_path
        self.f = h5py.File(path, "r") # Load for Sanity Check. See below!
        
        if self.output == True:
            print("Loaded HDF5")
            print("Loaded %i variants" % np.shape(self.f["calldata/GT"])[0])
            print("Loaded %i individuals" % np.shape(self.f["calldata/GT"])[1])
            print(list(self.f["calldata"].keys()))
            print(list(self.f["variants"].keys()))
            #self.f["samples"] # Samples Vector
        
        ### Sanity Check whether both Genotypes are there and nothing else
        assert(np.min(self.f["calldata/GT"]) == 0)
        assert(np.max(self.f["calldata/GT"]) == 1)

    def save_data(self, gt, ad, ref, alt, pos, 
                  rec, allfreq, samples, path, gp=[],
                  compression="gzip", ad_group=True, gt_type="int8"):
        """Create a new HDF5 File with Input Data.
        gt: Genotype data [l,k,2]
        ad: Allele depth [l,k,2]
        ref: Reference Allele [l]
        alt: Alternate Allele [l]
        pos: Position  [l]
        m: Map position [l]
        allfreq: allele frequency [l]
        samples: Sample IDs [k].
        Save genotype data as int8, readcount data as int16.
        ad_group: whether to save allele depth
        gt_type: What genotype data type save"""

        l, k, _ = np.shape(gt)  # Nr loci and Nr of Individuals

        if os.path.exists(path):  # Do a Deletion of existing File there
            os.remove(path)

        dt = h5py.special_dtype(vlen=str)  # To have no problem with saving
        with h5py.File(path, 'w') as f0:
            ### Create all the Groups
            f_map = f0.create_dataset("variants/MAP", (l,), dtype='f')
            if ad_group:
                f_ad = f0.create_dataset("calldata/AD", (l, k, 2), dtype='int8', compression=compression)
            f_ref = f0.create_dataset("variants/REF", (l,), dtype=dt)
            f_alt = f0.create_dataset("variants/ALT", (l,), dtype=dt)
            f_pos = f0.create_dataset("variants/POS", (l,), dtype='int32')
            f_af = f0.create_dataset('variants/AF_ALL', (l,), dtype='f')
            f_gt = f0.create_dataset("calldata/GT", (l, k, 2), dtype=gt_type, compression=compression)
            if len(gp)>0:
                f_gp = f0.create_dataset("calldata/GP", (l, k, 3), dtype="f", compression=compression)     
            f_samples = f0.create_dataset("samples", (k,), dtype=dt)

            ### Save the Data
            f_map[:] = rec
            if ad_group:
                f_ad[:] = ad
            f_ref[:] = ref.astype("S1")
            f_alt[:] = alt.astype("S1")
            f_pos[:] = pos
            f_af[:] = allfreq
            f_gt[:] = gt
            if len(gp)>0:
                f_gp[:] = gp
            f_samples[:] = np.array(samples).astype("S10")
        if self.output == True:
            print(f"Successfully saved {k} individuals to: {path}")

    def create_error_gt(self, freq_flips=0.01, gp=True, cty=0.99):
        """Create Error on the HDF5 of genotypes.
        freq_flips: How often to do flip of genotyps"""
        f = self.f
        gt = f["calldata/GT"]
        switch = np.random.random(np.shape(gt)) < freq_flips
        if self.output:
            print(f"Swapping frac of SNPs: {np.mean(switch):.6f}")
        gt_new = (gt + switch) %2 # Switch the Genotypes
        
        if gp:
            gp = gp_from_gts(gt_new, cty=cty)  
        else:
            gp = []

        self.save_data(gt_new, f["calldata/AD"], f["variants/REF"][:], 
                       f["variants/ALT"][:], f["variants/POS"], 
                       f["variants/MAP"], f["samples"][:], 
                       self.save_path, gp=gp)
            
    def downsample_gt(self, frac=0.9, cty=0.99, shuffle_cm=0, error=0,
                      ad=True, gp=False, mult_alt=False, 
                      gt_type="int8", simulated_error=None, oneModern=False, compression="gzip"):
        """Downsample the HDF5 to fewer reads.
        Update also the recombination and position map if needed to remove missing values
        frac: To what fraction of markers one downsamples
        ad: Whether original HDF5 has AD field
        cty: Certainty of Genotype Probabilities
        error: If >0: Fraction of genotypes to shuffle (at random)
        shuffle_cm: If >0 shuffle cM with random exp. waiting times (shuffle_cm mean)
        mult_alt: Whether there are multiple alternative Allelels in the original HDF5
        simulated_error: A dictionary storing empirical imputation error. If this is provided,
        then frac and error will be ignored.
        oneModern: whether to mimic the scenario where one sample is ancient and the other is modern.
        """

        f = self.f
        gt = f["calldata/GT"][:] # Load everything in 1 go
        r_map = f["variants/MAP"][:]
        bps = f["variants/POS"][:]

        if simulated_error:
            frac = 1
            error = 0
        
        ### Downsample
        l,n,_ = np.shape(gt)
        if frac<1:
            survive = np.random.random(l) <= frac
            print(f"Fraction Loci surviving {np.mean(survive):.6f}")
            gt_new = gt[survive,:,:].astype(gt_type)
            r_map_new = r_map[survive]
        else:
            gt_new = gt
            r_map_new = r_map
            survive = np.ones(l, dtype="bool")
        
        if shuffle_cm>0:
            gt_new = shuffle_haplos(gt_new, r_map_new, scale_cm=shuffle_cm, oneModern=oneModern)
            
        if error>0:
            switch = (np.random.random(np.shape(gt_new)) < error) & (gt_new >= 0)
            gt_new = (gt_new + switch) %2 # Switch the Genotypes
            
        if ad:
            ad_new = f["calldata/AD"][survive,:,:]
        else:
            ad_new = np.zeros(np.shape(gt_new), dtype="int8")
            
        if gp and not simulated_error:
            gp = gp_from_gts(gt_new, cty=cty)
        elif gp and simulated_error:
            gt_new, gp = gp_from_empirical(gt_new, bps, simulated_error, oneModern=oneModern)
        else:
            gp = []
        
        ref_new = f["variants/REF"][survive]
        
        if mult_alt:
            alt_new = f["variants/ALT"][survive,0]   
        else:
            alt_new = f["variants/ALT"][survive]
        
        pos_new = f["variants/POS"][survive]
        allfreq_new = f['variants/AF_ALL'][survive]
        
        ### Downsample where needed  
        self.save_data(gt_new, ad_new, ref_new, alt_new, pos_new, r_map_new, allfreq_new,
                       f["samples"], self.save_path, gp=gp,
                       ad_group=ad, gt_type=gt_type, compression=compression)
        
    def generate_binomial_rc(self, mean_rc=1):
        """Generate Readcount Data from GT data.
        mean_rc: The Mean total Readcount per site"""
        
        f = self.f
        gt = f["calldata/GT"]
        
        ### Create the Poisson Readcounts with the right mean
        rc_full = poisson_readcounts(gt, mean_rc, output=self.output) 
        
        self.save_data(gt, rc_full, f["variants/REF"][:], f["variants/ALT"][:], f["variants/POS"], 
               f["variants/MAP"], f["samples"][:], self.save_path)
        
    def generate_lambda_rc(self, mean_rc = 1, norm_counts=True,
                           lambda_path = "./Data/1000Genomes/Coverage/mean_cov1240k_Marcus.csv"):
        """Generate Readcount Data from GT data.
        Use Table found at lambda_path for Lambdas 
        (relative. mean coverages, normed to 1 genome-wide)
        norm_counts: Whether to normalize on overlapping Readcounts"""
        
        df_lambda = load_lambda(lambda_path, output=self.output)  ### Load the Lambda Data
        
        f = self.f
        gt = f["calldata/GT"]
        l, n, _ = np.shape(gt)
        
        pos_f = f["variants/POS"][:]  # The Position of the Original 
        _, i1, i2 = np.intersect1d(pos_f, df_lambda["Pos"], return_indices=True)
        
        if self.output==True:
            print(f"Found {len(i1)} / {l} Loci in Lambda Table")
        
        lambdas = df_lambda["Lambda"].values[i2]
        if norm_counts == True: # Normalize to extracted lambdas
            lambdas = lambdas / np.mean(lambdas)
            
        mean_cov = lambdas * mean_rc  # Extract the Means that Intersect
        gt = gt[i1,:,:]  # Downsample to Loci intersecting the Lambda Table 
        
        ### Do the Binomial Readcount Sampling
        rc_full = poisson_readcounts(gt, mean_cov[:,None], output=self.output) 
        
        i1 = list(i1)  # So that it works with HDF5
        self.save_data(gt, rc_full, f["variants/REF"][i1], f["variants/ALT"][i1], f["variants/POS"][i1], 
               f["variants/MAP"][i1], f["samples"][:], self.save_path)
        
    def generate_ph(self, coverage = 1.0, error = 0.0):
        """Generate Pseudo-Haploid Data with fraction coverage sites covered,
        and then error thrown down.
        coverage: Fraction of sites covered
        error: Fraction of sites with error. If >0, flip error added at random"""
        
        f = self.f
        gt = f["calldata/GT"]
        l, _, _ = np.shape(gt)
        
        idx = np.random.random(l)<=coverage  # Which sites are covered
        gt = gt[idx, :, :]  # Extract downsampled SNPs
        
        switch = [0,]
        if error>0:
            switch = (np.random.random(np.shape(gt)) < error) & (gt >= 0)
            gt = (gt + switch) %2 # Switch the Genotypes
                
        if self.output:
            print(f"{np.sum(idx)} / {len(idx)} SNPs pseudohaploidized.")
            print(f"Added fraction errors to SNPs: {np.mean(switch):.6f}")
            print(f"Added sum errors: {np.sum(switch):.0f}")
        
        rc = pseudo_haploid(gt) # Generate Pseudo-Haploid Readcounts
        
        idx = np.array(idx)  # So that it works with HDF5 (Boolean Indexing)
        self.save_data(gt, rc, f["variants/REF"][idx], f["variants/ALT"][idx], f["variants/POS"][idx], 
                       f["variants/MAP"][idx], f["samples"][:], self.save_path)
        
    def copy_ibdinfo(self, load_path="", save_path="", file="ibd_info.csv"):
        """Copy in the IBD Info from folder of load path into folder of save_path.
        file: Which file to copy (ibd_info: default)"""
        if len(load_path) == 0:
            load_path = self.original_path
            
        if len(save_path) ==0 :
            save_path = self.save_path
            
        save_path = os.path.dirname(save_path) + "/" + file
        load_path = os.path.dirname(load_path) + "/" + file
        
        ### Copy the file
        command = f"cp {load_path} {save_path}"
        os.system(command) 
        
##########################################
#### Some Small Helper Functions

def load_lambda(loadpath, ch=3, output=True):
    """Load and return the Lambda Vector
    for Chromosome ch, and from path loadpath"""
    df_lambda = pd.read_csv(loadpath)
    mean = np.mean(df_lambda["Lambda"])
    assert(np.isclose(mean, 1))  # Sanity Check if Valid Lambda Vector
    l=len(df_lambda)
    df_lambda = df_lambda[df_lambda["Ch"]==ch]
    if output==True:
        print(f"Extracted {len(df_lambda)} / {l} Loci on Chr.{ch}")
    return df_lambda

def poisson_readcounts(gt, mean_rc, output=True):
    """Create and return Poisson Readcount array.
    gt: Underlying Genotype Matrix [l, n, 2]
    Return readcound array: [l, n, 2]"""
    l, n, _ = np.shape(gt)
    rc_tot = np.random.poisson(lam=mean_rc, size = (l,n))  # Draw Full Readcounts

    p = np.mean(gt, axis=2) # Get the Mean Allele Frequency per locus and individual
    assert(np.max(p)<=1) ### Sanity Check whether allele freqs are right
    assert(np.min(p)>=0)

    rc_der = np.random.binomial(n=rc_tot, p=p)  # The derived Readcount (Binomial Sampling)
    rc_ref = rc_tot - rc_der  # The Ref Readcount

    rc_full = np.stack([rc_ref, rc_der], axis=2)
    assert(np.shape(rc_full) == np.shape(gt))  # Check whether data was created properly

    if output == True:
        print(f"Mean Readcount: {np.mean(rc_tot):.4f}")
    
    return rc_full

def pseudo_haploid(gt):
    """Create and return Pseudo-Haploid Readcount array
    gt: Underlying Genotype Matrix [l, n, 2]
    Return readcound array: [l, n, 2]"""
    
    p = np.mean(gt, axis=2) # Get the Mean Allele Frequency per locus and individual
    
    rc_der = np.random.binomial(n=1, p=p)  # The derived Readcount (Binomial Sampling)
    rc_ref = 1 - rc_der  # The Ref Readcount
    rc_full = np.stack([rc_ref, rc_der], axis=2)
    
    assert(np.max(rc_full)<=1)
    assert(np.min(rc_full)==0)
    assert(np.shape(rc_full) == np.shape(gt))  # Check whether data was created properly
    return rc_full

def gp_from_gts(gts, cty=0.99):
    """Create GP [l,k,3] from
    genotypes [l,k,2], with prob.
    of genotype set to cty"""
    gs = np.sum(gts, axis=2)
    l, k = np.shape(gs)
    gp = 0.5 * (1 - cty * np.ones((l,k,3), dtype="f")) # 0.5 because 2 alternative allels
    gp[gs==0,0]=cty
    gp[gs==1,1]=cty
    gp[gs==2,2]=cty
    #(np.sum(gp, axis=2)==1).all()
    return gp

def gp_from_empirical(gts, bps, simulated_error, cty=0.99, verbose=False, oneModern=False):
    """Create GP [l,k,3] from
    genotypes [l,k,2], with prob.
    of genotype set to empirical error"""
    gs = np.sum(gts, axis=2)
    l, k = np.shape(gs)
    assert(l == len(bps))
    gp = np.zeros((l, k, 3))
    count_err_impute = 0
    count_err_hard = 0
    step = 1 if not oneModern else 2
    for i, bp in enumerate(bps):
        for j in np.arange(k, step=step):
            dosage_groundtruth = gs[i,j]
            # no empirical data available, use the simple model
            if len(simulated_error[bp][dosage_groundtruth]) == 0:
                gt1, gt2 = gts[i,j]
                gt1_hat = (gt1 + int(np.random.rand()<0.01))%2
                gt2_hat = (gt2 + int(np.random.rand()<0.01))%2
                if gt1 != gt1_hat or gt2 != gt2_hat:
                    count_err_hard += 1
                gts[i,j] = [gt1_hat, gt2_hat]
                gp[i,j,:] = 0.5*(1 - cty)
                gp[i,j,gt1_hat + gt2_hat] = cty
            else:
                gp_err = random.choice(simulated_error[bp][dosage_groundtruth])
                gp[i,j,:] = gp_err
                gt_err = np.argmax(gp_err)
                if gt_err != dosage_groundtruth:
                    count_err_impute += 1
                    if verbose:
                        print(f'GP is switched from {dosage_groundtruth} to {gt_err} according to empriical error at {bp}')
                    if gt_err == 0 or gt_err == 2:
                        gts[i,j] = [int(gt_err/2), int(gt_err/2)]
                    else:
                        gts[i,j] = random.choice([[0,1], [1,0]])
    if step == 2:
        for i, bp in enumerate(bps):
            for j in np.arange(1, k, step=2):
                gt1, gt2 = gts[i,j]
                gt1_hat = (gt1 + int(np.random.rand()<0.001))%2
                gt2_hat = (gt2 + int(np.random.rand()<0.001))%2
                gts[i,j] = [gt1_hat, gt2_hat]
                gp[i,j,:] = (1-0.999)/2
                gp[i,j, gt1_hat + gt2_hat] = 0.999
    

    print(f'{count_err_impute/(k/step)} error added per sample to simulate the effect of low coverage imputation...', flush=True)
    print(f'{count_err_hard/(k/step)} error added per sample for sites not in the imputation dictionary...', flush=True)
    return gts, gp



def shuffle_haplos(gts, rec, scale_cm=0.4, oneModern=False):
    """Shuffle genotypes with random waiting times 
    (mean scale in cM, exponentially distributed)"""
    m = np.max(rec)
    l = np.shape(gts)[0] # nr of loci
    n = np.shape(gts)[1] # Nr individuals
    scale_m = scale_cm/100.0
    k = int(m/scale_m * 5) # Five times overshooting to be sure...
    
    w = np.random.exponential(scale=scale_m, size=(n,k))
    w = np.cumsum(w, axis=1) # Get swapping points
    idx = np.searchsorted(rec, w) # Get swapping indices 
    
    gts_new = np.copy(gts)
    # Iterate over all individuals and all swap points
    step = 1 if not oneModern else 2
    for i in np.arange(n, step=step):
        for j in range(0,k-1,2):
            a,b= idx[i,j], idx[i,j+1]
            if a>=l:
                break
            gts_new[a:b,i,1] = gts[a:b,i,0]
            gts_new[a:b,i,0] = gts[a:b,i,1]
    assert((np.sum(gts, axis=2) == np.sum(gts_new, axis=2)).all()) # Sanity check
    return gts_new