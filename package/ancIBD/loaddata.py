"""
Class containing functions to load data needed for HMM IBD Run.
Returns all relevant data for HMM in standardized format:
1) Haplotype Probabilities (for ancestral allele)
2) Allele Frequencies 
3) Genetic Map
@ Author: Harald Ringbauer, September 2020
"""

import pickle as pickle
import numpy as np
import os as os
import h5py as h5py
import warnings

###############################
###############################

class LoadData(object):
    """Class to load data in uniform format"""
    path = "" # Path to load from (usually folder)
    output = True # Whether to print output
    
    def __init__(self, path=""):
        """Can remember path if given"""
        self.path=path
    
    def return_p(self, **kwargs):
        """Return array of Allele Frequencies [l]"""
        raise NotImplementedError("Implement This in specific subclass.")

    def return_map(self, **kwargs):
        """Return genetic map [l] in Morgan"""
        raise NotImplementedError("Implement This in specific subclass.")
        
    def return_haplotypes_ll(self, **kwargs):
        """Return haplotype likelihoods [4,l,2]"""
        raise NotImplementedError("Implement This in specific subclass.")
        
    def load_all_data(self, **kwargs):
        """Load all haplotype likelihoods
        haplotype likelihoods [2*n,l,2]
        derived allele frequencies [l]
        map in Morgan [l]"""
        htsl = self.return_haplotypes_ll()
        p = self.return_p()
        m = self.return_map()
        self.check_valid_data(htsl, p, m)
        return htsl, p, m
    
    def check_valid_data(self, htsl, p, m):
        """Check whether data in valid format"""
        assert(np.shape(htsl)[1]==len(p))
        assert(np.sum(np.isnan(htsl))==0) # Check that no NAN exists
        assert(len(p)==len(m))
        pass 
        
    def set_params(self, **kwargs):
        """Set the Parameters.
        Takes keyworded arguments"""
        for key, value in kwargs.items():
            setattr(self, key, value)
            
###############################
### Class for Simulated date

class LoadSimulated(LoadData):
    """Class to load simulated data saved in standard format."""
    r_gap = 1.0
    r_path = ""
    
    def load_all_data(self, **kwargs):
        ### Calculate Paths
        haplo_path = os.path.join(self.path, "haplo_ll.tsv")
        p_path = os.path.join(self.path, "p.tsv")
        
        ### Load the Data
        #htsl = pickle.load(open(haplo_path, "rb")) 
        htsl = np.loadtxt(haplo_path, delimiter="\t", dtype="float")
        p = np.loadtxt(p_path, delimiter="\t", dtype="float")
        if len(self.r_path)==0:
            #m = np.ones(len(p), dtype="float") * self.r_gap
            m = np.arange(len(p), dtype="float") * self.r_gap # Load absolute Positions
        
        self.check_valid_data(htsl, p, m)
        return htsl, p, m
    
###############################
### Class for HDF5 data

class LoadHDF5(LoadData):
    """Class to load HDF data for 2 individuals in standard format."""
    path_h5=""
    iids = []
    ch = 3   # Which chromosome to load
    min_error=1e-5 # Minimum Probability of genotyping error
    p_col = "variants/AF_ALL" # The hdf5 column with der. all freqs
    
    def return_map(self, f):
        """Return the recombination map"""
        m = f["variants/MAP"][:]
        return m
    
    def return_p(self, f):
        """Return array of Allele Frequencies [l]
        self.p_col: The HDF5 field with the allele frequencies
        If not given, use default p=0.5
        """
        if len(self.p_col)==0:
            p = 0.5 * np.ones(len(f["variants/MAP"]))
        else:
            p = f[self.p_col][:] # Load the Allele Freqs from HDF5
            if len(p.shape) > 1:
                p = p[:,0]
        return p
        
    def get_individual_idx(self, f, iid="", f_col="samples"):
        """Return index of individual iid"""
        samples = f[f_col].asstr()[:]
        idx = (samples == iid)
        
        if np.sum(idx)==0:  # Sanity Check wheter IID found at all (raise stop)
            raise RuntimeWarning(f"No entry in H5 found for iid: {iid}") 
        if np.sum(idx)>1:   # Sanity Check wether multiple IIDs found (warning only)
            warnings.warn(f"{np.sum(idx)} entries found for iid: {iid}", RuntimeWarning)    
        idx=np.where(idx)[0][0] # Get the first IID match index
        return idx  
    
    def get_haplo_prob(self, f, idx):
        """Get haploid ancestral probability for indivual [2,l]"""
        h1 = f["calldata/GT"][:,idx,:].T
        m = np.max(f["calldata/GP"][:,idx,:], axis=1)
        m = np.minimum(m,  1 - self.min_error)
        h1 = (1 - h1) * m + h1 * (1 - m) # Probability of being ancestral
        return h1
        
    def load_all_data(self, **kwargs):
        """ Return haplotype likelihoods [4,l] for anc. allele
        derived allele frequencies [l]
        map in Morgan [l]"""
        path_h5_ch = f"{self.path}{self.ch}.h5"
        with h5py.File(path_h5_ch, "r") as f:
            m = self.return_map(f)
            p = self.return_p(f)
            
            idx1 = self.get_individual_idx(f, self.iids[0])
            idx2 = self.get_individual_idx(f, self.iids[1])
            h1 = self.get_haplo_prob(f, idx1)
            h2 = self.get_haplo_prob(f, idx2)
            htsl = np.concatenate((h1,h2), axis=0)
        
        self.check_valid_data(htsl, p, m)
        return htsl, p, m, self.iids
    
class LoadHDF5Multi(LoadHDF5):
    """Class to load HDF5 data for multiple individuals.
    Default now also for only pairs of individuals."""
    path_h5=""
    iids = []
    ch = 3   # Which chromosome to load
    min_error = 1e-5 # Minimum Probability of genotyping error  
    p_col = "" # The hdf5 column with der. all freqs. If given, use the field
    # If "default" use p=0.5 everywhere
    # default value for p_col in my hdf5s: variants/AF_ALL
    ploidy = 2 # ploidy of the individuals
    
    def get_haplo_prob(self, f, idcs):
        """Get haploid ancestral probability for n individuals 
        Return [n,l,2] array"""
        h1 = f["calldata/GT"][:,idcs,:] ### Get l, n, 2
        l,n,_ = np.shape(h1) # Get the nr loci
        m = np.max(f["calldata/GP"][:,idcs,:], axis=2)  
        m = np.minimum(m,  1 - self.min_error) # Max Cap of certainty
        h1 = (1-h1) * m[:,:,None] + h1 * (1 - m[:,:,None]) # Probability of being ancestral
        h1 = np.swapaxes(h1, 0, 1) # ->n,l,2
        h1 = np.swapaxes(h1, 1, 2) #-> n,2,l
        h1 = h1.reshape((2*n,l))
        #h1 = np.swapaxes(h1, 0, 1)
        return h1
    
    def filter_valid_data(self, hts, p, m, bp):
        """Filter to SNPs with fully valid data. 
        Return filtered data."""
        idx = ~(np.isnan(hts).any(axis=0)) # Flag all markers that are not null
        
        ### Output and Raise warning if too much data missing:
        if self.output:
            print(f"Filtering to {np.sum(idx)}/{len(idx)} SNPs with GP data (on target iids)")   
        if np.mean(idx)<1:  # Missing Data detected
            print(f"Attention: Some data in GP field is missing. Ideally, all GP entries are set.")   
        if np.mean(idx)<0.8:  # Less than 80 percent of data there
            print(f"Too much data missing: {np.sum(idx)}/{len(idx)} SNPs have GP ({np.mean(idx)*100}% there)")
            
        return hts[:,idx], p[idx], m[idx], bp[idx]
    
    def load_all_data(self, **kwargs):
        """ Return haplotype likelihoods [n*2,l] for anc. allele.
        along first axis: 2*i, 2*(i+1) haplotype of ind i
        derived allele frequencies [l]
        map in Morgan [l]
        bp positions [l]
        """
        path_h5_ch = f"{self.path}{self.ch}.h5"
        with h5py.File(path_h5_ch, "r") as f:
            m = self.return_map(f)            
            idcs = np.array([self.get_individual_idx(f, iid) for iid in self.iids])
            sort = np.argsort(idcs)   # Get the sorting Indices [has to go low to high]
            samples = self.iids[sort] # Get them in sorted order
            hts = self.get_haplo_prob(f, idcs[sort]) # Still does old loading
            bp = f['variants/POS'][:]
            if len(self.p_col)>0:
                p = self.get_p_hdf5(f, self.p_col)
            else:
                p = self.get_p(hts)  # Calculate Mean allele frequency from sample subset
                
        ### Filter to Valid data
        hts, p, m, bp = self.filter_valid_data(hts, p, m, bp)
        
        self.check_valid_data(hts, p, m)
        return hts, p, m, bp, samples
    
    def get_p(self, htsl):
        """Get Allele frequency from haplotype probabilities.
        Return array of derived allele freqs [l]"""
        p_anc = np.mean(htsl, axis=0) # Take the expected AF of ancestral.
        p_der = 1 - p_anc             # Get the derived allele frequency
        return p_der
    
    def get_p_hdf5(self, f, col):
        """Get allele frequs from HDF f 
        in dataset col.
        Return array of derived allele freqs [l]"""
        if col=="default":
            p = 0.5 * np.ones(len(f["variants/MAP"]))
        else:
            p = f[col][:] # Load the Allele Freqs from HDF5
            
        ### Filter if p is matrix (e.g. sometimes for three alt alleles)
        if len(np.shape(p))==2:
            p=p[:,0]
        elif len(np.shape(p))!=1:
            raise RuntimeWarning(f"Allele Frequency field in h5 {col} has invalid dimensions!")
            
        return p

class LoadH5Multi2(LoadHDF5Multi):
    """Update to more accurate genotype probabilities
    from diploid to haploid."""
    path_h5=""
    iids = []
    ch = 3   # Which chromosome to load
    min_error = 1e-3  # Minimum Probability of genotyping error  
    pph_error = 1e-2 #1e-3 # Point Phase Error
    p_col = "" # The hdf5 column with der. all freqs. If given, use the field
    # If "default" use p=0.5 everywhere
    # default value for p_col in my hdf5s: variants/AF_ALL
    
    def get_haplo_prob(self, f, idcs, ploidy=2):
        """Get haploid ancestral probability for n individuals 
        Return [n,l,2] array. Calculated from GP and GT
        in the proper way.
        ploidy: it can be either an integer or an array of integers
        when it's an integer, it must be either 1 or 2 and then we assume that all individuals have the same ploidy
        when it's an array, it must have the same length as idcs and then it specifies the ploidy of each individual
        """
            
        gt = f["calldata/GT"][:,idcs,:] ### Get l, n, 2 array of gt
        l,n,_ = np.shape(gt) # Get the nr loci
        gp = f["calldata/GP"][:,idcs,:] ### Get l, n, 3 array of gp
        
        if isinstance(ploidy, int) and ploidy == 1:
            # for haploid samples, the entry in h1[:,:,1] will never be used, so it's fine to leave it as all zeros
            h1 = np.zeros((l,n,2), dtype=np.float64)
            h1[:,:,0] = gp[:,:,0]
        else:
            # ploidy is either an integer or an array of integers
            if isinstance(ploidy, int):
                assert ploidy == 2
                ploidy = 2*np.ones(n, dtype=np.int64)
            ### The homozygotes
            g00 = gp[:,:,0]  # homo 00 prob.
            ### Heterozygote prob
            gp1 = gp[:,:,1]
            # set the het prob for haploid samples to be 0
            # for haploid samples the gp[:,:,1] actually encodes g11
            gp1[:, ploidy == 1] = 0.0 
            ## g11 not needed
        
            ### The max GT probabilities
            g01 = gp1[:,:] * (gt[:,:,0]==0) * (gt[:,:,1]==1) # het 01 prob.
            g10 = gp1[:,:] * (gt[:,:,0]==1) * (gt[:,:,1]==0) # het 10 prob.
        
            ### The default probability from the non-max GT ones
            # Assume it is split up 50/50 between 01 10 states
            idx = (gt[:,:,0]==0) * (gt[:,:,1]==0)
            g01[idx] = gp1[idx]/2
            g10[idx] = gp1[idx]/2
        
            idx = (gt[:,:,0]==1) * (gt[:,:,1]==1)
            g01[idx] = gp1[idx]/2
            g10[idx] = gp1[idx]/2
            ### The added probability from the non-max GT ones
        
            # Add point phasing error  here (swap g01 and g10 with prob. x)
            #g01t = g01 * (1 - self.pph_error) + g10 * self.pph_error
            #g10t = g10 * (1 - self.pph_error) + g01 * self.pph_error
        
            h1 = np.zeros((l,n,2), dtype=np.float64)
            h1[:,:,0] = g00 + g01
            h1[:,:,1] = g00 + g10
        
        h1 = np.swapaxes(h1, 0, 1) #  l,n,2->n,l,2
        h1 = np.swapaxes(h1, 1, 2) # -> n,2,l
        h1 = h1.reshape((2*n,l))
        
        ### Add minimum error to avoid overly extreme confidence genotypes
        h1 = np.maximum(h1, self.min_error)
        h1 = np.minimum(h1, 1-self.min_error)
        
        if self.output:
            print(f"Thresholding GP at {1 - self.min_error}")
            #print(f"Phase. Error added: {self.pph_error}")
        
        return h1


class LoadDict(LoadHDF5Multi):
    """Update to more accurate genotype probabilities
    from diploid to haploid."""
    iids = []
    min_error = 1e-3  # Minimum Probability of genotyping error  
    pph_error = 1e-2 #1e-3 # Point Phase Error
    p_col = "" # The hdf5 column with der. all freqs. If given, use the field
    # If "default" use p=0.5 everywhere
    # default value for p_col in my hdf5s: variants/AF_ALL
    
    def load_all_data(self, **kwargs):
        """ Return haplotype likelihoods [n*2,l] for anc. allele.
        along first axis: 2*i, 2*(i+1) haplotype of ind i
        derived allele frequencies [l]
        map in Morgan [l]
        bp positions [l]
        """
        f = self.path # dont provide a path, provide a dict 
        
        m = self.return_map(f)            
        idcs = np.array([self.get_individual_idx(f, iid) for iid in self.iids])
        sort = np.argsort(idcs)   # Get the sorting Indices [has to go low to high]
        samples = self.iids[sort] # Get them in sorted order
        hts = self.get_haplo_prob(f, idcs[sort]) # Still does old loading
        bp = f['variants/POS'][:]
        if len(self.p_col)>0:
            p = self.get_p_hdf5(f, self.p_col)
        else:
            p = self.get_p(hts)  # Calculate Mean allele frequency from sample subset                
        
        ### Filter to Valid data
        hts, p, m, bp = self.filter_valid_data(hts, p, m, bp)
        
        self.check_valid_data(hts, p, m)
        return hts, p, m, bp, samples


    def get_individual_idx(self, f, iid="", f_col="samples"):
        """Return index of individual iid"""
        samples = f[f_col][:]
        idx = (samples == iid)
        
        if np.sum(idx)==0:  # Sanity Check wheter IID found at all (raise stop)
            raise RuntimeWarning(f"No entry in H5 found for iid: {iid}") 
        if np.sum(idx)>1:   # Sanity Check wether multiple IIDs found (warning only)
            warnings.warn(f"{np.sum(idx)} entries found for iid: {iid}", RuntimeWarning)    
        idx=np.where(idx)[0][0] # Get the first IID match index
        return idx  


###############################    
###############################
### Factory Method
    
def load_loaddata(l_model="simulated", path="", **kwargs):
    """Factory method to return the right loading Model"""
    if l_model == "simulated":
        l_obj = LoadSimulated(path=path, **kwargs)
    elif l_model == "hdf5double":  # use overall allele frequency
        l_obj = LoadHDF5(path=path, **kwargs)
    elif l_model == "hdf5": # Original model, for legacy reasons
        l_obj = LoadHDF5Multi(path=path, **kwargs)
    elif l_model == "h5":  # Latest model, as described in paper.
        l_obj = LoadH5Multi2(path=path, **kwargs)
    elif l_model == "dict": # in memory python dict
        l_obj = LoadDict(path=path, **kwargs)
    else:
        raise NotImplementedError("Loading Model not found!")
    return l_obj