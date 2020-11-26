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
        return p
        
    def get_individual_idx(self, f, iid="", f_col="samples"):
        """Return index of individual iid"""
        samples = f[f_col].asstr()[:]
        idx = (samples == iid)
        if np.sum(idx)!=1:
            raise RuntimeWarning(f"{np.sum(idx)} entries found for {iid}")
        assert(np.sum(idx)>0) # Sanity Check
        idx=np.where(idx)[0][0]
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
    Developer Note: In progress of making this default even for pairs of individuals."""
    path_h5=""
    iids = []
    ch = 3   # Which chromosome to load
    min_error = 1e-5 # Minimum Probability of genotyping error  
    p_col = "" # The hdf5 column with der. all freqs. If given, use the field
    # If "default" use p=0.5 everywhere
    # default value for p_col in my hdf5s: variants/AF_ALL
    
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
    
    def load_all_data(self, **kwargs):
        """ Return haplotype likelihoods [n*2,l] for anc. allele.
        along first axis: 2*i, 2*(i+1) haplotype of ind i
        derived allele frequencies [l]
        map in Morgan [l]"""
        path_h5_ch = f"{self.path}{self.ch}.h5"
        with h5py.File(path_h5_ch, "r") as f:
            m = self.return_map(f)            
            idcs = np.array([self.get_individual_idx(f, iid) for iid in self.iids])
            sort = np.argsort(idcs)   # Get the sorting Indices [has to go low to high]
            samples = self.iids[sort] # Get them in sorted order
            hts = self.get_haplo_prob(f, idcs[sort])
        
            if len(self.p_col)>0:
                p = self.get_p_hdf5(f, self.p_col)  
            else:
                p = self.get_p(hts)  # Calculate Mean allele frequency from subset
            
        self.check_valid_data(hts, p, m)
        return hts, p, m, samples
    
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
        return p
    
###############################    
###############################
### Factory Method
    
def load_loaddata(l_model="simulated", path="", **kwargs):
    """Factory method to return the right loading Model"""
    if l_model == "simulated":
        l_obj = LoadSimulated(path=path, **kwargs)
    elif l_model == "hdf5double":  # use overall allele frequency
        l_obj = LoadHDF5(path=path, **kwargs)
    elif l_model == "hdf5":
        l_obj = LoadHDF5Multi(path=path, **kwargs)
    else:
        raise NotImplementedError("Loading Model not found!")
    return l_obj