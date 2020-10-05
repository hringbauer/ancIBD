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
        """Key Method. 
        haplotype likelihoods [4,l,2]
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
    
    def rmap_to_gaps(self, r_map=[], cm=False, min_gap=1e-10, max_gap=0.05):
        """Return the recombination map gaps [in Morgan]
        Input: Map Positions [l] units see cm below
        Return: Rec. Distance Array [l]
        cm: Whether input is in centimorgan or morgan
        min_cap: Minimum Map Gap between Loci to cap
        max_cap: Maximum Mapg Gap between Loci to cap"""
        gaps = np.zeros(len(r_map))  # Make new array
        gaps[1:] = r_map[1:] - r_map[:-1]  # Calculate Differences
        assert(np.min(gaps) >= 0)
        if cm == True:
            gaps = gaps / 100     # Normalize to Morgan if map in cM
         
        ### Extend the minimum gap where needed
        gaps = np.maximum(gaps, min_gap)
        gaps = np.minimum(gaps, max_gap)

        if self.output == True:
            max_g = np.max(gaps)
            print(f"Minimum Genetic Map: {np.min(r_map):.4f} Morgan")
            print(f"Maximum Genetic Map: {np.max(r_map):.4f} Morgan")
            print(f"Gaps bigger than 0.1 cM: {np.sum(gaps > 0.001)}")
            print(f"Maximum Gap: {max_g * 100:.4f} cM")
            print(f"Upper Gap Cutoff: {max_gap * 100:.4f} cM")
        return gaps
        
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
            m = np.ones(len(p), dtype="float") * self.r_gap
        
        self.check_valid_data(htsl, p, m)
        return htsl, p, m
    
###############################
### Class for HDF5 data

class LoadHDF5(LoadData):
    """Class to load simulated data saved in standard format."""
    path_h5=""
    iids = []
    ch = 3   # Which chromosome to load
    min_gap=1e-10 # Minimum Map Gap between two loci
    max_gap=0.05  # Maximum Map Gap between two loci
    min_error=1e-5 # Minimum Probability of genotyping erorr
    
    def return_map(self, f):
        """Return the recombination map"""
        r_map = f["variants/MAP"][:]
        m = self.rmap_to_gaps(r_map=r_map, cm=False)
        return m
    
    def return_p(self, f):
        """Return array of Allele Frequencies [l]"""
        p = 0.5 * np.ones(len(f["variants/MAP"]))
        return p
        
    def get_individual_idx(self, f, iid="", f_col="samples"):
        """Return index of individual iid"""
        samples = f[f_col][:]
        idx = (samples == iid)
        assert(np.sum(idx)==1) # Sanity Check
        idx=np.where(idx)[0][0]
        return idx  
    
    def get_haplo_prob(self, f, idx):
        """Get  haploid ancestral probability for indivual"""
        h1 = f["calldata/GT"][:,idx,:].T
        m = np.max(f["calldata/GP"][:,idx,:], axis=1)
        m = np.minimum(m,  1 - self.min_error)
        h1 = (1-h1) * m + h1 * (1 - m) # Probability of being ancestral
        return h1
        
    def load_all_data(self, **kwargs):
        """Key Method. 
        haplotype likelihoods [4,l,2]
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
        return htsl, p, m
    
###############################    
###############################
### Factory Method
    
def load_loaddata(l_model="simulated", path="", **kwargs):
    """Factory method to return the right loading Model"""
    if l_model == "simulated":
        l_obj = LoadSimulated(path=path, **kwargs)
    elif l_model == "hdf5":
        l_obj = LoadHDF5(path=path, **kwargs)
    else:
        raise NotImplementedError("Loading Model not found!")
    return l_obj