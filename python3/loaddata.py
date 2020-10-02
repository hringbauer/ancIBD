"""
Class containing functions to load data needed for HMM Run.
Returns all relevant data for HMM in unified format:
1) Haplotype Probabilities 
2) Allele Frequencies 
3) Genetic Map
@ Author: Harald Ringbauer, September 2020
"""

import pickle as pickle
import numpy as np
import os as os

###############################
###############################

class LoadData(object):
    """Class to load data in uniform format"""
    path = "" # Path to load from (usually folder)
    
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
        htsl = pickle.load(open(haplo_path, "rb")) 
        p = np.loadtxt(p_path, delimiter="\t", dtype="float")
        if len(self.r_path)==0:
            m = np.ones(len(p), dtype="float") * self.r_gap
        
        self.check_valid_data(htsl, p, m)
        return htsl, p, m
    
###############################    
###############################
### Factory Method
    
def load_loaddata(l_model="simulated", path="", **kwargs):
    """Factory method to return the right loading Model"""
    if l_model == "simulated":
        l_obj = LoadSimulated(path=path, **kwargs)
    else:
        raise NotImplementedError("Loading Model not found!")
    return l_obj