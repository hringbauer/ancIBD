"""
Class containing functions to calculate Emission Probabilities.
Contains Sub-Classes for various types of data, 
as well as a factory Method to return the 
correct subclass
@ Author: Harald Ringbauer, September 2020
"""

import numpy as np
import itertools as it

###############################
###############################

class Emissions(object):
    """Class for emission probabilities
    Has methods to return emission probabilities"""

    def give_emission_matrix(self, **kwargs):
        """Return Emission Matrix - for every possible set of states"""
        raise NotImplementedError("Implement This in specific subclass.")
    
    def give_emission_log(self, **kwargs):
        """Return the full emission Probability directly in Log Space. 
        By default just log the emissiom matrix. Can be overwritten for
        nicer computational properties"""
        m = self.give_emission_matrix(**kwargs)
        #assert(np.min(m)>0)
        return np.log(m)

    def set_params(self, **kwargs):
        """Set the Parameters.
        Takes keyworded arguments"""
        for key, value in kwargs.items():
            setattr(self, key, value)
            
class HaplotypeSharingEmissions(Emissions):
    """Class for emission probabilities of diploid haplotypes who
    possibly share one haplotype"""
    
    def hw_prob_haplo_pp(self, ht_p, p):
        """Calculate HW Probabilitiy of haplotype gt
        ht_p: [l] Array of haplotype prob for ancestral, l locis,
        p: [l] Array of (derived) allele frequencies
        returns [l] vector of HW prob of gt"""
        prob = (1 - ht_p[:]) * p +  ht_p[:] * (1-p)
        return prob

    def hw_prob_haplo_share_pp(self, ht_p1, ht_p2, p):
        """Calculate probability of sharing haplotypes
        ht_p1, ht_p2: Array of haplotype likelihood for ancestral, l locis
        p: [l] Array of (derived) allele frequencies
        returns [l] vector of prob that shared haplotype"""
        s1 = ht_p1 * ht_p2 * (1-p) # Prob both ancestral
        s2 = (1-ht_p1) * (1 -ht_p2) * p # Prob both derived
        prob = s1 + s2
        return prob

    def hw_prob_haplos_pp(self, hts_p, p):
        """Calculate HW Probabilitiy of haplotype gt
        gt: [k,l] Array of haplotype prob for ancestral
        p: [l] Array of (derived) allele frequencies
        returns [l] vector of HW prob of gt"""
        prob = (1-hts_p[:,:]) * p +  hts_p[:,:] * (1-p)
        prob_tot = np.prod(prob, axis=0)
        return prob_tot

    def hw_probs_shared(self, hts_p, p, shared=(0,2), dtype="float"):
        """Give emission probabilities for shared state.
        Assume 0/1 2/3 are diploid haplotypes, and 0/2 are shared
        hts_p: [4,l] Array of four ancestral haplotype probabilities
        p: [l] array of derived genotype probability.
        shared: tuple of length 2 giving the indices of the shared haplotypes"""
        not_shared = [i for i in range(0,4) if i not in shared]
        assert(len(not_shared)==2 & len(shared)==2)
        p_hw1 = self.hw_prob_haplo_pp(hts_p[not_shared[0],:], p=p)
        p_hw2 = self.hw_prob_haplo_pp(hts_p[not_shared[1],:], p=p)
        p_shared = self.hw_prob_haplo_share_pp(hts_p[shared[0],:],hts_p[shared[1],:], p=p)
        p = p_hw1*p_hw2*p_shared
        return p

    def give_emission_matrix(self, hts_p, p, dtype="float"):
        """Give emission Matrix for 5-state HMM.
        0th state: HW 1st-4th State: Haplotype Copying
        Input: p: [l] Array of (derived) allele frequencies
        hts_p: [4,l] Array of four ancestral haplotype probabilities.
        Return: emission matrix [5,l]."""
        l = np.shape(hts_p)[1]
        e_mat = np.zeros((5,l), dtype=dtype)
        e_mat[0,:] = self.hw_prob_haplos_pp(hts_p,p=p)
        for i, (j,k) in enumerate(it.product([0,1],repeat=2)):
            e_mat[i+1,:] = self.hw_probs_shared(hts_p, p=p, shared=[j,k+2])
        return e_mat
    
def load_emission_model(e_model="haploid_gl"):
    """Factory method to return the right Emission Model"""
    if e_model == "haploid_gl":
        e_obj = HaplotypeSharingEmissions()
    else:
        raise NotImplementedError("Emission Model not found!")
    return e_obj