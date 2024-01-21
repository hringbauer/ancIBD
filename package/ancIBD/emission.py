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
    
class HaplotypeSharingEmissions2(Emissions):
    """Class for emission probabilities of diploid haplotypes who
    possibly share one haplotype. Update to new emissions based
    on that output of Glimpse are posterior probabilities"""
    
    
    def hw_prob_haplo_share_pp(self, ht_p1, ht_p2, p, min_p = 0.001):
        """Calculate probability of sharing haplotypes
        ht_p1, ht_p2: Array of haplotype likelihood for ancestral, l locis
        p: [l] Array of (derived) allele frequencies
        returns [l] vector of prob that shared haplotype"""
        ### Cut off to minimal allele frequencies
        pt = np.maximum(p, min_p)
        pt = np.minimum(pt, 1-min_p)
        
        s1 = ht_p1 * ht_p2 / (1-pt) # Prob both ancestral
        s2 = (1-ht_p1) * (1-ht_p2) / pt # Prob both derived
        prob = s1 + s2
        return prob

    def hw_probs_shared(self, hts_p, p, shared=(0,2), dtype="float"):
        """Give emission probabilities for shared state.
        hts_p: [4,l] Array of four ancestral haplotype probabilities
        p: [l] array of derived genotype probability.
        shared: tuple of length 2 giving the indices of the shared haplotypes"""
        ### Sanity Check
        not_shared = [i for i in range(0,4) if i not in shared]
        assert(len(not_shared)==2 & len(shared)==2) # Sanity Check
        
        # Actual Calculation
        p_shared = self.hw_prob_haplo_share_pp(hts_p[shared[0],:],hts_p[shared[1],:], p=p)
        return p_shared

    def give_emission_matrix(self, hts_p, p, dtype="float"):
        """Give emission Matrix for 5-state HMM.
        0th state: HW 1st-4th State: Haplotype Copying
        Input: p: [l] Array of (derived) allele frequencies
        hts_p: [4,l] Array of four ancestral haplotype probabilities.
        Return: emission matrix [5,l]."""
        l = np.shape(hts_p)[1]
        e_mat = np.zeros((5,l), dtype=dtype)
        e_mat[0,:] = 1 # Because simply proportional to P(D)
        for i, (j,k) in enumerate(it.product([0,1],repeat=2)):
            e_mat[i+1,:] = self.hw_probs_shared(hts_p, p=p, shared=[j,k+2])
        return e_mat

class HaplotypeSharingEmissions3(HaplotypeSharingEmissions2):
    """
    Basically the same as HaplotypeSharingEmissions2, but with two additional states to account for IBD2 region
    """


    def give_emission_matrix(self, hts_p, p, p_min=1e-3, dtype='float'):
        """
        Give emission matrix for 7-state HMM.
        0th state: non-IBD state
        1st-4th state: IBD1 state
        5-7th state: IBD2 state
        """
        # to avoid division by zero error, set p to be at least p_min and at most 1 - p_min
        p = np.maximum(p, p_min)
        p = np.minimum(p, 1-p_min)
        
        l = np.shape(hts_p)[1]
        e_mat = np.zeros((7,l), dtype=dtype)
        e_mat[:5,:] = super().give_emission_matrix(hts_p, p, dtype)

        x1A = 1 - hts_p[0,:]
        x1B = 1 - hts_p[1,:]
        x2A = 1 - hts_p[2,:]
        x2B = 1 - hts_p[3,:]
        # IBD2 situation 1: 1A=2A, 1B=2B
        e_mat[5,:] = (1-x1A)*(1-x1B)*(1-x2A)*(1-x2B)/((1-p)**2) + \
                (1-x1A)*x1B*(1-x2A)*x2B/(p*(1-p)) + \
                x1A*(1-x1B)*x2A*(1-x2B)/(p*(1-p)) + \
                x1A*x1B*x2A*x2B/(p**2)
        
        # IBD2 situation 2: 1A=2B, 1B=2A
        e_mat[6,:] = (1-x1A)*(1-x1B)*(1-x2A)*(1-x2B)/((1-p)**2) + \
                (1-x1A)*x1B*x2A*(1-x2B)/(p*(1-p)) + \
                x1A*(1-x1B)*(1-x2A)*x2B/(p*(1-p)) + \
                x1A*x1B*x2A*x2B/(p**2)

        return e_mat

############################################ Experimental with X chromosomes ############################################
class HaplotypeSharing_twoHaploid(HaplotypeSharingEmissions2):
    """Class for emission probabilities of two haploid genomes (e.g, two male X chromosome) 
    who possibly share one haplotype"""
    
    def give_emission_matrix(self, hts_p, p, dtype="float"):
        """Give emission Matrix for 2-state HMM.
        0th state: HW 1st State: Haplotype Copying
        Input: p: [l] Array of (derived) allele frequencies
        hts_p: [4,l] Array of two ancestral haplotype probabilities. 
            Note that for haploid samples, hts_p[1,:] and hts_p[4,:] are not used and the values there are meaningless.
            I keep these two rows only for consistency with the diploid case.
        Return: emission matrix [2,l]."""
        l = np.shape(hts_p)[1]
        e_mat = np.zeros((2,l), dtype=dtype)
        e_mat[0,:] = 1 # Because simply proportional to P(D)
        e_mat[1,:] = self.hw_probs_shared(hts_p, p=p, shared=[0,2])
        return e_mat
    
class HaplotypeSharing_hapVSdiploid(HaplotypeSharingEmissions2):
    """Class for emission probabilities of one haploid and one diploid genomes (e.g, between a male and a female chromosome)
    who possibly share one haplotype"""

    def give_emission_matrix(self, hts_p, p, dtype="float"):
        """Give emission Matrix for 3-state HMM. 
        Note that the code assumes that the first sample is haploid and the second sample is diploid.
        0th state: HW 
        1st-2nd State: Haplotype Copying
        Input: p: [l] Array of (derived) allele frequencies
        hts_p: [4,l] Array of two ancestral haplotype probabilities. 
            Note that for hts_p[1,:] is not used and the values there are meaningless.
            I keep these it only for consistency with the diploid case.
        Return: emission matrix [3,l]."""

        l = np.shape(hts_p)[1]
        e_mat = np.zeros((3,l), dtype=dtype)
        e_mat[0,:] = 1 # Because simply proportional to P(D)
        e_mat[1,:] = self.hw_probs_shared(hts_p, p=p, shared=[0,2])
        e_mat[2,:] = self.hw_probs_shared(hts_p, p=p, shared=[0,3])
        return e_mat


def load_emission_model(e_model="haploid_gl"):
    """Factory method to return the right Emission Model"""
    if e_model == "haploid_gl":
        e_obj = HaplotypeSharingEmissions()
    elif e_model == "haploid_gl2":
        e_obj = HaplotypeSharingEmissions2()
    elif e_model == 'IBD2':
        e_obj = HaplotypeSharingEmissions3()
    elif e_model == "twoHaploid":
        e_obj = HaplotypeSharing_twoHaploid()
    elif e_model == "hapVSdiploid":
        e_obj = HaplotypeSharing_hapVSdiploid()
    else:
        raise NotImplementedError("Emission Model not found!")
    return e_obj