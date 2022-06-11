"""
Class for calculating Transition Probabilities, i.e. 
infitesimal transition Matrices.
Contains Sub-Classes, as well as factory Method.
@ Author: Harald Ringbauer, 2021, All rights reserved
"""
import numpy as np


class Transitions(object):
    """Class for transition probabilities.
    Has methods to return them"""
    trans_mat = []    # The full transition Matrix
    r_map = []        # Legacy Variable
    output=True

    def __init__(self):
        """Initialize Class"""
        #self.calc_transitions()
        pass

    def calc_transition_rate(self, n=0):
        """Calulate and return Transition Matrix"""
        raise NotImplementedError("Implement This!")
        
    def give_transition_rate(self):
        """Give the transition_matrix"""
        return self.trans_mat

    def set_params(self, **kwargs):
        """Set the Values."""
        for key, value in kwargs.items():
            setattr(self, key, value)


class FiveStateTransitions(Transitions):
    """Implements the Model Transitions [units generally in Morgan]"""
    ibd_in = 0.0005     # The rate of jumping into IBD copying state
    ibd_out = 0.001     # The rate of jumping out of IBD state
    ibd_jump = 0.05     # The rate of jumping within IBD to other haplotype pair
    
    min_gap=1e-10 # Minimum Map Gap between two loci [Morgan]
    max_gap=0.05  # Maximum Map Gap between two loci [Morgan]

    def calc_transition_rate(self, submat33=True, n=4):
        """Return Transition Rate Matrix [k,k] to exponate.
        n: Number of symetric IBD states. Usually four (2x2 copying possibilities)
        submat33: Whether to only fill in the first 3 states """

        if submat33 == True:
            # Initialize Transition Matrix Only do 3 States (bc symmetry)
            t_mat = -np.ones((3, 3))
        else:
            t_mat = -np.ones((n + 1, n + 1))  # Initialize Transition Matrix

        t_mat[1:, 0] = self.ibd_out  # The rate of jumping out roh
        t_mat[0, 1:] = self.ibd_in / n  # Jumping into any ROH State
        t_mat[1:, 1:] = self.ibd_jump / n  # Jumping between ROH State

        # Do the Diagonal (do the usual model - for inf. substract 1)
        di = np.diag_indices(np.shape(t_mat)[0])
        t_mat[di] =  - self.ibd_out - self.ibd_jump + \
                       self.ibd_jump / (n)  # To account for jump to self
        t_mat[0, 0] = - self.ibd_in   # The rate of staying in diploid

        # Sanity Check if everything was filled correctly
        if (submat33 == False):
            assert(np.all(np.sum(t_mat, axis=1) > -0.0001))
            assert(np.all(np.sum(t_mat, axis=1) < 0.0001))
            #print(np.sum(t_mat, axis=1))  # Output

        self.trans_mat = t_mat
        return t_mat
    
    def full_transition_matrix(self, r_vec, t=[], n=4, submat33=True):
        """Compute and return the full transition Matrix.
        Calculates the first 3 states (not more needed by symmetry)
        t full Transition Matrix [k,k]. NO LOG STATE. If not given caluclate
        r_vec Map Length of Jumps [l] in Morgan
        n: Number of symmetric, non-background states"""
        ### infinitesimal rate
        if len(t)==0:
            t = self.calc_transition_rate(submat33=submat33, n=n)
            
        ### Full matrix
        r_vec = self.rmap_to_gaps(r_map=r_vec)
        if submat33:
            t = self.prep_3x3matrix(t, n)
        t_mat = self.exponentiate_r(t, r_vec)

        # Normalize to transition rate into non-collapsed state
        if submat33:
            t_mat[:, :2, 2] = t_mat[:, :2, 2] / (n - 1)
        return t_mat
    
    def prep_3x3matrix(self, t, n):
        """Prepare and return the grouped 3x3 Matrix 
        (grouped 3rd State: Everything in other IBD states)
        t: Origianl transition matrix (only first three entry important)
        n: Number of symmetric states
        """
        if self.output:
            print(f"HMM State Number: {n}")
        # Initiate to -1 (for later Sanity Check if everything is filled)
        t_simple = -np.ones((3, 3), dtype="float")
        t_simple[:2, :2] = t[:2, :2]
        # The probability of staying when in diff. ROH State:
        t_simple[2, 2] = -np.sum(t[2, :2])  # Minus the rates of Jumping
        # Jumping into 3rd state: Sum over all reference states
        t_simple[:2, 2] = t[:2, 2] * (n - 1) # Collapse the  other states
        t_simple[2, :2] = t[2, :2]  # The jumping out probability is the same
        return t_simple


    def exponentiate_r(self, rates, rec_v):
        """Calculates exponentiation of the rates matrix with rec_v
        rates: 2D Matrix of transitions
        rec_v: Array of length l"""
        eva, evec = np.linalg.eig(rates)  # Do the Eigenvalue Decomposition
        assert(np.max(eva) <= 1)   # Sanity Check whether rate Matrix
        evec_r = np.linalg.inv(evec)    # Do the Inversion
        # Create vector of the exponentiated diagonals
        d = np.exp(rec_v[:, None] * eva)
        # Use some Einstein Sum Convention Fun (C Speed):
        res = np.einsum('...ik, ...k, ...kj ->...ij', evec, d, evec_r)
        # Make sure that all transition rates are valuable
        assert(0 <= np.min(res))
        return res
    
    def rmap_to_gaps(self, r_map=[], cm=False):
        """Return the recombination map gaps [in Morgan]
        Input: Map Positions [l] (units see cm parameter below)
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
        gaps = np.maximum(gaps, self.min_gap)
        gaps = np.minimum(gaps, self.max_gap)

        if self.output == True:
            max_g = np.max(gaps)
            print(f"Minimum Genetic Map: {np.min(r_map):.4f} Morgan")
            print(f"Maximum Genetic Map: {np.max(r_map):.4f} Morgan")
            print(f"Gaps bigger than 0.1 cM: {np.sum(gaps > 0.001)}")
            print(f"Maximum Gap: {max_g * 100:.4f} cM")
            print(f"Upper Gap Cutoff: {self.max_gap * 100:.4f} cM")
        return gaps

############################################
# Factory method to return Transition Object

def load_transition_model(t_model="standard"):
    """Load the Transition Model"""
    if t_model == "standard":
        t_obj = FiveStateTransitions()
    else:
        raise NotImplementedError("Transition Model not found!")

    return t_obj
