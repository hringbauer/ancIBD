import os as os
import numpy as np

from emission import load_emission_model
from transition import load_transition_model
from loaddata import load_loaddata
from hmm import load_fwd_bwd_func

class HMM_Full(object):
    """Analyze Class for HMMs. Wrapper for various subclasses, making them
    work together.
    Contains objects as field for pre-processing, transition as well as
    emission probabilities, loading and post-processing
    Contains most Parameters (but some of them like the output folders are decided
    by pre-processing subclass)"""
    folder_in = ""
    l_model = "" # simulated
    e_model = "" # haploid_gl
    t_model = "" # standard
    p_model = ""
    h_model = "" # FiveState FiveStateFast
    post_model = ""
    l_obj, t_obj, e_obj = 0, 0, 0 # Placeholder for the objects
    fwd_bwd = 0 # Placeholder Forward Backward Function
    submat33 = True
    in_val = 1e-4 # The initial prob. in each first/last IBD state


    def __init__(self, folder_in="", l_model="simulated", t_model="standard", 
                 e_model="haploid_gl", h_model = "FiveStateFast",
                 output=True, load=True):
        """Initialize Class.
        Keywords for which subclasses to load (via factory functions from there.
        load: Whether to load all models at intialization. True per Default."""
        self.set_params(folder_in=folder_in, l_model=l_model, 
                        t_model=t_model, e_model=e_model, 
                        h_model=h_model, output=output)
        
        if load:
            self.load_objects()

    def load_objects(self):
        """Load all the required Objects in right order"""
        self.l_obj = load_loaddata(l_model = self.l_model, path=self.folder_in)
        self.t_obj = load_transition_model(t_model = self.t_model)
        self.e_obj = load_emission_model(e_model = self.e_model)
        self.fwd_bwd = load_fwd_bwd_func(h_model=self.h_model)
        
    def run_fwd_bwd(self, output=True, full=True):
        """Run Forward Backward algorithm."""
        htsl, p, r_vec =  self.l_obj.load_all_data()
        e_mat = self.e_obj.give_emission_matrix(htsl, p)
        t_mat = self.t_obj.full_transition_matrix(r_vec, n=4, submat33=self.submat33)
        
        if full:
            post, fwd, bwd, tot_ll = self.fwd_bwd(np.log(e_mat), t_mat, in_val = self.in_val, full=full)
            return post, fwd, bwd, tot_ll
        else:
            post = self.fwd_bwd(np.log(e_mat), t_mat, in_val = self.in_val, full=full)
            return post                                                                    
                 
    def set_params(self, **kwargs):
        """Set the Parameters.
        Takes keyworded arguments"""
        for key, value in kwargs.items():
            setattr(self, key, value)