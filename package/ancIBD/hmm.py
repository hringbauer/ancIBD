"""
Methods for fwd_bwd calculations. With Factory methods.
Be careful with c functions, they need exact type input.
Contains Sub-Classes, as well as factory Method.
@ Author: Harald Ringbauer, 2019, All rights reserved
"""
import numpy as np
from scipy.special import logsumexp as logsumexp
from ancIBD.cfunc import fwd_bkwd_fast, fwd_bkwd_lowmem, fwd_bkwd_scaled

####################################################
### Python function FWD BWD. (C FUNCTIONS GET IMPORTED)

def fwd_bkwd_p(e_prob, t_mat,  in_val = 1e-4, 
               full=True, output=True):
    """Simple Python Implementation of forward / bacdkward to compare 
    optimized c versions to. Uses standard definitions.
    Takes emission and transition probabilities
    Input: 
    Emission Matrix e_prob0 [l,k] NORMAL SPACE
    Transition Matrix t_mat [l,k,k] NORMALIZED
    and initialized fwd and bwd probabilities. All in log Space
    output: Whether to print out output useful for debugging/monitoring
    Output: post, fwd0, bwd0, tot_ll; if full then only post"""
    e_prob0 = np.log(e_prob)
    n_states, n_loci = np.shape(e_prob0)
    
    ### Initialize FWD BWD and Posterior
    fwd0 = np.zeros((n_states, n_loci), dtype="float")
    t_mat0 = np.log(t_mat)
    
    fwd0[:, 0] = np.log(in_val)  # Initial Probabilities
    fwd0[0, 0] = np.log(1 - (n_states - 1) * in_val)

    bwd0 = np.zeros((n_states, n_loci), dtype="float")
    bwd0[:, -1] = np.log(in_val)
    bwd0[0, -1] = np.log(1 - (n_states - 1) * in_val)
    
    ### Forward Run
    for i in range(1, n_loci): 
        for j in range(n_states):
            trans_ll = fwd0[:, i - 1] + t_mat0[i, :, j]
            fwd0[j, i] = e_prob0[j, i] + logsumexp(trans_ll)
    
    ### Backward Run
    for i in range(n_loci - 1, 0, -1):  # Do the backward recursion
        for j in range(n_states):
            trans_ll = t_mat0[i, j, :] + e_prob0[:, i] + bwd0[:, i]
            bwd0[j, i - 1] = logsumexp(trans_ll)

    tot_ll = logsumexp(fwd0[:, -1] + bwd0[:, -1])  # Get total log likelihood from last entry
    
    if output:
        print(f"Total Log likelihood fwd: {tot_ll: .3f}")
    # Combine the forward and backward calculations
    post = fwd0 + bwd0 - tot_ll  # The formulat is f*b/tot_l
    
    if full:
        return post, fwd0, bwd0, tot_ll
    else:
        return post        

####################################################
####################################################
### Additional Helper Functions

def print_memory_usage():
    """Print the current Memory Usage in mB"""
    process = psutil.Process(os.getpid())
    mb_usage = process.memory_info().rss / 1e6
    print(f"Memory Usage: {mb_usage} mB")
    
    
####################################################
### Factory method to load the right function

def load_fwd_bwd_func(h_model="FiveState"):
    """Return fwd_bwd function"""
    if h_model == "FiveState":
        func = fwd_bkwd_p
    elif h_model == "FiveStateFast":
        func = fwd_bkwd_fast
    elif h_model == "FiveStatLowMem":
        func = fwd_bkwd_lowmem
    elif h_model == "FiveStateScaled":
        func = fwd_bkwd_scaled
    else:
        raise NotImplementedError("HMM Model not found!")
    return func