"""
Methods for fwd_bwd calculations. With Factory methods.
Be careful with c functions, they need exact type input.
Contains Sub-Classes, as well as factory Method.
@ Author: Harald Ringbauer, 2019, All rights reserved
"""
import numpy as np
from scipy.special import logsumexp as logsumexp_p
from cfunc import fwd_bkwd_fast, fwd_bkwd_lowmem

def fwd_bkwd(e_mat0, t_mat, in_val = 1e-4, full=False):
    """Takes emission and transition probabilities, and calculates posteriors.
    Uses speed-up specific for symmetric states 1...n (pooling same transition rates)
    Low-Mem: Do no save the full FWD BWD and Posterior. Use temporary
    Arrays for saving.
    Input:
    e_mat0: Emission probabilities [k x l] (log space)       (log space)
    t_mat: Transition Matrix: [l x 3 x 3]                     (normal space)
    in_val: Initial Probability of being in IBD state
    full: Boolean whether to return everything"""
    n_states = e_mat0.shape[0]
    n_loci = e_mat0.shape[1]
    stay, tot_ll = 0.,0.  #e Probablility of Staying

    # Initialize Posterior and Transition Probabilities
    post = np.empty((n_loci,n_states), dtype=np.float) # Array of 0 State Posterior
    trans_ll = np.empty(n_states-1, dtype=np.float) # Array for pre-calculations

    three_v = np.empty(3, dtype=np.float)     # Array of size three
    two_v = np.empty(2, dtype=np.float)       # Array of size two

    ### Initialize FWD BWD Arrays
    fwd0 = np.zeros(n_states, dtype=np.float)
    fwd0[:] = np.log(in_val)  # Initial Probabilities
    fwd0[0] = np.log(1 - (n_states - 1) * in_val)
    #cdef double[:] fwd = fwd0

    bwd0 = np.zeros(n_states, dtype=np.float)
    bwd0[:] = np.log(in_val)
    bwd0[0] = np.log(1 - (n_states - 1) * in_val)
    #cdef double[:] bwd = bwd0

    tmp = np.zeros(n_states, dtype=np.float)
    #cdef double[:] tmp = tmp0
    
    # Do transform to Log Space:
    t0 = np.log(t_mat)      

    #############################
    ### Do the Forward Algorithm

    post[0,:] = fwd0 # Add to first locus 0 Posterior
    for i in range(1, n_loci):  # Run forward recursion
        stay = np.log(t_mat[i, 1, 1] - t_mat[i, 1, 2])  # Do the log of the Stay term

        for k in range(1, n_states): # Calculate logsum of ROH states:
            trans_ll[k-1] = fwd0[k]
        f_l = logsumexp_p(trans_ll) # Logsum of ROH States

        # Do the 0 State:
        two_v[0] = fwd0[0] + t0[i, 0, 0]   # Staying in 0 State
        two_v[1] = f_l + t0[i, 1, 0]             # Going into 0 State
        tmp[0] = e_mat0[0, i] + logsumexp_p(two_v)

        ### Do the other states
        # Preprocessing:
        three_v[0] = fwd0[0] + t0[i, 0, 1]   # Coming from 0 State
        three_v[1] = f_l + t0[i, 1, 2]             # Coming from other ROH State

        for j in range(1, n_states):  # Do the final run over all states
            three_v[2] = fwd0[j] + stay
            tmp[j] = e_mat0[j, i] + logsumexp_p(three_v)

        ### Make tmp new fwd vec:
        fwd0 = tmp
        post[i,:] = fwd0  # Add to 0-State Posterior

    ### Get total log likelihood
    tot_ll = logsumexp_p(fwd0+bwd0)

    #############################
    ### Do the Backward Algorithm
    ## last0-State Posterior
    post[n_loci-1,:] = post[n_loci-1,:] + bwd0[:] - tot_ll

    for i in range(n_loci-1, 0, -1):  # Run backward recursion
        stay = np.log(t_mat[i, 1, 1] - t_mat[i, 1, 2])

        for k in range(1, n_states): # Calculate logsum of ROH states:
            trans_ll[k-1] = bwd0[k] + e_mat0[k, i]
        f_l = logsumexp_p(trans_ll) # Logsum of ROH States

        # Do the 0 State:
        two_v[0] = bwd0[0] + t0[i, 0, 0] + e_mat0[0, i]   # Staying in 0 State
        two_v[1] = f_l + t0[i, 0, 1]                         # Going into 0 State
        tmp[0] = logsumexp_p(two_v)

        ### Do the other states
        # Preprocessing:
        three_v[0] = e_mat0[0, i] + bwd0[0] + t0[i, 1, 0]
        three_v[1] = f_l + t0[i, 1, 2]    # Coming from other ROH State

        for j in range(1, n_states):  # Do the final run over all states
            three_v[2] = e_mat0[j, i] + bwd0[j] +  stay
            tmp[j] = logsumexp_p(three_v)  # Fill in the backward Probability

        ### Make tmp new bwd vec:
        bwd0 = tmp

        ### Finalize the 0 Posterior
        post[i-1,:] = post[i-1,:] + bwd0[:] - tot_ll

    print(f"Total Log likelihood: {tot_ll: .3f}")
    #print_memory_usage()   ## For MEMORY_BENCH

    if full==False:
        return post[:,:]  # For "fake" axis

    elif full==True:   # Return everything
        return post[:,:], fwd0, bwd0, tot_ll

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
        func = fwd_bkwd
    elif h_model == "FiveStateFast":
        func = fwd_bkwd_fast
    elif h_model == "FiveStatLowMem":
        func = fwd_bkwd_lowmem
    else:
        raise NotImplementedError("Transition Model not found!")

    return func