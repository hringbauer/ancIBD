import numpy as np
import pandas as pd
import socket as socket
import os as os
import sys as sys
import h5py


def get_af(f, min_gp=0.99):
    """Get Allele Frequency"""
    gp_max = np.max(f["calldata/GP"], axis=2)
    gp_good = (gp_max>=min_gp) # The decent genotype probabilitiies
    gp_max = 0 # Delete GP max (unnecessary now)
    
    gt1 = np.sum(f["calldata/GT"], axis=2)/2.0 # Get the genotype sum
    gp_good_c = np.sum(gp_good, axis=1)
    af = np.sum(gt1 * gp_good, axis=1) / gp_good_c
    return af

def get_af1000G(f):
    """Get Allele Frequency - ASSUME ALL GT are set!"""
    gt1 = np.sum(f["calldata/GT"], axis=2)/2.0 # Get the genotype sum
    af = np.mean(gt1, axis=1)
    return af

def merge_in_af(path_h5, af, col_af="AF_ALL"):
    """Merge in AF into hdf5 file. Save modified h5 in place 
    af: Array of allele frequencies to save"""
    
    ### Now create the new column in hdf5
    print("Adding map to HDF5...")
    with h5py.File(path_h5, 'a') as f0:
        group = f0["variants"]
        l = len(f0["variants/POS"]) # Get number of markers
        print(f"Loaded {len(af)} variants.")
        assert(l==len(af)) # Sanity Checks
        assert(np.min(af)>=0)
        assert(np.max(af)<=1)
        
        group.create_dataset(col_af, (l,), dtype='f')   
        f0[f"variants/{col_af}"][:] = af[:]
    print(f"Finshed merged in allele frequencies into {path_h5}")
    
#################################################################
### Bring over AF

def lift_af(h5_target, h5_original, field="variants/AF_ALL", 
            match_col="variants/POS", dt=np.float):
    """Bring over field from one h5 to another. Assume field does not exist in target
    h5_original: The original hdf5 path
    h5_target: The target hdf5 path
    field: Which fielw to copy over 
    """
    g= h5py.File(h5_original, "r") # To read ref data
    f = h5py.File(h5_target, 'a') # To append data
    
    l = len(f[match_col]) #  nr target loci
    p = 0.5 * np.ones(l, dtype=np.float) # Default value for p
    
    ### Match on position and lift over.
    its, i1, i2 = np.intersect1d(g[match_col][:], f[match_col][:], return_indices=True)
    p[i2] = g[field][:][i1]
    print(f"Intersection {len(i2)} out of {l} target HDF5 SNPs")
    
    ### Crate and write the field in target h5
    if not field in f:
        f.create_dataset(field, (l,), dtype='f')  
    else:
        print(f"Found pre-existing field.")
    f[field][:] = p
    
    print("We did it. Finished.")
    return