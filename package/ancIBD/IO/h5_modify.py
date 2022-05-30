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
    print("Adding Allele Frequencies to HDF5...")
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
            match_col="variants/POS", dt=np.float, p_def=0.5):
    """Bring over field from one h5 to another. Assume field does not exist in target
    h5_original: The original hdf5 path
    h5_target: The target hdf5 path
    field: Which fielw to copy over 
    p_def: Default Value of allele frequency
    """
    g= h5py.File(h5_original, "r") # To read ref data
    f = h5py.File(h5_target, 'a') # To append data
    
    l = len(f[match_col]) #  nr target loci
    p = p_def * np.ones(l, dtype=dt) # Default value for p
    
    ### Match on position and lift over.
    its, i1, i2 = np.intersect1d(g[match_col][:], f[match_col][:], return_indices=True)
    p[i2] = g[field][:][i1]
    print(f"Intersection {len(i2)} out of {l} target HDF5 SNPs")
    
    ### Crate and write the field in target h5
    if not field in f:
        f.create_dataset(field, (l,), dtype=dt)  
    else:
        print(f"Found pre-existing field.")
    f[field][:] = p
    
    f.close()   ### Close hdf5 files
    g.close()
    
def lift_af_df(h5_target, path_df, field="variants/AF_ALL",
               match_col="variants/POS", dt=np.float, p_def=0.5):
    """Load allele frequencies from dataframe at path_df [string]
    and merge into hdf5 file at h5_target [string] at field [string]
    Match positions on match_col [string]"""
    
    df = pd.read_csv(path_df, sep="\t") # Load Dataframe
    
    ### Load Target HDF5 
    f = h5py.File(h5_target, 'a') # Load target hdf5
    l = len(f[match_col]) #  nr target loci
    p = p_def * np.ones(l, dtype=dt) # Default value for p
    
    its, i1, i2 = np.intersect1d(df["pos"][:], f[match_col][:], return_indices=True)
    p[i2] = df["af"][:][i1]
    print(f"Intersection {len(i2)} out of {l} target HDF5 SNPs. {l-len(i2)} SNPs set to AF={p_def}")
    
    ### Crate and write the field in target h5
    if not field in f:
        f.create_dataset(field, (l,), dtype=dt)  
    else:
        print(f"Found pre-existing field.")
    f[field][:] = p
    
    f.close()  
    
    
#################################################################

def merge_in_ld_map(path_h5, path_snp1240k, chs=range(1,23), write_mode="a"):
    """Merge in MAP from eigenstrat .snp file into
    hdf5 file. Save modified h5 in place 
    path_h5: Path to hdf5 file to modify.
    path_snp1240k: Path to Eigenstrat .snp file whose map to use
    chs: Which Chromosomes to merge in HDF5 [list].
    write_mode: Which mode to use on hdf5. a: New field. r+: Change Field"""
    with h5py.File(path_h5, "r") as f:
        print("Lifting LD Map from eigenstrat to HDF5...")
        print("Loaded %i variants." % np.shape(f["calldata/GT"])[0])
        print("Loaded %i individuals." % np.shape(f["calldata/GT"])[1])

        ### Load Eigenstrat
        df_snp = pd.read_csv(path_snp1240k, header=None, sep=r"\s+", engine="python")
        df_snp.columns = ["SNP", "chr", "map", "pos", "ref", "alt"]

        rec = np.zeros(len(f["variants/POS"]))  # Create the array for vector

        for ch in chs:
            df_t = df_snp[df_snp["chr"] == ch]
            print(f"Loaded {len(df_t)} Chr.{ch} 1240K SNPs.")

            idx_f = f["variants/CHROM"][:].astype("str")==str(ch)
            if np.sum(idx_f)==0:  # If no markers found jump to next chromosome
                print("Did not find any markers...")
                continue
            rec_ch = np.zeros(len(idx_f), dtype="float")

            ### Intersect SNP positions
            its, i1, i2 = np.intersect1d(f["variants/POS"][idx_f], df_t["pos"], return_indices=True)

            l = np.sum(idx_f)
            print(f"Intersection {len(i2)} out of {l} HDF5 SNPs")

            ### Extract Map positions
            rec_ch[i1] = df_t["map"].values[i2]  # Fill in the values in Recombination map

            ### Interpolate if Needed (map position still 0)
            itp_idx = (rec_ch == 0)
            if np.sum(itp_idx) > 0:   # In case we have to interpolate
                print(f"Interpolating {np.sum(itp_idx)} variants.")
                x = df_t["pos"] 
                y = df_t["map"]   
                x1 = f["variants/POS"][:][idx_f]  # Extract all positions of interest
                rec_ch = np.interp(x1, x, y) 
            
            ### Make sure that sorted
            assert(np.all(np.diff(rec_ch)>=0))  # Assert the Recombination Map is sorted! (no 0 left and no funky stuff)
            rec[idx_f]=rec_ch # Set the Map position for chromosome indices
            print(f"Finished Chromosome {ch}.")
    
    ### Now create the new column in hdf5
    print("Adding map to HDF5...")
    with h5py.File(path_h5, write_mode) as f0:
        group = f0["variants"]
        l = len(f0["variants/POS"])
        if write_mode == "a":  # If appending new data
            group.create_dataset('MAP', (l,), dtype='f')   
        f0["variants/MAP"][:] = rec[:]

#################################################################
### Save HDF5


def save_h5(gt, ad, ref, alt, pos, 
            rec, samples, path, gp=[],
            compression="gzip", ad_group=True, gt_type="int8"):
    """Create a new HDF5 File with Input Data.
    gt: Genotype data [l,k,2]
    ad: Allele depth [l,k,2]
    ref: Reference Allele [l]
    alt: Alternate Allele [l]
    pos: Position  [l]
    m: Map position [l]
    samples: Sample IDs [k].
    Save genotype data as int8, readcount data as int16.
    ad: whether to save allele depth
    gt_type: What genotype data type save"""

    l, k, _ = np.shape(gt)  # Nr loci and Nr of Individuals

    if os.path.exists(path):  # Do a Deletion of existing File there
        os.remove(path)

    dt = h5py.special_dtype(vlen=str)  # To have no problem with saving
    with h5py.File(path, 'w') as f0:
        ### Create all the Groups
        f_map = f0.create_dataset("variants/MAP", (l,), dtype='f')
        if ad_group:
            f_ad = f0.create_dataset("calldata/AD", (l, k, 2), dtype='int8', compression=compression)
        f_ref = f0.create_dataset("variants/REF", (l,), dtype=dt)
        f_alt = f0.create_dataset("variants/ALT", (l,), dtype=dt)
        f_pos = f0.create_dataset("variants/POS", (l,), dtype='int32')
        f_gt = f0.create_dataset("calldata/GT", (l, k, 2), dtype=gt_type, compression=compression)
        if len(gp)>0:
            f_gp = f0.create_dataset("calldata/GP", (l, k, 3), dtype="f", compression=compression)     
        f_samples = f0.create_dataset("samples", (k,), dtype=dt)

        ### Save the Data
        f_map[:] = rec
        if ad_group:
            f_ad[:] = ad
        f_ref[:] = ref.astype("S1")
        f_alt[:] = alt.astype("S1")
        f_pos[:] = pos
        f_gt[:] = gt
        if len(gp)>0:
            f_gp[:] = gp
        f_samples[:] = np.array(samples).astype("S10")

    print(f"Successfully saved {k} individuals to: {path}")