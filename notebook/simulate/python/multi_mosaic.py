"""
Class for preparing 1000 Genome Mosaic Data
Wrapper for createMosaic - to simulate several indiviuals with several IBD
Save the created Data to HDF 5, and the ROH block info a csv in the same folder.
@ Author: Harald Ringbauer, 2020, All rights reserved
"""

import h5py                             # For Processing HDF5s
import numpy as np
import pandas as pd
import os                               # To delete File
from python.create_mosaic import Mosaic_1000G  # custom import

class Mosaic_1000G_Multi(object):
    """Class for Preparing Multiple 1000G Mosaics. And Saving them.
    A wrapper for createMosaic"""

    # Important Parameters:
    ch = 3                              # Which Chromosome to analyze
    pop_list = ["TSI"]                  # Which pop to copy from ["TSI"]
    chunk_length = 0.005                # Chunk Length of Chromosomes 0.005
    ibd_lengths = np.ones(5) * 0.05     # IBD Lengths per Individual
    iid = "iid"                         # Prefix of Artificial Individual Names
    n = 3       # Nr of individuals to simulate
    gp = True   # Whether to save genotype probabilities in hdf5

    path1000G = "/n/groups/reich/hringbauer/git/hapBLOCK/data/hdf5/1240k_1000G/chr"
    pop_path = "/n/groups/reich/hringbauer/git/hapBLOCK/data/hdf5/1240k_1000G/meta_df_all.csv"
    # Where to save the new HDF5 to by default
    save_path = "/n/groups/reich/hringbauer/git/hapBLOCK/output/simulated/TSI/"

    output = True  # whether to Print Output
    m_object = 0  # The Mosaic object

    def __init__(self):
        """Initialize"""
        pass  # Just go on

    def load_m_object(self):
        """Load the Mosaic Object"""
        print(self.save_path)  # Create Save folder if needed
        if not os.path.exists(self.save_path):
            os.makedirs(self.save_path)
        self.m_object = Mosaic_1000G(ch=self.ch, path1000G=self.path1000G,
                                     pop_path=self.pop_path, save_path=self.save_path)

    def save_hdf5(self, gt, ad, ref, alt, pos, rec, samples, path, gp=[], compression="gzip"):
        """Create a new HDF5 File with Input Data.
        gt: Genotype data [l,k,2]
        ad: Allele depth [l,k,2]
        ref: Reference Allele [l]
        alt: Alternate Allele [l]
        pos: Position  [l]
        m: Map position [l]
        samples: Sample IDs [k]"""

        l, k, _ = np.shape(gt)  # # loci and # of Individuals

        if os.path.exists(path):  # Delete existing File if there
            os.remove(path)

        dt = h5py.special_dtype(vlen=str)  # To have no problem with saving

        with h5py.File(path, 'w') as f0:
            # Create all the Groups
            f_map = f0.create_dataset("variants/MAP", (l,), dtype='f')
            f_ad = f0.create_dataset("calldata/AD", (l, k, 2), dtype='i', compression=compression)
            f_ref = f0.create_dataset("variants/REF", (l,), dtype=dt)
            f_alt = f0.create_dataset("variants/ALT", (l,), dtype=dt)
            f_pos = f0.create_dataset("variants/POS", (l,), dtype='i')
            f_gt = f0.create_dataset("calldata/GT", (l, k, 2), dtype='i', compression=compression)
            if len(gp)>0:
                f_gp = f0.create_dataset("calldata/GP", (l, k, 3), dtype="f", compression=compression)    
            f_samples = f0.create_dataset("samples", (k,), dtype=dt)

            # Save the Data
            f_map[:] = rec
            f_ad[:] = ad
            f_ref[:] = ref.astype("S1")
            f_alt[:] = alt.astype("S1")
            f_pos[:] = pos
            f_gt[:] = gt
            if len(gp)>0:
                f_gp[:] = gp
            f_samples[:] = np.array(samples).astype("S10")

        if self.output:
            print(f"Successfully saved {k} individuals to: {path}")

    def create_individuals(self):
        """Creates and saves Genotypes for sets of individuals.
        Parameters are loaded from the Class
        Save Individuals as well as IBD positions (in Pandas Dataframe)"""
        m = self.m_object # The mosaic object
        chunk_length = self.chunk_length
        pop_list = self.pop_list
        n = self.n  # number of IBD pairs to simulate
        iid = self.iid
        ibd_lengths = self.ibd_lengths

        l = len(m.f["variants/MAP"])
        c_min, c_max = np.min(m.f["variants/MAP"]), np.max(m.f["variants/MAP"])

        gts = -np.ones((l, 2*n, 2), dtype="int8") # Twice the individuals
        
        iids =np.empty(2*n, dtype="<U10")
        iids[::2] = [(iid + str(i) + "A") for i in range(n)]  # Create IIDs
        iids[1::2] = [(iid + str(i)+ "B") for i in range(n)]  # Create IIDs

        # Make lists for IBD dataframe
        ibd_begins, ibd_ends = [], []
        iid_list1, iid_list2, copy_iids = [], [], []

        for i in range(n):  # Iterate over all individuals
            if self.output:
                print(f"\nDoing Individual {iids[2*i]} and {iids[2*i+1]}")
            
            ### Create IID List and Genotypes
            ibd_list = self.create_ibd_list(ibd_lengths, c_min, c_max)
            gts_ibd, copy_ids = m.create_chunked_ibd_individuals(
                chunk_length=chunk_length, ibd_list=ibd_list, pop_list=pop_list)
            gts[:, 2*i, :] = gts_ibd[0]
            gts[:, 2*i+1, :] = gts_ibd[1] # Save in batches of two
            

            # Append IBD information to save
            ibd_begins += list(ibd_list[:, 0])
            ibd_ends += list(ibd_list[:, 1])
            iid_list1 += [iids[2*i], ] * len(ibd_list)
            iid_list2 += [iids[2*i+1], ] * len(ibd_list)
            copy_iids += list(copy_ids)
        
        assert(np.min(gts)>=0) # Sanity Check.
        if self.output:
            print(f"Finished creating genotypes. Dimension: {np.shape(gts)}")

        ### Save data:
        self.save_ibdlist(ibd_begins, ibd_ends, iid_list1, iid_list2,
                          copy_iids, ch=self.ch)

        self.save_genotypes(m.f, gts, iids)

    def extract_individuals(self):
        """Extracts and saves Genotypes as new hdf5 of a number of Individuals,
        with exactly the same genotypes
        identical state
        """
        m = self.m_object
        pop_list = self.pop_list

        gts, iids = m.get_gts_pop(pop_list)

        self.save_genotypes(m.f, gts, iids)  # Save the genotypes

    def create_ibd_list(self, ibd_lengths, min_ch, max_ch):
        """Create List of IBD positions to copy in [i,2]
        Evenly spaces Chromosome and places IBD randomly in there
        (enforcing gaps at least length of ibds)
        ibd_lengths: How long the IBDs are [In Morgan]
        min_ch: Minimum Map Position
        max_ch: Maximum Max Position"""
        k = len(ibd_lengths)

        if k == 0:   # If no blocks, return empty list
            return np.array([[], []]).T

        l = (max_ch - min_ch) / k  # Length Batch

        # Factor 2: Make sure that there is some gap
        pos = np.random.random(size=k) * (l - ibd_lengths * 2)
        ibd_begin = pos + np.arange(k) * l + min_ch
        ibd_end = ibd_begin + ibd_lengths

        ibd_list = np.array([ibd_begin, ibd_end]).T
        return ibd_list
    
    def gp_from_gts(self, gts):
        """Create GP [l,k,3] from
        genotypes [l,k,2], with prob.
        of genotype set to 1"""
        gs = np.sum(gts, axis=2)
        l, k = np.shape(gs)
        gp = np.zeros((l,k,3), dtype="f")
        gp[gs==0,0]=1
        gp[gs==1,1]=1
        gp[gs==2,2]=1
        #(np.sum(gp, axis=2)==1).all()
        return gp

    def save_genotypes(self, f, gts, samples, path=""):
        """Save the full genotype Matrix
        f: HDF5 File with Matching entries for Meta Data per Locus
        gts: Genotype Matrix [l,k,2]
        samples: List of sample IIDs [k]
        Path: Where to save to"""
        if len(path) == 0:
            path = self.save_path
        ch  = self.ch
        path = self.save_path + f"sim_ch{ch}.h5"

        gt = gts
        ad = gts

        ref, alt = f["variants/REF"][:], f["variants/ALT"][:, 0]
        pos = f["variants/POS"]
        rec = f["variants/MAP"]
        
        if self.gp:
            gp = self.gp_from_gts(gts) # Get gentoype probabilities
        else:
            gp = [] # Don't save it

        # Maybe filter for Loci here
        self.save_hdf5(gt, ad, ref, alt, pos, rec, samples, path, gp=gp)

    def save_ibdlist(self, ibd_beg, ibd_end, iids1, iids2, copyiids, ch=0, path=""):
        """Save the ROH List to Path.
        iids: Individual IIDs to save
        copyiid: Which Individuals the ibd was copied from"""
        if len(path) == 0:
            path = self.save_path

        if ch == 0:
            ch = self.ch

        # Create the Saving Dataframe
        df_save = pd.DataFrame({"IBD_Begin": ibd_beg,
                                "IBD_End": ibd_end,
                                "iid1": iids1,
                                "iid2": iids2,
                                "copyiid": copyiids,
                                "chr": ch})

        path = path + "ibd_info.csv"
        df_save.to_csv(path, sep="\t", index=False)

        if self.output == True:
            print(f"Successfully saved to {path}")

    def save_metafile(self, path=""):
        """Save the population File"""
        # if len(path)==0:
        #    path = self.path

        save_df = self.m_object.meta_df
        save_df.to_csv(path, sep="\t", index=False)
        print(f"Successfully saved to {path}")

####################################################################
####################################################################
### Code to simulate copied in IBD

def multi_run_lengths(base_path="./Simulated/1000G_Mosaic/TSI/", pop_list=["TSI"], n=20,
                      ch=3, chunk_length=0.005, lengths=[1.0, 3.0, 5.0, 10.0], n_blocks=5):
    """Create Multiple IBD runs and saves combined data into base_path hdf5 and ibd_info df
    base_path:  Start of SavePaths
    pop_list: The Reference Populations for Mosaic
    n: Number of IBD pairs to simulate [2x individuals]
    chunk_length: Lenths of the Chunks to mosaic
    ch: Chromosome to use
    lengths: Length of the Blocks to copy in [in cM]
    n_blocks: Nr of the Blocks to copy in"""

    t = Mosaic_1000G_Multi()  # Create the MultiRUn Object
    t.pop_list = pop_list
    t.n = n
    t.chunk_length = chunk_length
    t.ch = ch  # The Chromosome

    for l in lengths:
        t.ibd_lengths = np.ones(n_blocks) * 0.01 * l  # Set the Lengths [in Morgan]
        # Where to save the new HDF5 to
        t.save_path = base_path + "ch" + str(t.ch) + "_" + str(int(l)) + "cm/"
        t.load_m_object()
        t.create_individuals()


def copy_population(base_path="./Simulated/1000G_Mosaic/TSI0/",
                    path1000G="./Data/1000Genomes/HDF5/1240kHDF5/Eur1240chr",
                    pop_list=["TSI", ], ch=3):
    """Extract one or more population(s)
    base_path: Where to save extracted data to
    path1000G: Where to take the reference hdf5 from
    pop_list: The copied Populations"""
    t = Mosaic_1000G_Multi()  # Create the MultiRUn Object
    t.pop_list = pop_list
    t.path1000G = path1000G
    t.ch = ch
    t.save_path = base_path + "ch" + str(t.ch) + "/"
    t.load_m_object()  # Load the inner Mosaic Object
    t.extract_individuals()