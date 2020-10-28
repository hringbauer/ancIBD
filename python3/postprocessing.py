"""
Class for calling IBD from Posterior Data. Saves results as tab seperated .tsv
created by Pandas and Numpy
Contains Sub-Classes, as well as factory Method.
Pls always change parameters with set_params method!
@ Author: Harald Ringbauer, 2020
"""

import numpy as np
import pandas as pd
import os as os

class PostProcessing(object):
    """Class that can do PostProcessing of HAPSBURG output.
    (for one individual). Sometimes post-processing is done outside that,
    Has Methods to save the output. Saves using standard hapROH format"""
    folder = ""          # The Folder to operate in.
    cutoff_post = 0.9  # Cutoff Probability for ROH State
    ibd_min_l = 0.04  # Cutoff [in Morgan]
    max_gap = 0.01  # The Maximum Gap Length to be Merged [in Morgan]
    merge = True  # Whether to Merge ROH Blocks
    save = 0      # What to save. 0: Nothing 1: Save post-processed IBD. 2: Save 0-posterior. 3: Save full posterior
    output = True # Whether to plot output
    ch=1
    iid="iid0"

    def __init__(self, folder=""):
        """Initialize Class.
        Load: Whether to immediately Load the Posterior Data"""
        self.folder=folder 
        pass

    def set_params(self, **kwargs):
        """Set the Parameters.
        Takes keyworded arguments"""
        for key, value in kwargs.items():
            setattr(self, key, value)
            
    def roh_posterior(self, posterior0):
        """Load and return the posterior.
        Input: Log space [l]
        Output: Normal space [l]"""
        roh_post = 1 - np.exp(posterior0)  # Go to non-logspace probability
        return roh_post
    
    def ibd_stat_to_block(self, ibd):
        """Convert IBD status per marker
        into list of ibd.
        Input: IBD stats [l] boolean.
        Return start and end indexes of IBD blocks"""
        x1 = np.hstack([[False], ibd, [False]]).astype("int")  # padding
        d = np.diff(x1)
        starts = np.where(d == 1)[0]
        ends = np.where(d == -1)[0]
        return starts, ends
    
    def create_df(self, starts, ends, starts_map, ends_map, 
              l, l_map, iid, ch, ibd_min_l):
        """Create and returndthe hapBLOCK/hapROH dataframe."""

        full_df = pd.DataFrame({'Start': starts, 'End': ends,
                                'StartM': starts_map, 'EndM': ends_map, 'length': l,
                                'lengthM': l_map, 'iid': iid, "ch": ch})
        df = full_df[full_df["lengthM"] > ibd_min_l]  # Cut out long blocks
        return df
    
    def merge_called_blocks(self, df, max_gap=0):
        """Merge Blocks in Dataframe df and return merged Dataframe"""
        if len(df) == 0:
            return df  # In case of empty dataframe don't do anything

        if max_gap == 0:
            max_gap = self.max_gap

        df_n = df.drop(df.index)  # Create New Data frame with all raws removed
        row_c = df.iloc[0, :].copy()

        # Iterate over all rows, update blocks if gaps small enough
        for index, row in df.iterrows():
            if row["StartM"] - row_c["EndM"] < max_gap:
                row_c["End"] = row["End"]
                row_c["EndM"] = row["EndM"]
                row_c["length"] = row_c["End"] - row_c["Start"]
                row_c["lengthM"] = row_c["EndM"] - row_c["StartM"]

            else:  # Save and go to next row
                df_n.loc[len(df_n)] = row_c  # Append a row to new df
                row_c = row.copy()

        df_n.loc[len(df_n)] = row_c   # Append the last row

        if self.output == True:
            print(f"Merged n={len(df) - len(df_n)} gaps < {max_gap} M")
        return df_n
    
    def save_output(self, df, r_map=[], post=[], save_folder=""):
        """Save hapBLOCK output in standardized format."""
        if len(save_folder)==0:
            save_folder = self.folder
            
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
            if self.output:
                print(f"Created {save_folder}.")
            
        path_ibd = os.path.join(save_folder, "ibd.tsv")
        df.to_csv(path_ibd, sep="\t", index=False)
        if len(r_map)>0:
            path_map = os.path.join(save_folder, "map.tsv")
            np.savetxt(path_map, r_map, delimiter="\t")
        if len(post)>0:
            path_posterior = os.path.join(save_folder, "posterior.tsv")
            np.savetxt(path_posterior, post, delimiter="\t")
                
        if self.output:
            print(f"Successfully saved output to {save_folder}")
        
    def call_roh(self, r_map, post0):
        """Call ROH of Homozygosity from Posterior Data
        bigger than cutoff.
        post0: posterior in format [5,l], log space"""
        ibd_post = self.roh_posterior(post0[0,:])
        ibd = ibd_post > self.cutoff_post

        if self.output == True:
            frac_ibd = np.mean(ibd)
            print(f"Fraction Markers above IBD cutoff: {frac_ibd:.4f}")

        # Identify Stretches by difference (up and down)
        starts, ends = self.ibd_stat_to_block(ibd)
        l = ends - starts
        ends_map = r_map[ends - 1]  # -1 to stay within bounds
        starts_map = r_map[starts]
        l_map = ends_map - starts_map

        # Create hapROH Dataframe
        df = self.create_df(starts, ends, starts_map, ends_map, 
                            l, l_map, self.iid, self.ch, ibd_min_l=self.ibd_min_l)

        # Merge Blocks in Postprocessing Step
        if self.merge:
            df = self.merge_called_blocks(df)

        if self.output:
            print(f"Called n={len(df)} IBD Blocks > {self.ibd_min_l * 100} cM")
            l = np.max(df["lengthM"])
            print(f"Longest Block: {l *100:.2f} cM")

        if self.save==1:
            self.save_output(df)
        elif self.save==2:
            self.save_output(df, r_map=r_map, post=post0[0,:])
        elif self.save==3:
            self.save_output(df, r_map=r_map, post=post0)
        return df, r_map, post0
    
class NoPostProcessing(PostProcessing):
    """Class that does no Postprocessing"""
    
    def call_roh(*prms, **kwd):
        """Does nothing"""
        pass
    
def load_Postprocessing(p_model="hapROH"):
    """Factory Method for PostProcessing class"""
    if p_model == "hapROH":
        pp = PostProcessing()
    elif p_model == "None":
        pp = NoPostProcessing()
    else:
        raise RuntimeError(f"Postprocessing method {p_model} not available!")
    return pp