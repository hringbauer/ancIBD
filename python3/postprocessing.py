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
    cutoff_post = 0.99    # Cutoff Probability for ROH State
    min_cm1 = 6      # Cutoff of block length [in cM]
    min_cm2 = 2
    max_gap = 0.01  # The Maximum Gap Length to be Merged [in Morgan]
    snp_cm = 0.0 # the maximum gap between two consecutive SNPs [in Morgan]
    save = 0        # What to save. 0: Nothing 1: Save post-processed IBD. 2: Save 0-posterior. 3: Save full posterior
    output = True   # Whether to plot output
    ch = 3            # Chromosome to analyze
    iid = "iid0"

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
        """Calculate the IBD posterior.
        Input: Posterior 0 [l]
        Output: 1 - Posterior 0 [l]"""
        #roh_post = 1 - np.exp(posterior0)  # Go to non-logspace probability
        roh_post = 1 - posterior0
        return roh_post
    
    def ibd_stat_to_block(self, ibd, r_map):
        """Convert IBD status per marker
        into list of ibd.
        Input: IBD stats [l] boolean.
        Return start and end indexes of IBD blocks"""
        x1 = np.hstack([[False], ibd, [False]]).astype("int")  # padding
        d = np.diff(x1)
        starts = np.where(d == 1)[0]
        ends = np.where(d == -1)[0]
        
        # if self.snp_gap > 0:
            # print(f'filter segments with SNP density less than {self.snp_cm} per cM.')
            # assert(len(ibd) == len(r_map))
            # gap = np.diff(r_map)
            # assert(np.min(gap) >= 0)
            # big_gap_loc = np.where(gap > self.snp_gap)[0]
            # isIBD = np.where(ibd == 1)[0]
            # big_gap_loc = big_gap_loc[np.logical_and(np.isin(big_gap_loc, isIBD), np.isin(big_gap_loc + 1, isIBD))] # the span of the two locus is in IBD state
            # if len(big_gap_loc) > 0:
            #     print(f'found large gaps in IBD region in ch{self.ch}')
            #     print(big_gap_loc)
            #     starts = np.sort(np.append(starts, 1 + big_gap_loc))
            #     ends = np.sort(np.append(ends, big_gap_loc))

        return starts, ends
    
    def create_df(self, ch, starts, ends, starts_map, ends_map, 
              l, l_map, starts_bp, ends_bp, min_cm, iid1, iid2):
        """Create and returndthe hapBLOCK/hapROH dataframe."""

        full_df = pd.DataFrame({'iid1': iid1, "iid2": iid2, "ch": ch, 
                                'Start': starts, 'End': ends, 'length': l,
                                'StartM': starts_map, 'EndM': ends_map, 'lengthM': l_map, 
                                'StartBP': starts_bp, 'EndBP': ends_bp})
        df = full_df[full_df["lengthM"] > min_cm/100.0]  # Cut out long blocks
        return df

    def filter_by_snp_density(self, df):
        """remove segments with low SNP density"""
        print(f'filter segments with SNP density less than {self.snp_cm} per cM')
        prevN = len(df)
        df = df[(df['End'] - df['Start'])/(100*df['lengthM']) >= self.snp_cm]
        currN = len(df)
        print(f'{prevN - currN} segments removed due to low SNP density...')
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

        if self.output:
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
        self.save_ibd_df(df_ibd=df, save_path = path_ibd)

        if len(r_map)>0:
            #path_map = os.path.join(save_folder, "map.tsv")
            #np.savetxt(path_map, r_map, delimiter="\t")
            path_map = os.path.join(save_folder, "map.npy")
            np.save(path_map, r_map)
        if len(post)>0:
            path_posterior = os.path.join(save_folder, "posterior.npy")
            np.save(path_posterior, post)
                
            path_posterior = os.path.join(save_folder, "posterior.tsv")
            np.savetxt(path_posterior, post, delimiter="\t")
        if self.output:
            print(f"Successfully saved all output to {save_folder}")
            
    def save_ibd_df(self, df_ibd, save_path):
        """Save IBD Dataframe to path save_path"""
        df_ibd.to_csv(save_path, sep="\t", index=False)
        if self.output:
            print(f"Saved IBD output to: {save_path}")
        
    def call_roh(self, r_map, bp_map, post0, iid1="", iid2=""):
        """Call ROH of Homozygosity from Posterior Data
        bigger than cutoff.
        post0: posterior in format [5,l], log space"""
        ibd_post = self.roh_posterior(post0[0,:])
        ibd = ibd_post > self.cutoff_post
        
        if len(iid1)==0:
            iid1=self.iid

        if self.output:
            frac_ibd = np.mean(ibd)
            print(f"Fraction Markers above IBD cutoff: {frac_ibd:.4f}")

        # Identify Stretches by difference (up and down)
        starts, ends = self.ibd_stat_to_block(ibd, r_map)
        l = ends - starts
        ends_map = r_map[ends - 1]  # -1 to stay within bounds
        starts_map = r_map[starts]
        l_map = ends_map - starts_map
        starts_bp = bp_map[starts]
        ends_bp = bp_map[ends - 1]

        # Create hapROH Dataframe
        df = self.create_df(self.ch, starts, ends, starts_map, ends_map, 
                            l, l_map, starts_bp, ends_bp, self.min_cm1, iid1, iid2)

        # Merge Blocks in Postprocessing Step
        if self.max_gap>0:
            df = self.merge_called_blocks(df)
        # remove segments with low SNP density (likely false positives)
        if self.snp_cm>0:
            df = self.filter_by_snp_density(df)

        if self.output:
            print(f"Called n={len(df)} IBD Blocks > {self.min_cm1} cM")
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


class IBD2Postprocessing(PostProcessing):
    def create_df(self, starts, ends, starts_map, ends_map, 
              l, l_map, ch, starts_bp, ends_bp, min_cm, iid1, iid2, segment_type="IBD1"):
        """Create and returndthe hapBLOCK/hapROH dataframe."""
        segment_types = [segment_type]*len(starts)
        full_df = pd.DataFrame({'iid1': iid1, "iid2": iid2, "ch": ch,
                                'Start': starts, 'End': ends, 'length': l,
                                'StartM': starts_map, 'EndM': ends_map, 'lengthM': l_map,
                                'StartBP': starts_bp, 'EndBP': ends_bp,
                                'segment_type': segment_types})
        df = full_df[full_df["lengthM"] > min_cm/100.0]  # Cut out long blocks
        return df

    def get_posterior(self, post0):
        ibd1 = np.sum(post0[1:, :], axis=0)
        ibd2 = np.sum(post0[5:, :], axis=0)
        return ibd1, ibd2


    def call_roh(self, r_map, bp_map, post0, iid1="", iid2=""):
        """Call ROH of Homozygosity from Posterior Data
        bigger than cutoff.
        post0: posterior in format [7,l], log space"""
        assert(post0.shape[0] == 7)
        ibd1_post, ibd2_post = self.get_posterior(post0)
        ibd1 = ibd1_post > self.cutoff_post
        ibd2 = ibd2_post > self.cutoff_post
        
        if len(iid1)==0:
            iid1=self.iid

        if self.output:
            frac_ibd1 = np.mean(ibd1)
            print(f"Fraction Markers above IBD1 cutoff: {frac_ibd1:.4f}")
            frac_ibd2 = np.mean(ibd2)
            print(f'Fraction Markers above IBD2 cutoff: {frac_ibd2:.4f}')

    ############################ Writing IBD1 blocks to pandas dataframe ############################
        # Identify Stretches by difference (up and down)
        starts_ibd1, ends_ibd1 = self.ibd_stat_to_block(ibd1, r_map)
        l_ibd1 = ends_ibd1 - starts_ibd1
        ends_map_ibd1 = r_map[ends_ibd1 - 1]  # -1 to stay within bounds
        starts_map_ibd1 = r_map[starts_ibd1]
        l_map_ibd1 = ends_map_ibd1 - starts_map_ibd1
        starts_bp_ibd1 = bp_map[starts_ibd1]
        ends_bp_ibd1 = bp_map[ends_ibd1 - 1]
        # Create hapROH Dataframe
        df1 = self.create_df(starts_ibd1, ends_ibd1, starts_map_ibd1, ends_map_ibd1, 
                            l_ibd1, l_map_ibd1, self.ch, starts_bp_ibd1, ends_bp_ibd1, self.min_cm1,
                            iid1, iid2, segment_type='IBD1')
        if self.output:
            print(f"Called n={len(df1)} IBD1 Blocks > {self.min_cm1} cM")
            l = np.max(df1["lengthM"])
            print(f"Longest Block: {l *100:.2f} cM")
        # Merge Blocks in Postprocessing Step
        if self.max_gap>0:
            df1 = self.merge_called_blocks(df1)

    ########################### Writing IBD2 blocks to pandas dataframe ###############################
        starts_ibd2, ends_ibd2 = self.ibd_stat_to_block(ibd2, r_map)
        l_ibd2 = ends_ibd2 - starts_ibd2
        ends_map_ibd2 = r_map[ends_ibd2 - 1]  # -1 to stay within bounds
        starts_map_ibd2 = r_map[starts_ibd2]
        l_map_ibd2 = ends_map_ibd2 - starts_map_ibd2
        starts_bp_ibd2 = bp_map[starts_ibd2]
        ends_bp_ibd2 = bp_map[ends_ibd2 - 1]
        df2 = self.create_df(starts_ibd2, ends_ibd2, starts_map_ibd2, ends_map_ibd2, 
                            l_ibd2, l_map_ibd2, self.ch, starts_bp_ibd2, ends_bp_ibd2, self.min_cm2,
                            iid1, iid2, segment_type='IBD2')
        if self.output:
            print(f"Called n={len(df2)} IBD2 Blocks > {self.min_cm2} cM")
            l = np.max(df2["lengthM"])
            print(f"Longest Block: {l *100:.2f} cM")
        # Merge Blocks in Postprocessing Step
        if self.max_gap>0:
            df2 = self.merge_called_blocks(df2)

        # merge IBD2 blocks to the same dataframe as IBD1 blocks
        df = pd.concat((df1, df2), ignore_index=True)

        if self.save==1:
            self.save_output(df)
        elif self.save==2:
            self.save_output(df, r_map=r_map, post=post0[0,:])
        elif self.save==3:
            self.save_output(df, r_map=r_map, post=post0)
        return df, r_map, post0
    
def load_Postprocessing(p_model="hapROH"):
    """Factory Method for PostProcessing class"""
    if p_model == "hapROH":
        pp = PostProcessing()
    elif p_model == "None":
        pp = NoPostProcessing()
    elif p_model == "IBD2":
        pp = IBD2Postprocessing()
    else:
        raise RuntimeError(f"Postprocessing method {p_model} not available!")
    return pp