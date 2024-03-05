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
    min_cm = 4      # Cutoff of block length [in cM]
    max_gap = 0.01  # The Maximum Gap Length to be Merged [in Morgan]
    save = 0        # What to save. 0: Nothing 1: Save post-processed IBD. 2: Save 0-posterior. 3: Save full posterior
    output = True   # Whether to plot output
    ch = 3            # Chromosome to analyze
    iid = "iid0"
    mask = ""

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
              l, l_map, ch, starts_bp, ends_bp, min_cm, iid1, iid2=""):
        """Create and returndthe hapBLOCK/hapROH dataframe."""

        full_df = pd.DataFrame({'Start': starts, 'End': ends,
                                'StartM': starts_map, 'EndM': ends_map, 'length': l,
                                'lengthM': l_map, "ch": ch, 'StartBP': starts_bp, 'EndBP':ends_bp, \
                                'iid1': iid1, "iid2": iid2})
        df = full_df[full_df["lengthM"] > min_cm/100.0]  # Cut out long blocks
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
                row_c['EndBP'] = row['EndBP']

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
            path_map = os.path.join(save_folder, "map.tsv")
            np.savetxt(path_map, r_map, delimiter="\t")
        if len(post)>0:
            path_posterior = os.path.join(save_folder, "posterior.tsv")
            np.savetxt(path_posterior, post, delimiter="\t")
        if self.output:
            print(f"Successfully saved all output to {save_folder}")
            
    def save_ibd_df(self, df_ibd, save_path):
        """Save IBD Dataframe to path save_path"""
        df_ibd.to_csv(save_path, sep="\t", index=False)
        if self.output:
            print(f"Saved IBD output to: {save_path}")
    
    def mask_regions(self, df):
        df_maskTrack = pd.read_csv(self.mask, sep='\t')
        df_maskTrack = df_maskTrack[df_maskTrack['ch'] == self.ch]
        if len(df_maskTrack) == 0 or len(df) == 0:
            return df
        df2mask = df.copy()
        df_masked = []
        for _, row in df_maskTrack.iterrows():
            mask_start_bp, mask_end_bp = int(row['start_bp']), int(row['end_bp'])
            mask_start_cm, mask_end_cm = row['start_cm'], row['end_cm']
            for index, row in df2mask.iterrows():
                if row['StartBP'] >= mask_end_bp or row['EndBP'] <= mask_start_bp:
                    # no action needed
                    df_masked.append(pd.DataFrame([df2mask.loc[index]]))
                elif row['StartBP'] >= mask_start_bp and row['EndBP'] <= mask_end_bp:
                    # the whole segment is masked, drop it
                    continue
                elif row['StartBP'] < mask_start_bp and row['EndBP'] > mask_end_bp:
                    # the mask is inside the segment, split the segment into two
                    df_masked.append(pd.DataFrame([{'Start': row['Start'], 'End': row['End'],
                                      'StartM': row['StartM'], 'EndM': mask_start_cm/100,
                                      'length': row['length'],
                                      'lengthM': mask_start_cm/100 - row['StartM'], "ch": row['ch'],
                                      'StartBP': row['StartBP'], 'EndBP':mask_start_bp, \
                                      'iid1': row['iid1'], "iid2": row['iid2']}], index=[0]))
                    df_masked.append(pd.DataFrame([{'Start': row['Start'], 'End': row['End'],
                                      'StartM': mask_end_cm/100, 'EndM': row['EndM'],
                                      'length': row['length'],
                                      'lengthM': row['EndM'] - mask_end_cm/100, "ch": row['ch'],
                                      'StartBP': mask_end_bp, 'EndBP':row['EndBP'], \
                                      'iid1': row['iid1'], "iid2": row['iid2']}], index=[0]))
                elif row['StartBP'] >= mask_start_bp and row['StartBP'] < mask_end_bp:
                    # the start of the segment is inside the mask
                    df_masked.append(pd.DataFrame([{'Start': row['Start'], 'End': row['End'],
                                      'StartM': mask_end_cm/100, 'EndM': row['EndM'],
                                      'length': row['length'],
                                      'lengthM': row['EndM'] - mask_end_cm/100, "ch": row['ch'],
                                      'StartBP': mask_end_bp, 'EndBP':row['EndBP'], \
                                      'iid1': row['iid1'], "iid2": row['iid2']}], index=[0]))
                elif row['EndBP'] > mask_start_bp and row['EndBP'] <= mask_end_bp:
                    # the end of the segment is inside the mask
                    df_masked.append(pd.DataFrame([{'Start': row['Start'], 'End': row['End'],
                                      'StartM': row['StartM'], 'EndM': mask_start_cm/100,
                                      'length': row['length'],
                                      'lengthM': mask_start_cm/100 - row['StartM'], "ch": row['ch'],
                                      'StartBP': row['StartBP'], 'EndBP':mask_start_bp, \
                                      'iid1': row['iid1'], "iid2": row['iid2']}], index=[0]))
            if len(df_masked) > 0:
                #df_masked = pd.DataFrame(df_masked)
                df_masked = pd.concat(df_masked, ignore_index=True)
            else:
                df_masked = pd.DataFrame(columns=df.columns)
            print(f'size of df_masked: {len(df_masked)}')
            df2mask = df_masked.copy()
            df_masked = [] # why need this complex set-up? because one segment can be potentially masked multiple times (bug prone!!!)
        return df2mask


        
    def call_roh(self, r_map, bp, post0, iid1="", iid2=""):
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
        starts, ends = self.ibd_stat_to_block(ibd)
        l = ends - starts
        ends_map = r_map[ends - 1]  # -1 to stay within bounds
        starts_map = r_map[starts]
        l_map = ends_map - starts_map
        ends_bp = bp[ends - 1]  # -1 to stay within bounds
        starts_bp = bp[starts]

        # Create hapROH Dataframe
        df = self.create_df(starts, ends, starts_map, ends_map, 
                            l, l_map, self.ch, starts_bp, ends_bp, min_cm=self.min_cm,
                            iid1=iid1, iid2=iid2)

        # Merge Blocks in Postprocessing Step
        if self.max_gap>0:
            df = self.merge_called_blocks(df)
        # mask out regions if applicable
        
        if len(self.mask) > 0:
            print('Applying mask to IBD segments...')
            print(df)
            df = self.mask_regions(df)
            print(df)
            df = df[df['lengthM'] > self.min_cm/100.0] # remove short segments created by intersecting with masks

        if self.output:
            print(f"Called n={len(df)} IBD Blocks > {self.min_cm} cM")
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
    min_cm2_init = 1.0
    min_cm2_after_merge = 2.0
    cutoff_post1 = 0.99    # Cutoff Probability for IBD1 state
    cutoff_post2 = 0.975    # Cutoff Probability for IBD2 state

    def create_df(self, starts, ends, starts_map, ends_map, 
              l, l_map, starts_bp, ends_bp, ch, min_cm, iid1, iid2, segment_type="IBD1"):
        """Create and returndthe hapBLOCK/hapROH dataframe."""
        segment_types = [segment_type]*len(starts)
        full_df = pd.DataFrame({'iid1': iid1, "iid2": iid2, "ch": ch,
                                'Start': starts, 'End': ends, 'length': l,
                                'StartM': starts_map, 'EndM': ends_map, 'lengthM': l_map,
                                'StartBP':starts_bp, 'EndBP':ends_bp, \
                                'segment_type': segment_types})
        df = full_df[full_df["lengthM"] > min_cm/100.0]  # Cut out long blocks
        return df

    def get_posterior(self, post0):
        ibd1 = np.sum(post0[1:, :], axis=0)
        ibd2 = np.sum(post0[5:, :], axis=0)
        return ibd1, ibd2


    def call_roh(self, r_map, bp, post0, iid1="", iid2=""):
        """Call ROH of Homozygosity from Posterior Data
        bigger than cutoff.
        post0: posterior in format [7,l], log space"""
        assert(post0.shape[0] == 7)
        ibd1_post, ibd2_post = self.get_posterior(post0)
        ibd1 = ibd1_post > self.cutoff_post1
        ibd2 = ibd2_post > self.cutoff_post2
        
        if len(iid1)==0:
            iid1=self.iid

        if self.output:
            frac_ibd1 = np.mean(ibd1)
            print(f"Fraction Markers above IBD1 cutoff {self.cutoff_post1}: {frac_ibd1:.4f}")
            frac_ibd2 = np.mean(ibd2)
            print(f'Fraction Markers above IBD2 cutoff {self.cutoff_post2}: {frac_ibd2:.4f}')

    ############################ Writing IBD1 blocks to pandas dataframe ############################
        # Identify Stretches by difference (up and down)
        starts_ibd1, ends_ibd1 = self.ibd_stat_to_block(ibd1)
        l_ibd1 = ends_ibd1 - starts_ibd1
        ends_map_ibd1 = r_map[ends_ibd1 - 1]  # -1 to stay within bounds
        starts_map_ibd1 = r_map[starts_ibd1]
        starts_bp_ibd1 = bp[starts_ibd1]
        ends_bp_ibd1 = bp[ends_ibd1 - 1]  # -1 to stay within bounds
        l_map_ibd1 = ends_map_ibd1 - starts_map_ibd1

        # Create hapROH Dataframe
        df1 = self.create_df(starts_ibd1, ends_ibd1, starts_map_ibd1, ends_map_ibd1, 
                            l_ibd1, l_map_ibd1, starts_bp_ibd1, ends_bp_ibd1, self.ch, self.min_cm,
                            iid1, iid2, segment_type='IBD1')
        if self.output:
            print(f"Called n={len(df1)} IBD1 Blocks > {self.min_cm} cM")
            l = np.max(df1["lengthM"])
            print(f"Longest Block: {l *100:.2f} cM")
        # Merge Blocks in Postprocessing Step
        if self.max_gap>0:
            df1 = self.merge_called_blocks(df1)

    ########################### Writing IBD2 blocks to pandas dataframe ###############################
        starts_ibd2, ends_ibd2 = self.ibd_stat_to_block(ibd2)
        l_ibd2 = ends_ibd2 - starts_ibd2
        ends_map_ibd2 = r_map[ends_ibd2 - 1]  # -1 to stay within bounds
        starts_map_ibd2 = r_map[starts_ibd2]
        starts_bp_ibd2 = bp[starts_ibd2]
        ends_bp_ibd2 = bp[ends_ibd2 - 1]  # -1 to stay within bounds
        l_map_ibd2 = ends_map_ibd2 - starts_map_ibd2

        # here this is to call all segments before merge, so we use min_cm2_init
        df2 = self.create_df(starts_ibd2, ends_ibd2, starts_map_ibd2, ends_map_ibd2, 
                            l_ibd2, l_map_ibd2, starts_bp_ibd2, ends_bp_ibd2, self.ch, self.min_cm2_init,
                            iid1, iid2, segment_type='IBD2')
        if self.output:
            print(f"Called n={len(df2)} IBD2 Blocks > {self.min_cm2} cM")
            l = np.max(df2["lengthM"])
            print(f"Longest Block: {l *100:.2f} cM")
        # Merge Blocks in Postprocessing Step
        if self.max_gap>0:
            df2 = self.merge_called_blocks(df2)
        # after merge, filter to segments longer than min_cm2_after_merge
        df2 = df2[100*df2['lengthM'] > self.min_cm2_after_merge]

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
    elif p_model == 'IBD2':
        pp = IBD2Postprocessing()
    else:
        raise RuntimeError(f"Postprocessing method {p_model} not available!")
    return pp