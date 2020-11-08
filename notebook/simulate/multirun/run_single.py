import pandas as pd
import numpy as np
import sys as sys
import os as os
sys.path.insert(0,"/n/groups/reich/hringbauer/git/hapBLOCK/notebook/simulate/")  # hack to get specific path
from python.multi_mosaic import multi_run_lengths

### parameters are defined below in main function            
#########################################################
#########################################################

cms = [4, 8, 12, 16, 20]
n_blocks=4
n = 100 # How many Individuals to simulate
chunk_length=0.0025
ch=3
base_path = "/n/groups/reich/hringbauer/git/hapBLOCK/output/simulated/TSI/"
pop_list=["TSI", ]


#########################################################
#########################################################

if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Script needs exactly 1 argument.")
    
    ### Load correct set of IIDs
    i = int(sys.argv[1])
    cm = cms[i]
    
    #########################################################
    multi_run_lengths(base_path="/n/groups/reich/hringbauer/git/hapBLOCK/output/simulated/TSI/", pop_list=["TSI", ],
                      lengths=[cm], n_blocks=n_blocks, n=n, ch=ch, chunk_length=chunk_length)