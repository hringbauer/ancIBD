import numpy as np
import pandas as pd
import os as os
import sys as sys

sys.path.insert(0, "/n/groups/reich/hringbauer/git/hapBLOCK/python3/prepare")  # hack to get local package first in path
from prepare_h5 import vcf_to_1240K_hdf

#####################################################
### Set Version. Manually make sure that Folders exist (some need creation)
vrs = "53.1"
v0 = vrs.split(".")[0]

#####################################################
### Run the Main Program
    
if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Script needs argument (which chromosome to run)")
    
    ### Load Input Parameters
    ch = int(sys.argv[1])  # Get Parameter of python function
    
    ## Run the Script
    ### Set Script Parameters
    base_path = f"/n/groups/reich/hringbauer/git/hapBLOCK"    
    in_vcf_path = f"/n/data1/hms/genetics/reich/1000Genomes/ali/WholeGenomeImputation/imputed_r2/v{vrs}/chr{ch}.bcf"
    path_vcf = f"{base_path}/data/vcf/1240k_v{vrs}/ch{ch}.vcf.gz"
    path_h5 = f"{base_path}/data/hdf5/1240k_v{vrs}/ch{ch}.h5"
    marker_path = f"{base_path}/data/filters/1240K_1000G/snps_bcftools_ch{ch}.csv"
    map_path = f"/n/groups/reich/DAVID/V{v0}/V{vrs}/v{vrs}_1240k_all.snp"
    
    ### Actual Run.
    vcf_to_1240K_hdf(in_vcf_path = in_vcf_path,
                     path_vcf = path_vcf,
                     path_h5 = path_h5,
                     marker_path = marker_path,
                     map_path = map_path,
                     buffer_size=20000, chunk_width=8, chunk_length=20000,
                     ch=ch)

    print(f"Finished processing chromosome: {ch}")