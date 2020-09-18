import numpy as np
import pandas as pd
import os as os
import sys as sys
import h5py
import allel

sys.path.insert(0, "/n/groups/reich/hringbauer/git/HAPSBURG/package/")  # hack to get local package first in path
from hapsburg.PackagesSupport.h5_python.h5_functions import merge_in_ld_map

#####################################################

def bctools_filter_vcf(in_vcf_path="", out_vcf_path="", marker_path=""):
    """Same as PLINK, but with bcftools and directly via Marker Positions.
    filter_iids: Whether to use the .csv with Indivdiduals to extract"""
    print(f"Running: bcftools view -Oz -o {out_vcf_path} -T {marker_path} -m2 -M2 -v snps {in_vcf_path}")
    os.system(f"bcftools view -Oz -o {out_vcf_path} -T {marker_path} -m2 -M2 -v snps {in_vcf_path}")
    print("Finished BCF tools filtering.")
    
def vcf_to_1240K_hdf(in_vcf_path = "", path_vcf = "",
                     path_h5 = "", marker_path="", ch=3,
                     map_path = "/n/groups/reich/DAVID/V43/V43.5/v43.5.snp"):
    """Convert Ali's vcf to 1240K hdf5"""    
    bctools_filter_vcf(in_vcf_path = in_vcf_path,
                       out_vcf_path= path_vcf,
                       marker_path = marker_path)
    print("Finished downsampling to 1240K")
    
    allel.vcf_to_hdf5(input=path_vcf, output=path_h5, 
                  fields = ['variants/*', 'calldata/*', "samples"], compression="gzip") # Do the conversion to hdf5. Takes hours
    print("Finished conversion to hdf5!")
    
    merge_in_ld_map(path_h5=path_h5, 
                path_snp1240k=map_path,
                chs=[ch])
    
    
if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError("Script needs argument (which chromosome to run)")
    
    ### Load correct set of IIDs
    ch = int(sys.argv[1])  # Get Parameter of python function
    
    in_vcf_path = f"/n/groups/reich/ali/WholeGenomeImputation/imputed/v43.4/chr{ch}.bcf"
    path_vcf = f"/n/groups/reich/hringbauer/git/hapBLOCK/data/vcf/1240k_v43/ch{ch}.vcf.gz"
    path_h5 = f"/n/groups/reich/hringbauer/git/hapBLOCK/data/hdf5/1240k_v43/ch{ch}.h5"
    marker_path = f"/n/groups/reich/hringbauer/git/hapBLOCK/data/filters/ho_snps_bcftools_ch{ch}.csv"
    
    vcf_to_1240K_hdf(in_vcf_path=in_vcf_path, path_vcf=path_vcf,
                     marker_path=marker_path, path_h5=path_h5, ch=ch)