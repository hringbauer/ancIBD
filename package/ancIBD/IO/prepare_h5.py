"""
Functions to prepare HDF5 file from imputed VCFs
@ Author: Harald Ringbauer, 2021
"""

import numpy as np
import pandas as pd
import socket as socket
import os as os
import sys as sys
import h5py
import allel
from ancIBD.IO.h5_modify import merge_in_af, get_af, lift_af_df, merge_in_ld_map # get_af1000G, lift_af
#sys.path.insert(0, "/n/groups/reich/hringbauer/git/HAPSBURG/package/")  # hack to get local package first in path
#from hapsburg.PackagesSupport.h5_python.h5_functions import merge_in_ld_map
#sys.path.append("/n/groups/reich/hringbauer/git/hapBLOCK/python3/")

def save_1240kmarkers(snp1240k_path="", marker_path="", ch=0):
    """Save all 1240 Markers of .snp eigenstrat file.
    to marker_path.
    ch: Chromosome. If null filter all of them"""
    df_snp = pd.read_csv(snp1240k_path, header=None, 
                         sep=r"\s+", engine="python")
    df_snp.columns = ["SNP", "chr", "map", 
                      "pos", "ref", "alt"]
    if ch>0:
        df_snp = df_snp[df_snp["chr"] == ch]
    print(f"Loaded {len(df_snp)} Chr.{ch} SNPs.")

    df_save = df_snp[["chr", "pos"]]
    df_save.to_csv(marker_path, sep="\t", header=None, index=False)
    print(f"Saved {len(df_save)} 1240k Markers on Chr. {ch} to {marker_path}")
    
def save_1240_1000g_kmarkers(ch=3, snp_path="", marker_path=""):
    """Save all 1240 and 1000G Markers of .snp eigenstrat file.
    to marker_path. Loads Ali Path file
    snp_path: Where to find the SNPs plus their types"""
    df1=pd.read_csv(snp_path, sep="\t", header=None)
    df1.columns=["chr", "pos", "pos1", "ref", "alt", "type"]
    print(f"Loaded {len(df1)} SNPs on Chr. {ch}")
    df2 = df1[df1["type"]=="1kg+1240k"]
    print(f"Loaded {len(df2)} Chr.{ch} SNPs in both 1240K and 1000G.")
    df_save = df2[["chr", "pos1"]]
    df_save.to_csv(marker_path, sep="\t", header=None, index=False)
    print(f"Saved {len(df_save)} 1240k+1000G Markers on Chr. {ch} to {marker_path}")
    
def bctools_filter_vcf(in_vcf_path="", out_vcf_path="", marker_path=""):
    """Same as PLINK, but with bcftools and directly via Marker Positions.
    filter_iids: Whether to use the .csv with Indivdiduals to extract"""
    command = f"bcftools view -Oz -o {out_vcf_path} -T {marker_path} -M2 -v snps {in_vcf_path}"
    os.system(command)
    
    #!bcftools view -Oz -o $out_vcf_path -T $marker_path -M2 -v snps $in_vcf_path
    print("Finished BCF tools filtering.")
    
def bctools_filter_vcf_allvariants(in_vcf_path="", out_vcf_path="", marker_path=""):
    """Same as PLINK, but with bcftools and directly via Marker Positions.
    filter_iids: Whether to use the .csv with Indivdiduals to extract"""
    command = f"bcftools view -Oz -o {out_vcf_path} -T {marker_path} -v snps {in_vcf_path}"
    os.system(command)
    print("Finished BCF tools filtering.")
    
def merge_vcfs(in_vcf_paths=[], out_vcf_path=""):
    """Merges Set of VCFs into one VCF. 
    in_vcf_paths: List of VCFs to merge
    out_vcf_path: Output of VCF"""
    paths_merge = " ".join(in_vcf_paths)
    command = f"bcftools concat -n -o {out_vcf_path} {paths_merge}"
    os.system(command)
    print("Finished BCF tools filtering.")
    
##############################################################
### Function Identical to vcf_to_hdf5.py in cluster/ folder

def vcf_to_1240K_hdf(in_vcf_path = "/n/groups/reich/ali/WholeGenomeImputation/imputed/v43.4/chr3.bcf",
                     path_vcf = "./data/vcf/1240k_v43/ch3.vcf.gz",
                     path_h5 = "./data/hdf5/1240k_v43/ch3.h5",
                     marker_path="./data/filters/ho_snps_bcftools_ch3.csv",
                     map_path="/n/groups/reich/DAVID/V43/V43.5/v43.5.snp",
                     af_path="",
                     col_sample_af="AF_SAMPLE",
                     chunk_length=10000, chunk_width=8, buffer_size=20000,
                     ch=3):
    """Convert Ali's vcf to 1240K hdf5. 
    marker_path: Path to file containing SNPs to downsample to
    If marker_path empty, no SNP filtering done.
    map_path: Path to eigenstrat SNP file containing genetic map.
    These are merged into a hdf5 field "variants/MAP"
    If map_path empty, no genetic map is merged in.
    af_path: Path to tab-seperated table containing allele frequencies
    There are merged into the hdf5 field "variants/AF_ALL"
    If no such path given, no allele frequencies are merged in.
    col_sample_af: The hdf5 column name for the allele frequency calculated from sample.
    If left empty, no such column will be calculated or added. 
    """ 
    print("Print downsampling to 1240K...")
    if len(marker_path)>0:
        bctools_filter_vcf(in_vcf_path = in_vcf_path,
                           out_vcf_path= path_vcf,
                           marker_path = marker_path)
    else:
        path_vcf = in_vcf_path  ### If no temporary VCF available
    
    if os.path.exists(path_h5):  # Do a Deletion of existing hdf5 File there
        print(f"Deleting previous HDF5 file at path_h5: {path_h5}...")
        os.remove(path_h5)
    
    print("Converting to HDF5...")
    allel.vcf_to_hdf5(input=path_vcf, output=path_h5, 
                      fields = ['variants/*', 'calldata/*', "samples"], 
                      types = {"samples":"S60", "calldata/GT":np.int8,
                               "calldata/GP":np.float32, "calldata/PL":np.float32}, 
                      buffer_size=buffer_size,
                      chunk_length = chunk_length, chunk_width=chunk_width,
                      compression="gzip") # Do the conversion to hdf5. Takes hours
    print("Finished conversion to hdf5!")
    
    print("Merging in LD Map..")
    if len(map_path)>0:
        merge_in_ld_map(path_h5=path_h5, 
                    path_snp1240k=map_path,
                    chs=[ch])
        
    ### Calculate and Merge in Allele Frequency
    if len(col_sample_af)>0:
        print(f"Calculating in sample allele frequencies and saving at hdf5 column {col_sample_af}")
        with h5py.File(path_h5, "r") as f:
            af = get_af(f) 
        merge_in_af(path_h5, af, col_af=col_sample_af)
        
    if len(af_path)>0:
        lift_af_df(h5_target=path_h5, path_df=af_path,
                   field='variants/AF_ALL') # The field to use for Allele Frequencies

    print(f"Transformation complete! Find new hdf5 file at: {path_h5}\n")