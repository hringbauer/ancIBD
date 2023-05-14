# run ancIBD given an input vcf file

import argparse
from ancIBD.run import hapBLOCK_chroms
from ancIBD.IO.prepare_h5 import vcf_to_1240K_hdf
from pathlib import Path
import h5py
import os


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run ancIBD.')
    parser.add_argument('--vcf', action="store", dest="vcf", type=str, required=True,
                        help="path to the imputed vcf file")
    parser.add_argument('--ch', action="store", dest="ch", type=int, required=True, help='chromosome number (1-22).')
    parser.add_argument('--marker_path', action="store", dest="marker_path", type=str, required=True, help='path to the marker file')
    parser.add_argument('--map_path', action="store", dest="map_path", type=str, required=True, help='path to the map file')
    parser.add_argument('--af_path', action="store", dest="af_path", type=str, required=False, help='path to the allele frequency file (optional)')
    parser.add_argument('--out', action="store", dest="out", type=str, required=False, help='output folder to store IBD results and the intermediary .hdf5 file. If not specified, the results will be stored in the same folder as the input vcf file.')
    parser.add_argument('--prefix', action="store", dest="prefix", type=str, required=False,
                        help="prefix of output file. If not specified, the prefix will be the same as the input vcf")
    parser.add_argument('--min', action="store", dest="min", type=float, required=False, default=8, help="minimum length of IBD segment in cM. Default is 8.")
    parser.add_argument('--iid', action="store", dest="iid", type=str, required=False, help="A list of sample iids to run ancIBD on (each line contains one sample IID). The sample list must match the sample name in the provided vcf file. If unspecified, ancIBD will run on all samples in the vcf file")
    parser.add_argument('--pair', action="store", dest="pair", type=str, required=False, help="A list of sample pairs to run ancIBD on (each line contains two sample IIDs separated by a whitespace). The sample list must match the sample name in the provided vcf file, and, if --iid is specified, all samples must also appear in the iid file. If unspecified, ancIBD will run on all pairs of samples in the vcf file")
    args = parser.parse_args()

    ch = args.ch
    path2vcf = args.vcf

    # determine prefix of output file
    prefix = args.prefix
    if prefix is None:
        prefix = Path(path2vcf).stem
        if prefix.endswith('.vcf'):
            prefix = prefix[:-4]

    # determine output folder
    oDir = args.out
    if oDir is None:
        oDir = Path(path2vcf).parent.absolute()
    else:
        # create the output folder if it does not already exist
        if os.path.isdir(oDir) is False:
            # create the output folder and its parental folders if they do not exist
            Path(oDir).mkdir(parents=True)

    path_vcf_1240k = os.path.join(f"{oDir}", f"{prefix}.1240k.vcf")
    path_h5 = os.path.join(f"{oDir}", f"{prefix}.ch{ch}.h5")
    vcf_to_1240K_hdf(in_vcf_path = path2vcf,
                 path_vcf = path_vcf_1240k, path_h5 = path_h5,
                 marker_path = f"/mnt/archgen/users/yilei/bin/ancIBD_data/filters/snps_bcftools_ch{ch}.csv",
                 map_path = "/mnt/archgen/users/yilei/bin/ancIBD_data/afs/v51.1_1240k.snp",
                 af_path = f"/mnt/archgen/users/yilei/bin/ancIBD_data/afs/v51.1_1240k_AF_ch{ch}.tsv",
                 col_sample_af = "",
                 buffer_size=20000, chunk_width=8, chunk_length=20000,
                 ch=ch)

    if args.iid is None:
        with h5py.File(path_h5, 'r') as f:
            iids = f['samples'][:].astype('str')
    else:
        iids = []
        with open(args.iid, 'r') as f:
            for line in f:
                iids.append(line.strip())
    
    run_iids = []
    if not args.pair is None:
        with open(args.pair, 'r') as f:
            for line in f:
                id1, id2 = line.strip().split()
                run_iids.append((id1, id2))
    
    # I set folder_out to an empty string because I save the output file separately in the code below
    df_ibd = hapBLOCK_chroms(folder_in=os.path.join(f"{oDir}", f"{prefix}.ch"),
                             iids=iids, run_iids=run_iids,
                             ch=ch, folder_out="",
                             output=False, prefix_out='', logfile=False,
                             l_model='hdf5', e_model='haploid_gl', h_model='FiveStateScaled', t_model='standard',
                             ibd_in=1, ibd_out=10, ibd_jump=400,
                             min_cm=args.min, cutoff_post=0.99, max_gap=0.0075)
    
    df_ibd.to_csv(os.path.join(f"{oDir}", f"{prefix}.ch{ch}.tsv"), sep='\t', index=False)
    os.remove(path_vcf_1240k) # remove intermediary vcf files
