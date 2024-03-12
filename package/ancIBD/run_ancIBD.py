# run ancIBD given an input vcf file

import argparse
from ancIBD.run import hapBLOCK_chroms
from ancIBD.IO.prepare_h5 import vcf_to_1240K_hdf
from pathlib import Path
import h5py
import os


def main():
    parser = argparse.ArgumentParser(description='Run ancIBD (on Autosomes).')
    parser.add_argument('--vcf', action="store", dest="vcf", type=str, required=False,
                        help="path to the imputed vcf file")
    parser.add_argument('--h5', action="store", dest="h5", type=str, required=False,
                        help="path to hdf5 file. If specified, ancIBD will skip the vcf to hdf5 conversion step. \
                        Only one of --vcf and --h5 should be specified.\
                        But please make sure that the hdf5 file has suffix ch{chromosome number}.h5 (e.g, test.ch20.h5).")
    parser.add_argument('--ch', action="store", dest="ch", type=int, required=True, help='chromosome number (1-22).')
    parser.add_argument('--marker_path', action="store", dest="marker_path", type=str, required=False, default="", help='path to the marker file')
    parser.add_argument('--map_path', action="store", dest="map_path", type=str, required=False, help='path to the map file')
    parser.add_argument('--af_path', action="store", dest="af_path", type=str, required=False, default="", help='path to the allele frequency file (optional)')
    parser.add_argument('--af_column', action='store', dest='af_column', type=str, required=False, default='', help='column name of the allele frequency in the hdf5. For example, "variants/AF_ALL" or "variants/AF_SAMPLE".')
    parser.add_argument('--out', action="store", dest="out", type=str, required=False, help='output folder to store IBD results and the intermediary .hdf5 file. If not specified, the results will be stored in the same folder as the input vcf file.')
    parser.add_argument('--prefix', action="store", dest="prefix", type=str, required=False,
                        help="prefix of output file. If not specified, the prefix will be the same as the input vcf")
    parser.add_argument('--ibd-in', action="store", dest="ibd_in", type=float, required=False, default=1, help="IBD in rate. Default is 1.")
    parser.add_argument('--ibd-out', action="store", dest="ibd_out", type=float, required=False, default=10, help="IBD out rate. Default is 10.")
    parser.add_argument('--min', action="store", dest="min", type=float, required=False, default=8, help="minimum length of IBD segment in cM. Default is 8.")
    parser.add_argument('--iid', action="store", dest="iid", type=str, required=False, help="A list of sample iids to run ancIBD on (each line contains one sample IID). The sample list must match the sample name in the provided vcf file. If unspecified, ancIBD will run on all samples in the vcf file")
    parser.add_argument('--pair', action="store", dest="pair", type=str, required=False, help="A list of sample pairs to run ancIBD on (each line contains two sample IIDs separated by a whitespace). The sample list must match the sample name in the provided vcf file, and, if --iid is specified, all samples must also appear in the iid file. If unspecified, ancIBD will run on all pairs of samples in the vcf file")
    parser.add_argument('--IBD2', action="store_true", dest="IBD2", required=False, help="If specified, ancIBD will enable the detection of IBD2 segments. Default is false.")
    parser.add_argument('--post2', action='store', dest='post2', type=float, required=False, default=0.8, help='posterior probability cutoff for IBD2 segments. Default is 0.975.')
    parser.add_argument('--min_cm2_init', action='store', dest='min_cm2_init', type=float, required=False, default=0.25, help='minimum length of IBD2 segment in cM before merging. Default is 1.0.')
    parser.add_argument('--min_cm2_after_merge', action='store', dest='min_cm2_after_merge', type=float, required=False, default=4.0, help='minimum length of IBD2 segment in cM after merging. Default is 2.0.')
    parser.add_argument('--mask', action='store', dest='mask', type=str, required=False, default="", help='Mask file to mask out regions for IBD calling.')
    parser.add_argument('-v', '--verbose', action='store_true', dest='verbose', required=False, help='turn on verbose mode')
    args = parser.parse_args()

    if args.vcf is None and args.h5 is None:
        raise ValueError("One of --vcf and --h5 must be specified.")
    elif args.vcf is not None and args.h5 is not None:
        raise ValueError("Only one of --vcf and --h5 should be specified.")

    ch = args.ch

    # determine prefix of output file
    prefix = args.prefix
    if prefix is None:
        prefix = f'ch{ch}'
    else:
        prefix += f'.ch{ch}'

    # determine output folder
    oDir = args.out
    if oDir is None:
        oDir = Path(args.vcf).parent.absolute() if args.vcf is not None else Path(args.h5).parent.absolute()
    else:
        # create the output folder if it does not already exist
        if os.path.isdir(oDir) is False:
            # create the output folder and its parental folders if they do not exist
            Path(oDir).mkdir(parents=True, exist_ok=True)

    if args.vcf is not None:
        path_vcf_1240k = os.path.join(f"{oDir}", f"{prefix}.1240k.vcf")
        path_h5 = os.path.join(f"{oDir}", f"{prefix}.h5")
        col_sample_af = "" if len(args.af_path) > 0 else "AF_SAMPLE"
        if len(args.af_path) == 0:
            print("WARNING: allele frequency file not provided. The allele frequency will be calculated from the input vcf file. Make sure you have enough sample size to have good estimates of allele frequencies ")
            print('The sample allele frequency will be stored in the "variants/AF_SAMPLE" column of the output hdf5 file.')
        vcf_to_1240K_hdf(in_vcf_path = args.vcf,
                 path_vcf = path_vcf_1240k, path_h5 = path_h5,
                 marker_path = args.marker_path,
                 map_path = args.map_path,
                 af_path = args.af_path,
                 col_sample_af = col_sample_af,
                 buffer_size=20000, chunk_width=8, chunk_length=20000,
                 ch=ch)
        if len(args.marker_path) > 0:
            os.remove(path_vcf_1240k) # remove intermediary vcf files
    else:
        path_h5 = args.h5

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
    p_col = args.af_column
    if len(p_col) == 0:
        p_col = 'variants/AF_ALL' if len(args.af_path) > 0 else 'variants/AF_SAMPLE'
    df_ibd = hapBLOCK_chroms(folder_in=path_h5[:2+path_h5.rfind('ch')],
                             iids=iids, run_iids=run_iids,
                             ch=ch, folder_out="",
                             output=args.verbose, prefix_out='', logfile=False,
                             l_model='h5', e_model='haploid_gl2', h_model='FiveStateScaled', t_model='standard', p_col=p_col,
                             ibd_in=args.ibd_in, ibd_out=args.ibd_out, ibd_jump=400,
                             min_cm=args.min, cutoff_post=0.99, max_gap=0.0075, 
                             IBD2=args.IBD2, cutoff_post2=args.post2, min_cm2_init=args.min_cm2_init, min_cm2_after_merge=args.min_cm2_after_merge,
                             mask=args.mask)    
    df_ibd.to_csv(os.path.join(f"{oDir}", f"{prefix}.tsv"), sep='\t', index=False)
    
