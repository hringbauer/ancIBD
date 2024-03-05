# run ancIBDX given an input vcf file

import argparse
from ancIBD.run import hapBLOCK_chroms_mixedPloidy
from ancIBD.IO.prepare_h5 import vcf_to_1240K_hdf
from ancIBD.IO.ind_ibd import create_ind_ibd_df
from pathlib import Path
import h5py
import os
import warnings

def main():
    parser = argparse.ArgumentParser(description='Run ancIBD (on Autosomes).')
    parser.add_argument('--vcf', action="store", dest="vcf", type=str, required=False,
                        help="path to the imputed vcf file")
    parser.add_argument('--h5', action="store", dest="h5", type=str, required=False,
                        help="path to hdf5 file. If specified, ancIBD will skip the vcf to hdf5 conversion step. \
                        Only one of --vcf and --h5 should be specified.\
                        But please make sure that the hdf5 file has suffix ch{chromosome number}.h5 (e.g, test.ch20.h5).")
    parser.add_argument('--ch', action="store", dest="ch", type=str, required=False, default='X', help='chromosome name, default by X.')
    parser.add_argument('--marker_path', action="store", dest="marker_path", type=str, required=False, help='path to the marker file')
    parser.add_argument('--map_path', action="store", dest="map_path", type=str, required=False, help='path to the map file')
    parser.add_argument('--af_path', action="store", dest="af_path", type=str, required=False, default="", help='path to the allele frequency file (optional)')
    parser.add_argument('--af_column', action='store', dest='af_column', type=str, required=False, default='', help='column name of the allele frequency in the hdf5. For example, "variants/AF_ALL" or "variants/AF_SAMPLE".')
    parser.add_argument('--out', action="store", dest="out", type=str, required=False, help='output folder to store IBD results and the intermediary .hdf5 file. If not specified, the results will be stored in the same folder as the input vcf file.')
    parser.add_argument('--prefix', action="store", dest="prefix", type=str, required=False,
                        help="prefix of output file. If not specified, the prefix will be the same as the input vcf")
    parser.add_argument('--ibd-in', action="store", dest="ibd_in", type=float, required=False, default=1, help="IBD in rate. Default is 1.")
    parser.add_argument('--ibd-out', action="store", dest="ibd_out", type=float, required=False, default=10, help="IBD out rate. Default is 10.")
    parser.add_argument('--ibd-jump', action="store", dest="ibd_jump", type=float, required=False, default=400, help="IBD jump rate. Default is 400.")
    parser.add_argument('--min', action="store", dest="min", type=float, required=False, default=8, help="minimum length of IBD segment in cM. Default is 8.")
    parser.add_argument('--ploidy', action="store", dest="ploidy", type=str, required=True, default='1', 
                        help="path to a .txt file listing the ploidy of each sample. \
                        Each sample takes one line, the first column is sample iid and the second column is the ploidy. \
                        The two columns can be separted by any white space.\
                        ancIBD will only run on individuals listed in this file, even if others exist in the vcf/hdf5 file.")
    parser.add_argument('--pair', action="store", dest="pair", type=str, required=False, help="A list of sample pairs to run ancIBD on (each line contains two sample IIDs separated by a whitespace). The sample list must match the sample name in the provided vcf file, and, if --iid is specified, all samples must also appear in the iid file. If unspecified, ancIBD will run on all pairs of samples in the vcf file")
    parser.add_argument('--mask', action='store', dest='mask', type=str, required=False, default="", help='Mask file to mask out regions for IBD calling.')
    parser.add_argument('--bin', action="store", dest="bin", type=str, required=False, default='8,12,16,20', help='length bin over which IBD sharing summary statistics for pairs of samples will be calculated. Default is 8,12,16,20.')
    parser.add_argument('--snpcm', action="store", dest="snpcm", type=int, required=False, default=180, help='minimum number of SNPs per cM. Default is 180.')
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
            Path(oDir).mkdir(parents=True)

    if args.vcf is not None:
        path_vcf_1240k = os.path.join(f"{oDir}", f"{prefix}.1240k.vcf") # temporary vcf file
        path_h5 = os.path.join(f"{oDir}", f"{prefix}.h5")
        col_sample_af = "" if len(args.af_path) > 0 else "AF_SAMPLE"
        if len(args.af_path) == 0 and len(args.af_column) == 0:
            print("WARNING: allele frequency file or entries in VCF that stores AF not provided. \
                The allele frequency will be calculated from the input vcf file. \
                  Make sure you have enough sample size to have good estimates of allele frequencies ")
            print('The sample allele frequency will be stored in the "variants/AF_SAMPLE" column of the output hdf5 file.')
        vcf_to_1240K_hdf(in_vcf_path = args.vcf,
                 path_vcf = path_vcf_1240k, path_h5 = path_h5,
                 marker_path = args.marker_path,
                 map_path = args.map_path,
                 af_path = args.af_path,
                 col_sample_af = col_sample_af,
                 buffer_size=20000, chunk_width=8, chunk_length=20000,
                 ch=ch)
        os.remove(path_vcf_1240k) # remove intermediary vcf files
    else:
        path_h5 = args.h5

    with h5py.File(path_h5, 'r') as f:
        iids = f['samples'][:].astype('str')
    iids2run = []
    ploidy = []
    with open(args.ploidy, 'r') as f:
        for line in f:
            iid, p = line.strip().split()
            if iid in iids:
                iids2run.append(iid)
                ploidy.append(int(p))
            else:
                warn = f'{iid} in the ploidy file is not in the hdf5 file. It will be ignored...'
                warnings.warn(warn)
    assert(len(iids2run) == len(ploidy))
    # check if all samples in iids are in iids2run
    for iid in iids:
        if iid not in iids2run:
            warn = f'{iid} in the hdf5 file is not in the ploidy file. It will be ignored...'
            warnings.warn(warn)
    
    run_iids = []
    if not args.pair is None:
        with open(args.pair, 'r') as f:
            for line in f:
                id1, id2 = line.strip().split()
                if (not id1 in iids2run) or (not id2 in iids2run):
                    warn = f'{id1} or {id2} not in the ploidy/hdf5 file. It will be ignored...'
                    warnings.warn(warn)
                    continue
                run_iids.append((id1, id2))
    
    # I set folder_out to an empty string because I save the output file separately in the code below
    p_col = args.af_column
    if len(p_col) == 0:
        p_col = 'variants/AF_ALL' if len(args.af_path) > 0 else 'variants/AF_SAMPLE'
    df_ibd = hapBLOCK_chroms_mixedPloidy(folder_in=path_h5[:2+path_h5.rfind('ch')],
                             iids=iids2run, ploidy=ploidy, run_iids=run_iids,
                             ch=ch, folder_out="",
                             output=args.verbose, prefix_out='', logfile=False, p_col=p_col,
                             ibd_in=args.ibd_in, ibd_out=args.ibd_out, ibd_jump=args.ibd_jump,
                             min_cm=args.min, cutoff_post=0.99, max_gap=0.0075,
                             mask=args.mask)
    df_ibd.to_csv(os.path.join(f"{oDir}", f"{prefix}.tsv"), sep='\t', index=False)
    # create a IBD summary statistics file
    min_cms = [float(l) for l in args.bin.split(',')]
    create_ind_ibd_df(ibd_data = os.path.join(oDir, f"{prefix}.tsv"),
                      min_cms = min_cms, snp_cm = args.snpcm, min_cm = min_cms[0], sort_col = 0,
                      savepath = os.path.join(oDir, "ibd_ind.tsv"))
    
