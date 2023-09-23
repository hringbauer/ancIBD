## a wrapper script to summarize the results of ancIBD (combine results from all 22 autosomes)
from ancIBD.IO.ind_ibd import combine_all_chroms
from ancIBD.IO.ind_ibd import create_ind_ibd_df, create_ind_ibd_df_IBD2
import argparse
import os
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='Run ancIBD.')
    parser.add_argument('--tsv', action="store", dest="tsv", type=str, required=True,
                        help="base path to the individual IBD files.")
    parser.add_argument('--ch', action="store", dest="ch", type=str, required=False, default='1-22', help='chromosome number, expressed in the format chrom-chrom, e.g, 1-22). The default is 1-22.')
    parser.add_argument('--bin', action="store", dest="bin", type=str, required=False, default='8,12,16,20', help='length bin over which IBD sharing summary statistics for pairs of samples will be calculated. Default is 8,12,16,20.')
    parser.add_argument('--snp_cm', action="store", dest="snp_cm", type=float, required=False, default=220, help='minimum number of SNPs per centimorgan for a segment to be considered. The default is 220 to reduce false positive rates.')
    parser.add_argument('--out', action="store", dest="out", type=str, required=False, help='output folder to store results. If not specified, the results will be stored in the current directory.')
    parser.add_argument('--IBD2', action='store_true', dest='IBD2', default=False, help='if specified, the script will also summarize IBD2 segments. The default is False.')
    args = parser.parse_args()

    ch1, ch2 = args.ch.split('-')
    chs = range(int(ch1), int(ch2)+1)
    oDir = args.out
    if oDir is None:
        # set output directory to current directory
        oDir = os.getcwd()

    combine_all_chroms(chs=chs, folder_base=args.tsv,
                   path_save=os.path.join(oDir, 'ch_all.tsv'))
    
    #### now save a summary of IBD sharing between all pairs of individuals
    min_cms = [float(l) for l in args.bin.split(',')]
    if not args.IBD2:
        create_ind_ibd_df(ibd_data = os.path.join(oDir, 'ch_all.tsv'),
                      min_cms = min_cms, snp_cm = args.snp_cm, min_cm = min_cms[0], sort_col = 0,
                      savepath = os.path.join(oDir, "ibd_ind.tsv"))
    else:
        create_ind_ibd_df_IBD2(ibd_data = os.path.join(oDir, 'ch_all.tsv'),
                      min_cms = min_cms, snp_cm = args.snp_cm, min_cm = min_cms[0], sort_col = 0,
                      savepath = os.path.join(oDir, "ibd_ind.tsv"))