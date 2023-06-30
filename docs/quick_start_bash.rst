Quick Start for ancIBD
========================

This notebook is a quick start guide for running ancIBD. It uses wrapper scripts for various functions introduced in section :ref:`preparing input <create_hdf5_from_vcf>`
and :ref:`calling IBD with ancIBD <run_ancIBD>`. 
Writing your own wrapper script for these functions provides more flexibility, while using the command line interface 
to be introduced in this quick starting guide is easier. 
We have created two command-line interfaces (``ancIBD-run`` and ``ancIBD-summary``)for running ancIBD quickly on your imputed data. 
The test data used to run these tutorials can be downloaded from https://www.dropbox.com/sh/q18yyrffbdj1yv1/AAC1apifYB_oKB8SNrmQQ-26a?dl=0. 
It contains imputed vcf of a subset of samples from early Neolithic Britain that belong to an extended pedigree 
(`Fowler et al. <https://www.nature.com/articles/s41586-021-04241-4>`_). 


calling IBD
***************

In addition to the imputed vcf files, you need additionally three files, all of which are provided in the same dropbox link as indicated above.

* marker_path: Path of the 1240k SNPs to use (you can find those in `./filters/snps_bcftools_ch*.csv` from the download link)
* map_path: Path of the map file to use (eigenstrat .snp file, you can find it in `./afs/v51.1_1240k.snp` from the download link)
* af_path (optional): Path of allele frequencies to merge into hdf5 file (you can find it in `./afs/v51.1_1240k_AF_ch*.tsv` from the download link. If not provided, allele frequencies calculated from samples themselves will be used)

We now run ancIBD on ch20 as an example. To run the following command, change the path to the above three files according to your own environment if needed. The file path in the following tutorial has assumed that the folder downloaded from dropbox link is in the same directory as this jupyter notebook.


