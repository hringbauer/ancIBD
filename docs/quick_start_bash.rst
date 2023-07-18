Running ancIBD via bash (available in the next release)
======================

Users not familiar with using Python functions can run ``ancIBD`` via bash, we have created two command-line interfaces (``ancIBD-run`` and ``ancIBD-summary``). Here we describe how to run those. 


Calling IBD
~~~~~~~~~~~
The test data used to run these tutorials can be downloaded from https://www.dropbox.com/sh/q18yyrffbdj1yv1/AAC1apifYB_oKB8SNrmQQ-26a?dl=0.

In addition to the imputed .vcf files, you need additionally three files,
all of which are provided in the same Dropbox link as indicated above.

-  marker_path: Path of the 1240k SNPs to use (you can find those in
   ``./filters/snps_bcftools_ch*.csv`` from the download link)
-  map_path: Path of the map file to use (eigenstrat .snp file, you can
   find it in ``./afs/v51.1_1240k.snp`` from the download link)
-  af_path (optional): Path of allele frequencies to merge into hdf5
   file (you can find it in ``./afs/v51.1_1240k_AF_ch*.tsv`` from the
   download link. If not provided, allele frequencies calculated from
   samples themselves will be used)


We showcase how to run ``ancIBD`` on ch20 as an example. The file path in the following tutorial assumes that the folder downloaded from the Dropbox link is in the same directory as this jupyter notebook. You will have to change the path to the above three files if needed.


.. code:: bash

    # Modify file paths according to your own environment if needed
    ancIBD-run --vcf ./data/vcf.raw/example_hazelton_chr20.vcf.gz --ch 20 --out test --marker_path ./data/filters/snps_bcftools_ch20.csv --map_path ./data/afs/v51.1_1240k.snp --af_path ./data/afs/v51.1_1240k_AF_ch20.tsv --prefix example_hazelton

.. literalinclude:: qs_output1.txt
   :language: console


If you already have the appropriate hdf5 file for your samples, you can also supply the command line with the hdf5 file directly. Make sure that the hdf5 file has the suffix “ch{chromosome number}.h5” (e.g. “test.ch20.h5”).

.. code:: bash

    ancIBD-run --h5 ./test/example_hazelton.ch20.h5 --ch 20 --out test --marker_path ./data/filters/snps_bcftools_ch20.csv --map_path ./data/afs/v51.1_1240k.snp --af_path ./data/afs/v51.1_1240k_AF_ch20.tsv --prefix example_hazelton

now we can do the same for all the 22 autosomes. This takes about 6min.

.. code:: bash
    
    for ch in {1..22};
    do
        marker_path=data/filters/snps_bcftools_ch$ch.csv
        af_path=data/afs/v51.1_1240k_AF_ch$ch.tsv
        vcf_path=data/vcf.raw/example_hazelton_chr$ch.vcf.gz
        ancIBD-run --vcf $vcf_path \
            --ch $ch --out test --marker_path $marker_path \
            --map_path ./data/afs/v51.1_1240k.snp \
            --af_path $af_path \
            --prefix example_hazelton
    done


.. note::


   For large sample sizes, we recommend that one parallizes over
   autosomes for speed-up (e.g., by submitting array jobs on a cluster).
   The above for-loop is efficient only for small sample sizes.

Combine IBD over 22 autosomes and generate summary statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we have individual IBD files for each of the autosomes, we can combine the information across chromosomes and obtain genome-wide summary statistics for all pairs of samples. Note that only pairs of samples that share at least one IBD passing the length cutoff are recorded).

.. code:: bash

    ancIBD-summary --tsv test/example_hazelton.ch --out test/

.. literalinclude:: qs_output2.txt
   :language: console


To view the complete options provided by the two command-line interfaces, you can use -h. 

.. code:: bash

    ancIBD-run -h


.. parsed-literal::

    usage: ancIBD-run [-h] [--vcf VCF] [--h5 H5] --ch CH --marker_path MARKER_PATH
                      --map_path MAP_PATH [--af_path AF_PATH] [--out OUT]
                      [--prefix PREFIX] [--min MIN] [--iid IID] [--pair PAIR]
    
    Run ancIBD.
    
    optional arguments:
      -h, --help            show this help message and exit
      --vcf VCF             path to the imputed vcf file
      --h5 H5               path to hdf5 file. If specified, ancIBD will skip the
                            vcf to hdf5 conversion step. Only one of --vcf and
                            --h5 should be specified.
      --ch CH               chromosome number (1-22).
      --marker_path MARKER_PATH
                            path to the marker file
      --map_path MAP_PATH   path to the map file
      --af_path AF_PATH     path to the allele frequency file (optional)
      --out OUT             output folder to store IBD results and the
                            intermediary .hdf5 file. If not specified, the results
                            will be stored in the same folder as the input vcf
                            file.
      --prefix PREFIX       prefix of output file. If not specified, the prefix
                            will be the same as the input vcf
      --min MIN             minimum length of IBD segment in cM. Default is 8.
      --iid IID             A list of sample iids to run ancIBD on (each line
                            contains one sample IID). The sample list must match
                            the sample name in the provided vcf file. If
                            unspecified, ancIBD will run on all samples in the vcf
                            file
      --pair PAIR           A list of sample pairs to run ancIBD on (each line
                            contains two sample IIDs separated by a whitespace).
                            The sample list must match the sample name in the
                            provided vcf file, and, if --iid is specified, all
                            samples must also appear in the iid file. If
                            unspecified, ancIBD will run on all pairs of samples
                            in the vcf file


.. code:: bash

    ancIBD-summary -h


.. parsed-literal::

    usage: ancIBD-summary [-h] --tsv TSV [--ch CH] [--bin BIN] [--snp_cm SNP_CM]
                          [--out OUT]
    
    Run ancIBD.
    
    optional arguments:
      -h, --help       show this help message and exit
      --tsv TSV        base path to the individual IBD files.
      --ch CH          chromosome number, expressed in the format chrom-chrom,
                       e.g, 1-22). The default is 1-22.
      --bin BIN        length bin over which IBD sharing summary statistics for
                       pairs of samples will be calculated. Default is 8,12,16,20.
      --snp_cm SNP_CM  minimum number of SNPs per centimorgan for a segment to be
                       considered. The default is 220 to reduce false positive
                       rates.
      --out OUT        output folder to store results. If not specified, the
                       results will be stored in the current directory.

