Quick Start for ancIBD
======================

This notebook is a quick start guide for running ancIBD. It uses wrapper
scripts for various functions introduced in section `preparing
input <create_hdf5_from_vcf.ipynb>`__ and `calling IBD with
ancIBD <run_ancIBD.ipynb>`__. Writing your own wrapper script for these
functions provides more flexibility, while using the command line
interface to be introduced in this quick starting guide is easier. We
have created two command-line interfaces (``ancIBD-run`` and
``ancIBD-summary``)for running ancIBD quickly on your imputed data. The
test data used to run these tutorials can be downloaded from
https://www.dropbox.com/sh/q18yyrffbdj1yv1/AAC1apifYB_oKB8SNrmQQ-26a?dl=0.
It contains imputed vcf of a subset of samples from early Neolithic
Britain that belong to an extended pedigree (`Fowler et
al. <https://www.nature.com/articles/s41586-021-04241-4>`__).

calling IBD
~~~~~~~~~~~

In addition to the imputed vcf files, you need additionally three files,
all of which are provided in the same dropbox link as indicated above.

-  marker_path: Path of the 1240k SNPs to use (you can find those in
   ``./filters/snps_bcftools_ch*.csv`` from the download link)
-  map_path: Path of the map file to use (eigenstrat .snp file, you can
   find it in ``./afs/v51.1_1240k.snp`` from the download link)
-  af_path (optional): Path of allele frequencies to merge into hdf5
   file (you can find it in ``./afs/v51.1_1240k_AF_ch*.tsv`` from the
   download link. If not provided, allele frequencies calculated from
   samples themselves will be used)

We now run ancIBD on ch20 as an example. To run the following command,
change the path to the above three files according to your own
environment if needed. The file path in the following tutorial has
assumed that the folder downloaded from dropbox link is in the same
directory as this jupyter notebook.

.. command-output:: python --version


.. code:: bash

    # Modify file paths according to your own environment if needed
    ancIBD-run --vcf ./data/vcf.raw/example_hazelton_chr20.vcf.gz --ch 20 --out test --marker_path ./data/filters/snps_bcftools_ch20.csv --map_path ./data/afs/v51.1_1240k.snp --af_path ./data/afs/v51.1_1240k_AF_ch20.tsv --prefix example_hazelton


.. parsed-literal::

    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch20.1240k.vcf -T ./data/filters/snps_bcftools_ch20.csv -M2 -v snps ./data/vcf.raw/example_hazelton_chr20.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 28940 variants.
    Loaded 6 individuals.
    Loaded 30377 Chr.20 1240K SNPs.
    Intersection 28827 out of 28940 HDF5 SNPs
    Interpolating 113 variants.
    Finished Chromosome 20.
    Adding map to HDF5...
    Intersection 28827 out of 28940 target HDF5 SNPs. 113 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch20.h5
    


If you already have the appropriate hdf5 file for your samples, you can
also supply the command line with the hdf5 file directly. But please
make sure that the hdf5 file has suffix “ch{chromosome number}.h5” (e.g,
“test.ch20.h5”).

.. code:: ipython3

    !ancIBD-run --h5 ./test/example_hazelton.ch20.h5 --ch $ch --out test --marker_path $marker_path --map_path $map_path --af_path $af_path --prefix example_hazelton

now we can do the same for the all the 22 autosomes. This takes about
6min.

.. code:: bash

    %%bash
    
    map_path='./data/afs/v51.1_1240k.snp'
    
    for ch in {1..22};
    do
        marker_path=data/filters/snps_bcftools_ch$ch.csv
        af_path=data/afs/v51.1_1240k_AF_ch$ch.tsv
        vcf_path=data/vcf.raw/example_hazelton_chr$ch.vcf.gz
        ancIBD-run --vcf $vcf_path \
            --ch $ch --out test --marker_path $marker_path --map_path $map_path --af_path $af_path --prefix example_hazelton
    done


.. parsed-literal::

    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch1.1240k.vcf -T data/filters/snps_bcftools_ch1.csv -M2 -v snps data/vcf.raw/example_hazelton_chr1.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 88408 variants.
    Loaded 6 individuals.
    Loaded 93166 Chr.1 1240K SNPs.
    Intersection 88115 out of 88408 HDF5 SNPs
    Interpolating 293 variants.
    Finished Chromosome 1.
    Adding map to HDF5...
    Intersection 88115 out of 88408 target HDF5 SNPs. 293 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch1.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch2.1240k.vcf -T data/filters/snps_bcftools_ch2.csv -M2 -v snps data/vcf.raw/example_hazelton_chr2.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 93875 variants.
    Loaded 6 individuals.
    Loaded 98657 Chr.2 1240K SNPs.
    Intersection 93471 out of 93875 HDF5 SNPs
    Interpolating 404 variants.
    Finished Chromosome 2.
    Adding map to HDF5...
    Intersection 93471 out of 93875 target HDF5 SNPs. 404 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch2.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch3.1240k.vcf -T data/filters/snps_bcftools_ch3.csv -M2 -v snps data/vcf.raw/example_hazelton_chr3.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 77345 variants.
    Loaded 6 individuals.
    Loaded 81416 Chr.3 1240K SNPs.
    Intersection 77013 out of 77345 HDF5 SNPs
    Interpolating 332 variants.
    Finished Chromosome 3.
    Adding map to HDF5...
    Intersection 77013 out of 77345 target HDF5 SNPs. 332 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch3.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch4.1240k.vcf -T data/filters/snps_bcftools_ch4.csv -M2 -v snps data/vcf.raw/example_hazelton_chr4.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 68518 variants.
    Loaded 6 individuals.
    Loaded 71634 Chr.4 1240K SNPs.
    Intersection 68254 out of 68518 HDF5 SNPs
    Interpolating 264 variants.
    Finished Chromosome 4.
    Adding map to HDF5...
    Intersection 68254 out of 68518 target HDF5 SNPs. 264 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch4.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch5.1240k.vcf -T data/filters/snps_bcftools_ch5.csv -M2 -v snps data/vcf.raw/example_hazelton_chr5.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 69063 variants.
    Loaded 6 individuals.
    Loaded 74004 Chr.5 1240K SNPs.
    Intersection 68899 out of 69063 HDF5 SNPs
    Interpolating 164 variants.
    Finished Chromosome 5.
    Adding map to HDF5...
    Intersection 68899 out of 69063 target HDF5 SNPs. 164 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch5.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch6.1240k.vcf -T data/filters/snps_bcftools_ch6.csv -M2 -v snps data/vcf.raw/example_hazelton_chr6.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 75347 variants.
    Loaded 6 individuals.
    Loaded 78867 Chr.6 1240K SNPs.
    Intersection 75059 out of 75347 HDF5 SNPs
    Interpolating 288 variants.
    Finished Chromosome 6.
    Adding map to HDF5...
    Intersection 75059 out of 75347 target HDF5 SNPs. 288 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch6.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch7.1240k.vcf -T data/filters/snps_bcftools_ch7.csv -M2 -v snps data/vcf.raw/example_hazelton_chr7.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 59603 variants.
    Loaded 6 individuals.
    Loaded 62595 Chr.7 1240K SNPs.
    Intersection 59324 out of 59603 HDF5 SNPs
    Interpolating 279 variants.
    Finished Chromosome 7.
    Adding map to HDF5...
    Intersection 59324 out of 59603 target HDF5 SNPs. 279 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch7.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch8.1240k.vcf -T data/filters/snps_bcftools_ch8.csv -M2 -v snps data/vcf.raw/example_hazelton_chr8.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 60828 variants.
    Loaded 6 individuals.
    Loaded 63916 Chr.8 1240K SNPs.
    Intersection 60530 out of 60828 HDF5 SNPs
    Interpolating 298 variants.
    Finished Chromosome 8.
    Adding map to HDF5...
    Intersection 60530 out of 60828 target HDF5 SNPs. 298 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch8.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch9.1240k.vcf -T data/filters/snps_bcftools_ch9.csv -M2 -v snps data/vcf.raw/example_hazelton_chr9.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 50546 variants.
    Loaded 6 individuals.
    Loaded 52765 Chr.9 1240K SNPs.
    Intersection 50307 out of 50546 HDF5 SNPs
    Interpolating 239 variants.
    Finished Chromosome 9.
    Adding map to HDF5...
    Intersection 50307 out of 50546 target HDF5 SNPs. 239 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch9.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch10.1240k.vcf -T data/filters/snps_bcftools_ch10.csv -M2 -v snps data/vcf.raw/example_hazelton_chr10.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 58610 variants.
    Loaded 6 individuals.
    Loaded 61131 Chr.10 1240K SNPs.
    Intersection 58364 out of 58610 HDF5 SNPs
    Interpolating 246 variants.
    Finished Chromosome 10.
    Adding map to HDF5...
    Intersection 58364 out of 58610 target HDF5 SNPs. 246 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch10.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch11.1240k.vcf -T data/filters/snps_bcftools_ch11.csv -M2 -v snps data/vcf.raw/example_hazelton_chr11.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 54590 variants.
    Loaded 6 individuals.
    Loaded 57163 Chr.11 1240K SNPs.
    Intersection 54365 out of 54590 HDF5 SNPs
    Interpolating 225 variants.
    Finished Chromosome 11.
    Adding map to HDF5...
    Intersection 54365 out of 54590 target HDF5 SNPs. 225 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch11.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch12.1240k.vcf -T data/filters/snps_bcftools_ch12.csv -M2 -v snps data/vcf.raw/example_hazelton_chr12.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 53737 variants.
    Loaded 6 individuals.
    Loaded 56133 Chr.12 1240K SNPs.
    Intersection 53528 out of 53737 HDF5 SNPs
    Interpolating 209 variants.
    Finished Chromosome 12.
    Adding map to HDF5...
    Intersection 53528 out of 53737 target HDF5 SNPs. 209 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch12.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch13.1240k.vcf -T data/filters/snps_bcftools_ch13.csv -M2 -v snps data/vcf.raw/example_hazelton_chr13.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 38927 variants.
    Loaded 6 individuals.
    Loaded 40441 Chr.13 1240K SNPs.
    Intersection 38774 out of 38927 HDF5 SNPs
    Interpolating 153 variants.
    Finished Chromosome 13.
    Adding map to HDF5...
    Intersection 38774 out of 38927 target HDF5 SNPs. 153 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch13.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch14.1240k.vcf -T data/filters/snps_bcftools_ch14.csv -M2 -v snps data/vcf.raw/example_hazelton_chr14.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 35885 variants.
    Loaded 6 individuals.
    Loaded 37903 Chr.14 1240K SNPs.
    Intersection 35744 out of 35885 HDF5 SNPs
    Interpolating 141 variants.
    Finished Chromosome 14.
    Adding map to HDF5...
    Intersection 35744 out of 35885 target HDF5 SNPs. 141 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch14.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch15.1240k.vcf -T data/filters/snps_bcftools_ch15.csv -M2 -v snps data/vcf.raw/example_hazelton_chr15.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 34280 variants.
    Loaded 6 individuals.
    Loaded 35991 Chr.15 1240K SNPs.
    Intersection 34159 out of 34280 HDF5 SNPs
    Interpolating 121 variants.
    Finished Chromosome 15.
    Adding map to HDF5...
    Intersection 34159 out of 34280 target HDF5 SNPs. 121 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch15.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch16.1240k.vcf -T data/filters/snps_bcftools_ch16.csv -M2 -v snps data/vcf.raw/example_hazelton_chr16.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 34335 variants.
    Loaded 6 individuals.
    Loaded 36000 Chr.16 1240K SNPs.
    Intersection 34138 out of 34335 HDF5 SNPs
    Interpolating 198 variants.
    Finished Chromosome 16.
    Adding map to HDF5...
    Intersection 34138 out of 34335 target HDF5 SNPs. 197 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch16.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch17.1240k.vcf -T data/filters/snps_bcftools_ch17.csv -M2 -v snps data/vcf.raw/example_hazelton_chr17.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 28892 variants.
    Loaded 6 individuals.
    Loaded 30733 Chr.17 1240K SNPs.
    Intersection 28794 out of 28892 HDF5 SNPs
    Interpolating 98 variants.
    Finished Chromosome 17.
    Adding map to HDF5...
    Intersection 28794 out of 28892 target HDF5 SNPs. 98 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch17.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch18.1240k.vcf -T data/filters/snps_bcftools_ch18.csv -M2 -v snps data/vcf.raw/example_hazelton_chr18.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 33846 variants.
    Loaded 6 individuals.
    Loaded 35327 Chr.18 1240K SNPs.
    Intersection 33720 out of 33846 HDF5 SNPs
    Interpolating 126 variants.
    Finished Chromosome 18.
    Adding map to HDF5...
    Intersection 33720 out of 33846 target HDF5 SNPs. 126 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch18.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch19.1240k.vcf -T data/filters/snps_bcftools_ch19.csv -M2 -v snps data/vcf.raw/example_hazelton_chr19.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 18092 variants.
    Loaded 6 individuals.
    Loaded 19273 Chr.19 1240K SNPs.
    Intersection 18018 out of 18092 HDF5 SNPs
    Interpolating 74 variants.
    Finished Chromosome 19.
    Adding map to HDF5...
    Intersection 18018 out of 18092 target HDF5 SNPs. 74 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch19.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch20.1240k.vcf -T data/filters/snps_bcftools_ch20.csv -M2 -v snps data/vcf.raw/example_hazelton_chr20.vcf.gz
    Finished BCF tools filtering to target markers.
    Deleting previous HDF5 file at path_h5: test/example_hazelton.ch20.h5...
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 28940 variants.
    Loaded 6 individuals.
    Loaded 30377 Chr.20 1240K SNPs.
    Intersection 28827 out of 28940 HDF5 SNPs
    Interpolating 113 variants.
    Finished Chromosome 20.
    Adding map to HDF5...
    Intersection 28827 out of 28940 target HDF5 SNPs. 113 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch20.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch21.1240k.vcf -T data/filters/snps_bcftools_ch21.csv -M2 -v snps data/vcf.raw/example_hazelton_chr21.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 15707 variants.
    Loaded 6 individuals.
    Loaded 16727 Chr.21 1240K SNPs.
    Intersection 15640 out of 15707 HDF5 SNPs
    Interpolating 67 variants.
    Finished Chromosome 21.
    Adding map to HDF5...
    Intersection 15640 out of 15707 target HDF5 SNPs. 67 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch21.h5
    
    Print downsampling to 1240K...
    Running bash command: 
    bcftools view -Ov -o test/example_hazelton.ch22.1240k.vcf -T data/filters/snps_bcftools_ch22.csv -M2 -v snps data/vcf.raw/example_hazelton_chr22.vcf.gz
    Finished BCF tools filtering to target markers.
    Converting to HDF5...
    Finished conversion to hdf5!
    Merging in LD Map..
    Lifting LD Map from eigenstrat to HDF5...
    Loaded 15483 variants.
    Loaded 6 individuals.
    Loaded 16420 Chr.22 1240K SNPs.
    Intersection 15408 out of 15483 HDF5 SNPs
    Interpolating 75 variants.
    Finished Chromosome 22.
    Adding map to HDF5...
    Intersection 15408 out of 15483 target HDF5 SNPs. 75 SNPs set to AF=0.5
    Transformation complete! Find new hdf5 file at: test/example_hazelton.ch22.h5
    


.. container:: alert alert-info

   Note

   For large sample sizes, we recommend that one parallizes over
   autosomes for speed-up (e.g, by submitting array jobs on a cluster).
   The above for-loop is efficient only for small sample sizes.

Combine IBD over 22 autosomes and generate summary statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now that we have individual IBD files for each of the autosome, we can
combine the information across chromosomes and obtain genome-wide
summary statistics for all pairs of samples (Only pairs of samples that
share at least one IBD passing the length cutoff are recorded).

.. code:: ipython3

    !ancIBD-summary --tsv test/example_hazelton.ch --out test/


.. parsed-literal::

    Chromosome 1; Loaded 10 IBD
    Chromosome 2; Loaded 9 IBD
    Chromosome 3; Loaded 6 IBD
    Chromosome 4; Loaded 9 IBD
    Chromosome 5; Loaded 8 IBD
    Chromosome 6; Loaded 7 IBD
    Chromosome 7; Loaded 9 IBD
    Chromosome 8; Loaded 7 IBD
    Chromosome 9; Loaded 6 IBD
    Chromosome 10; Loaded 7 IBD
    Chromosome 11; Loaded 5 IBD
    Chromosome 12; Loaded 5 IBD
    Chromosome 13; Loaded 8 IBD
    Chromosome 14; Loaded 6 IBD
    Chromosome 15; Loaded 3 IBD
    Chromosome 16; Loaded 6 IBD
    Chromosome 17; Loaded 4 IBD
    Chromosome 18; Loaded 5 IBD
    Chromosome 19; Loaded 8 IBD
    Chromosome 20; Loaded 6 IBD
    Chromosome 21; Loaded 6 IBD
    Chromosome 22; Loaded 6 IBD
    Saved 146 IBD to test/ch_all.tsv.
    > 8.0 cM: 146/146
    Of these with suff. SNPs per cM> 220:               113/146
    4     9
    2     8
    1     7
    13    7
    6     7
    8     7
    10    7
    21    6
    5     6
    7     6
    16    6
    11    5
    9     5
    12    4
    18    4
    20    4
    3     4
    14    3
    17    3
    22    3
    15    2
    Name: ch, dtype: int64
    Saved 9 individual IBD pairs to: test/ibd_ind.tsv


To view the complete options provided by the two command-line interface,
use -h. For power users or people interested in applying the method
beyond 1240k SNP set, keep in mind that one can obtain maximum
flexibility by writing one’s own wrappers (see section `prepare
input <create_hdf5_from_vcf.ipynb>`__, `run
ancIBD <run_ancIBD.ipynb>`__, and `visualization <plot_IBD.ipynb>`__)

.. code:: ipython3

    !ancIBD-run -h


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


.. code:: ipython3

    !ancIBD-summary -h


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

