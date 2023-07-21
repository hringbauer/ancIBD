Overview
============

The software package ``ancIBD`` detects Identity-by-descent (IBD) segments in typical human aDNA data. It takes as input imputed genotype data. The default parameters of ``ancIBD`` are optimized for imputation using the software `GLIMPSE <https://odelaneau.github.io/GLIMPSE/glimpse1/index.html>`, downsampled to the 1240k SNP set.

Scope
**********

``ancIBD`` can be applied to a substantial fraction of the aDNA record but as it relies on imputation is comparably data-hungry. Our tests showed that ``ancIBD`` requires 

- at least 0.25x average coverage depth for whole-genome-sequencing (WGS) data 
- 1.0x depth on dense target SNPs (corresponding broadly to at least 600k SNPs covered for 1240k or TWIST captured aDNA data, two popular SNP captures in human aDNA)

Close to that coverage limit, imputation starts to break down and false positive IBD rates quickly increase. Inferred IBD segments for data below that coverage limit have to be interpreted with extreme caution, as false positive and error rates become substantial. Generally, the shorter the IBD, the less robust the calls. The minimum output IBD length is 8 centimorgan (cM), but we note that already IBD shorter than 12 cM are enriched for false-positive IBD segments. Please always treat the output with necessary caution and not as a black box.

``ancIBD`` relies on imputation with a haplotype reference panel. We observed that using the present-day 1000 Genome reference panel results in robust IBD calls for up to several ten-thousands year-old human genomes for global  `homo sapiens` ancient genomes sharing the out-of-Africa bottleneck (i.e. from Eurasia, Oceania, and Americas), however, some Sub-Saharan ancestries can be problematic as they contain deeply diverged haplotypes that are not represented in the 1000 Genome reference panel. Currently, ``ancIBD`` is not applicable to other humans such as Neanderthals or Denisovans due to the lack of a suitable haplotype reference panel for imputation.

Preprocessing the data
**********
We recommend imputing ancient data using the software `GLIMPSE <https://odelaneau.github.io/GLIMPSE/glimpse1/index.html>`_, imputing ancient samples one by one as described `in its tutorial <https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_b38.html>`_. The default parameters of ``ancIBD`` are optimized for data imputed using the modern 1000 Genome reference panels and all SNPs in this reference panel, and then downsampling to the so-called 1240k SNP set widely used in human ancient DNA. The imputed 1240k SNP VCF **needs** to contain a field for the phased diploid genotypes (GT) as well as the three genotype probabilities (as GP field). This imputed VCF is then transformed into a so-called .hdf5 file - which is the input for ``ancIBD`` functions to call and visualize IBD.

Citing
**********

A pre-print that describes ``ancIBD`` and several applications is available here:

`ancIBD - Screening for identity by descent segments in human ancient DNA <https://doi.org/10.1101/2023.03.08.531671>`_

You can cite this article if you use ``ancIBD`` for your scientific work.

Example Data
**********

You can find the test aDNA data used to run the tutorials described throughout this documentation `here <https://www.dropbox.com/sh/q18yyrffbdj1yv1/AAC1apifYB_oKB8SNrmQQ-26a?dl=0>`_. This Dropbox folder contains GLIMPSE-imputed .vcf files of a subset of samples from early Neolithic Britain that are part of a published extended pedigree (`Fowler, Olalde et al. 2021 <https://www.nature.com/articles/s41586-021-04241-4>`__).

Using ``ancIBD`` via Python functions
**********

One can run ``ancIBD`` using Python functions that are imported from the package. We provide example Jupyter notebooks on how to:

-   `Prepare the data <create_hdf5_from_vcf.ipynb>`__
-   `Call IBD with ancIBD <run_ancIBD.ipynb>`__
-   `Visualize the IBD output <plot_IBD.ipynb>`__

The example notebooks and data can be also downloaded `here <https://www.dropbox.com/sh/q18yyrffbdj1yv1/AAC1apifYB_oKB8SNrmQQ-26a?dl=0Users1>`_ Users can modify hose functions and embed them into Python wrapper scripts or their own interactive Jupyter notebooks. 

Using ``ancIBD`` via the command line (available in the next release)
**********

Alternatively, one can also run ``ancIBD`` directly from the command line. These commands are automatically added during the installation of the Python package. You can find a detailed walk-through in the section `Running ancIBD via bash <quick_start_bash.rst>`__.

Contact
**********

If you have bug reports, suggestions, or any general comments please do not hesitate to reach out - we are happy to hear from you! Your suggestions will help us to improve this software.

You can report bugs as an issue on the ``ancIBD`` `GitHub page <https://github.com/hringbauer/ancIBD>`_

We are also happy to hear from you via email:

-   harald_ringbauer AT eva mpg de
-   yilei_huang AT eva mpg de

(fill in AT with @ and other blanks with dots)


Lead Authors:
Harald Ringbauer, Yilei Huang, 2023
