Overview
============

The package ``ancIBD`` detects Identity-by-descent (IBD) segments in typical human aDNA data. It takes as input imputed genotype data by `GLIMPSE <https://odelaneau.github.io/GLIMPSE/glimpse1/index.html>`_.

Scope
**********

``ancIBD`` is a comparably data-hungry method but can be applied to a substantial fraction of the aDNA record. Our tests showed that ``ancIBD`` requires at least 0.25x average coverage depth for whole-genome-sequencing (WGS) data and 1.0x average coverage on target SNPs (corresponding broadly to at least 600k SNPs covered for 1240k or TWIST captured aDNA data, two popular SNP captures in human aDNA). When one gets close to that coverage limit,  imputation starts to break down and false positive IBD rates increase, in particular for shorter IBD segments. Inferred IBD segments for data below that coverage limit have to be interpreted with extreme caution, as false positive and error rates become substantial and likely dominate any true signal in most demographic scenarios.

We recommend imputing ancient data using the software GLIMPSE, imputing ancient samples one by one. The default parameters of ``ancIBD`` are optimized for data imputed using the modern 1000 Genome reference panels and all SNPs in this reference panel, and then downsampling to the so-called 1240k SNP set widely used in human ancient DNA. 

As ``ancIBD`` relies on imputed data, it works well for up to several ten-thousands year old modern human genomes when using the present-day 1000 Genome reference panel. We have observed that ``ancIBD`` performs well for global ancient genomes sharing the out-of-Africa bottleneck (i.e. modern humans from Eurasia, Oceania, Americas), however, some Sub-Saharan ancestries can be problematic as they contain deeply diverged haplotypes that are not represented well in the 1000 Genome reference panel and are imputed exceptionally poorly.

Citing
**********

You can find a pre-print describing the method and several applications here:

`ancIBD - Screening for identity by descent segments in human ancient DNA <https://doi.org/10.1101/2023.03.08.531671>`_

You can cite this article if you use ``ancIBD`` for your scientific work.

Example Data
**********

You can find the test aDNA data used to run the tutorials described throughout this documentation `here <https://www.dropbox.com/sh/q18yyrffbdj1yv1/AAC1apifYB_oKB8SNrmQQ-26a?dl=0>`_. This Dropbox folder contains GLIMPSE-imputed .vcf files of a subset of samples from early Neolithic Britain that are part of a published extended pedigree (`Fowler, Olalde et al. 2021 <https://www.nature.com/articles/s41586-021-04241-4>`__).

Using ``ancIBD`` via Python functions
**********

One can run ``ancIBD`` using Python functions, and we provide example Jupyter notebooks on how to:

-   `Prepare the data <create_hdf5_from_vcf.ipynb>`__
-   `Call IBD with ancIBD <run_ancIBD.ipynb>`__
-   `Visualize the IBD output <plot_IBD.ipynb>`__

Users can embed the underlying Python functions into Python wrapper scripts and their own interactive Jupyter notebooks.

Using ``ancIBD`` via the command line (available in the next release)
**********

Alternatively, one can also run ``ancIBD`` directly from the command line. These commands are added during the installation of the Python package. You can find a detailed walk-through in the section `Running ancIBD via bash <quick_start_bash.rst>`__.

Contact
**********

If you have bug reports, suggestions, or any general comments please do not hesitate to reach out - we are happy to hear from you! Your suggestions will help us to improve this software.

You can report bugs as an issue on the official ``ancIBD`` `github development page <https://github.com/hringbauer/ancIBD>`_

We are also happy to hear from you via email:

-   harald_ringbauer AT eva mpg de
-   yilei_huang AT eva mpg de

(fill in AT with @ and other blanks with dots)


Lead Authors:
Harald Ringbauer, Yilei Huang, 2023
