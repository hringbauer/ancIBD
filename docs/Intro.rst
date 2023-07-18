Overview
============

``ancIBD`` is a Python package to detect Identity-by-descent (IBD) segments in typical human aDNA data. It takes as input imputed genotype data by `GLIMPSE <https://odelaneau.github.io/GLIMPSE/glimpse1/index.html>`_.

Scope
**********

``ancIBD`` is a comparably data-hungry method but can be applied to a substantial fraction of the aDNA record. Our tests showed that ``ancIBD`` requires at least 0.25x average coverage depth for whole-genome-sequencing (WGS) data and 1.0x average coverage on target SNPs (corresponding broadly to at least 600k SNPs covered for 1240k or TWIST captured aDNA data, two popular SNP captures in human aDNA). When one gets close to that coverage limit,  imputation starts to break down and false positive IBD rates increase, in particular for shorter IBD segments. Inferred IBD segments for data below that coverage limit have to be interpreted with extreme caution, as false positive and error rates become substantial and likely dominate any true IBD signal in most demographic scenarios.

Citing
**********

Please find a pre-print describing the method and several applications here:

`ancIBD - Screening for identity by descent segments in human ancient DNA <https://doi.org/10.1101/2023.03.08.531671>`_

You can cite this article if you use ``ancIBD`` for your scientific work.

Contact
**********

If you have bug reports, suggestions, or any general comments please do not hesitate to reach out - we are happy to hear from you! Your suggestions will help us to improve this software.

You can report bugs as an issue on the official ``ancIBD`` `github development page <https://github.com/hringbauer/ancIBD>`_

We are also happy to hear from you per email:

-harald_ringbauer AT eva mpg de
-yilei_huang AT eva mpg de

(fill in AT with @ and other blanks with dots)


Lead Authors:
Harald Ringbauer, Yilei Huang, 2023
