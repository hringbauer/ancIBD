Overview
============

``ancIBD`` is a Python package for detecting Identity-by-descent(IBD) segments for low-coverage aDNA data. It takes as input imputed genotype by `GLIMPSE <https://odelaneau.github.io/GLIMPSE/glimpse1/index.html>`_.

Scope
**********

According to our tests, ``ancIBD`` requires at least 0.25x average coverage depth for whole-genome-sequencing (WGS) data and for 1.0x average coverage (roughly equiavalent to 600k SNPs captured) on target SNPs for 1240k or TWIST captured DNA data (common SNP captures in human aDNA). IBD calls for aDNA data below that coverage limit have to be interpreted with extreme caution, as false positive and error rates become substantial and likely dominate any true IBD signal.

Citing
**********

A pre-print is in preparation and will be available soon.


Contact
**********

If you have bug reports, suggestions or any general comments please do not hesitate to reach out. We are happy to hear from you! Bug reports and user suggestions will help us to improve this software.

You can report bugs as an issue on the official  `github development page <https://github.com/hringbauer/ancIBD>`_

We are also happy to hear from you directly:
- harald_ringbauer AT eva mpg de
- yilei_huang AT eva mpg de

(fill in AT with @ and other blanks with dots)


Authors:
Harald Ringbauer, Yilei Huang, 2023
