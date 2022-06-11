# ancIBD
This software package screens ancient human DNA for long IBD blocks (Identity by Descent segments) shared between pairs of individuals.

## Scope of ancIBD
### Ancestry
As the method relies on imputation with a modern reference panel, it is not appicable for deeply diverged humans (such as Neanderthals or Denisovans.) However, tests showed that the method works for Eurasian ancient human DNA (tested up to 45,000 years ago). 

### Coverage
The method is relatively data-hungry, and needs imputation with Glimpse that only works well for samples with >0.5x coverage on 1240k SNPs. 

Ideally, your data has 1x coverage for 1240k SNPs, or 0.5x coverage for WGS data. For robust IBD calls of long IBD (8 cM or longer), at least 600,000 SNPs on the 1240k panel should be covered at least once. Note that there still can be occasional false positive IBD, so please always treat the output with necessary caution and not as a black box.

Generally, the shorter the IBD, the less robust the calls, and IBD shorter than 8 cM are prone to false-positive signals.


## Input Data
The starting point of the `ancIBD` pipeline is a VCF that has been inputed and phased with the software `Glimpse` (https://odelaneau.github.io/GLIMPSE). The default parameters of `ancIBD` are tuned to work well on imputed data that used the 1000G reference haplotype panel. Generally, imputation should be done on all 1000G SNPs, even if your data is 1240k capture data, and then downsampled to 1240k SNPs only, for which the parameters of ancIBD have been optimized. This VCF **needs** to contain a field for the phased genotype (GT) as well as the three genotype probabilities (as GP field). 

This full imputed VCF is then transformed into a so called .hdf5 file - which is the input for ancIBD functions to call and visualize IBD.

# Quick Start and Examples
For example uses cases that guide you through a typical data preparation and IBD screening, please walk through the Vignette notebooks in the `vignette` folder. The example notebooks and data can be downloaded from this link:
https://www.dropbox.com/sh/q18yyrffbdj1yv1/AAC1apifYB_oKB8SNrmQQ-26a?dl=0

To do so you can use `jupyter notebook` or `jupyter lab`, two powerful interfaces for interacting with Python code.

There are three main Vignette notebooks"
- How to prepare hdf5 data from a Glimpse imputed VCF
This Vignette Notebook runs you through the steps for producing the input that ancIBD needs, a so called hdf5 file.
- How to run ancIBD
This notebook walks you through the steps to call IBD segments with `ancIBD`.
It assumes one has data in hdf5 format, including genetic map and ideally also allele frequency data.
- How to visualize IBD
This notebook gives examples how to use ancIBD visualization code. Get the most out of your IBD calls by using crisp visualizations out of the box!

Beyond these three core Vignettes, there is one optional vignette notebook diving into the hdf5 file format and how to load and explore such files:
- h5_example.ipynb
Understanding HDF5 is not necessary for `ancIBD`, which loads the data automatically. This notebook is rather meant as a first starting point for the interested reader (as reading HDF5 files might be also interesting for other analyses) and for better trouble-shooting of HDF5 files.

## c Extension (for experts)
For performance reasons, the heavy lifting of the algorithm is coded into a c method (cfunc.c). This "extension" is built via cython from cfunc.pyx This should be done automatically via the package cython (as CYTHON=True in `setup.py` by default).

You can also set CYTHON=False, then the extension is compiled from cfunc.c directly (experimental, not tested on all platforms!).

## Development
The code used to develop this package is deposited at the github repository: 
https://github.com/hringbauer/ancIBD


## Citation
If you use `ancIBD` for a scientific publication, there will be a Preprint #soon.

## Contact
If you have bug reports, suggestions or general comments, please contact me. I am happy to hear from you. Bug reports and user suggestions will help me to improve this software - so please do not hesitate to reach out!

harald_ringbauer AT eva mpg de
yilei_huang AT eva.mpg.de
(fill in blanks with dots and AT with @)

Authors:
Harald Ringbauer, Yilei Huang, 2022