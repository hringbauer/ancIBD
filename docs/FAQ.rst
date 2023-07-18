Frequently Asked Questions (FAQ)
==================================


Can I call IBD on X chromosomes?
*********************************
In theory you can. Please follow `GLIMPSE's tutorial for chrX <https://odelaneau.github.io/GLIMPSE/glimpse1/tutorial_chrX.html>`_ for imputing the X chromosome. Then you can use the same pipeline as for the autosomes. 

However, we note that performance of ancIBD on X chromosome is not benchmarked and likely it will benefit from parameter fine-tuning as the X chromosome has a different demography from autosomes. And we note that the model in ancIBD is only suitable for calling IBD between females because the model assumes that phasing switch error occur at the same rate for both samples.
For calling IBD between males or between a male and a female, stay tuned for updates!


Do I need to mask certain genomic regions?
******************************************
In most cases you don't need to apply genomic mask as the SNP density filter (>220 SNPs per cM) can effectively remove most false positive IBD segments.

However, if you are interested in using shorter segments (e.g, <8cM), a mask file to remove regions with low SNP density is recommended. A mask file is provided in the folder “./data/map” in the `dropbox link <https://www.dropbox.com/sh/q18yyrffbdj1yv1/AAC1apifYB_oKB8SNrmQQ-26a?dl=0>`_.
These masks are based on false positive rates along the genome out of the 10,156 Eurasia samples described in our manuscript. But you can also use your own mask file that suits your purpose.


In the plotting function "plot_pde_individual_from_ibd_df", is it possible to plot the theoretical distribution for aunt-nephew and siblings?
***********************************************************************************************************************************************
Yes it is! You have to update the parameters according - and change the following default parameters accordingly:

comm_ancs =[4,4,2,2]
ms=[4,6,5,4]
labels=["First Cousins", "Second Cousins", "5 generations anc.", "4 generations and."]

`comma_ancs ` should be 4 (as in siblings and aunt nephew the first generation has 2x2 haplotypes). For siblings there are 2 meiosis seperating the target haplotypes from the common ancestor haplotypes (therefore `ms =2`; for aunt-nephew there are 3 (and therefore `ms=3`). You c should also update the labels and optionally also colors (cs) - matching the respective relationship.

Doing such updates, one can plot the expected IBD relationships of various kinds of genealogical relationships, you have to update  `comm_ancs `,  `ms ` (and  `labels ` and optionaly  `cs `) accordingly.
