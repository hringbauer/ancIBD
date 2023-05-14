# Frequently asked Questions

  In the plotting function "plot_pde_individual_from_ibd_df", is it possible to plot the theoretical distribution for aunt-nephew and siblings?

Yes it is! You have to update the parameters according - and change the following default parameters accordingly:

comm_ancs =[4,4,2,2]
ms=[4,6,5,4]
labels=["First Cousins", "Second Cousins", "5 generations anc.", "4 generations and."]

`comma_ancs ` should be 4 (as in siblings and aunt nephew the first generation has 2x2 haplotypes). For siblings there are 2 meiosis seperating the target haplotypes from the common ancestor haplotypes (therefore `ms =2`; for aunt-nephew there are 3 (and therefore `ms=3`). You c should also update the labels and optionally also colors (cs) - matching the respective relationship.

Doing such updates, one can plot the expected IBD relationships of various kinds of genealogical relationships, you have to update  `comm_ancs `,  `ms ` (and  `labels ` and optionaly  `cs `) accordingly.
