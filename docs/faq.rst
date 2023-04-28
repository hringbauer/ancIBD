# Frequently asked Questions

  In the plotting function "plot_pde_individual_from_ibd_df", is it possible to plot the theoretical distribution for aunt-nephew and siblings?

Yes it is! You have to update the parameters according - and change the following default parameters accordingly:

comm_ancs =[4,4,2,2]
ms=[4,6,5,4]
labels=["First Cousins", "Second Cousins", "5 generations anc.", "4 generations and."]

 `comma_ancs ` should be four (as in siblings and aunt nephew the first generation has 2x2 haplotypes). For siblings you have two meiosis ( `ms ` should be 4 then) and four aunt-nephew three (ms should be three). And update the labels and maybe also colors (cs) accordingly!

Equivalently, you can plot the expectations of all kind of relationships, you have to update  `comm_ancs `,  `ms ` (and  `labels ` and optionaly  `cs `) accordingly.
