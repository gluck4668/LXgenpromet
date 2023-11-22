install.packages("devtools")

library(devtools)

install_github("gluck4668/LXgenpromet")

library(LXgenpromet)

#--------------------------------
data(gene_meta_example)
data(pro_meta_example)

#--------------------------------

rm(list=ls())

if(!is.null(dev.list()))
  dev.off()

gene_meta ="gene_meta_ORA_results.csv"
pro_meta ="protein_meta_ORA_results.csv"

LXgenpromet (gene_meta,gene_meta)
