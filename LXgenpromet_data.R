
library(openxlsx)

gene_meta_example <- read.csv("gene_meta_ORA_results.csv")
pro_meta_example <- read.csv("protein_meta_ORA_results.csv")

usethis::use_data(gene_meta_example,overwrite = T)
usethis::use_data(pro_meta_example,overwrite = T)

rm(list=ls())

data(gene_meta_example) # a list of protein Uniprot ID
data(pro_meta_example) # getting from http://impala.molgen.mpg.de/

