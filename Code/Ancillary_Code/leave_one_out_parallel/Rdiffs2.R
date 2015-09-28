rm(list=ls())
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
print(args)
name <- args[1]

id <- as.integer(args[2])

rank <- 16 + id %/% 492
j <- 1 +  id %% 492

#Leave one out test for rank k and entry (i,j)
source('../Phylosig_functions.R',echo = TRUE)
for(i in 1:492){
	do_ij_rank(name,rank,i,j)
}

