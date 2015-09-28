library(vegan)
library(phytools)
library(geiger)
library(igraph)
library(magrittr)
library(Matrix)

source("RDPG_functions.R")

for(dir_name in list.dirs()){

adj <- get_Adj(dir_name, base=dir_name)
svds <- svd(adj)$d

rank <- rankMatrix(adj, method="tolNorm2")[1]

pdf(file=paste(dir_name,"ds.pdf",sep="_"),8,8)
plot(svds,main=paste("Rank = ",rank,sep=""))
dev.off()

}
