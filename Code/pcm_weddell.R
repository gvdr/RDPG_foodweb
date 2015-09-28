rm(list=ls())
setwd("/scratch-local/gvd16/wd/")

Name <- "Weddell" #Name of the Food Web to analyse
Mode <- "all" #Mode of the Jaccard Similarity, either "in", "out" or "all"
cumulative <- FALSE #Cumulative (TRUE) or single rank (FALSE) matrix approximation
id_column <- 3 #the column in the ids.csv file where the names are. 

mantel_jaccard <- function(Name,Mode,cumulative,id_column){

library(phytools)
library(vegan)
library(geiger)
library(igraph)
library(magrittr)
library(reshape2)
library(ggplot2)
library(Cairo)
source("RDPG_functions.R")

Tree <-read.newick(paste(Name,"/tree.phy",sep=""))
Tree <- collapse.singles(Tree)
Ids <- read.csv(paste(Name,"/ids.csv",sep=""), header=TRUE)[id_column]
Ids <- as.character(Ids[,1])

#Weddell has some non living-species entries (detritus and such) wich is not in the names list
#so, we have to insert them by hand
if(Name == "Weddell"){
Ids <- append(Ids,c("non-organic1","non-organic2","non-organic3","non-organic4"))
}

#Let's get the Adjacency matrix from "interactions.txt"
#and those juicy traits
Adj <- get_Adj(Name,base=Name)
Traits <- get_traits(Adj, names=Ids)

#match everything against the phylogeny
Datatree_In <-  treedata(Tree,Traits$In,sort=TRUE)
matched_taxa <- c(row.names(Datatree_In$data)) #we may need this record
Datatree_Out <- treedata(Tree,Traits$Out,sort=TRUE)



###Compute the variance-covariance phylogenetic matrix
vcv <- vcv.phylo(Datatree_In$phy)


###Compute Mantel correlation between rank-k fws and vcv
random_samples <- 99
max_dim <- 20

#let compute mantel test of Jaccard similarity against vcv
#for the selected Mode

SRankvsMSig <- matrix(NA,random_samples,max_dim)
SRankvsPSig <- matrix(NA,random_samples,max_dim)
colnames(SRankvsMSig) <- colnames(SRankvsPSig) <- as.character(1:max_dim)
row.names(SRankvsMSig) <- rep("Mantel test",random_samples)
row.names(SRankvsPSig) <- rep("p-values",random_samples)

for(rank in 1:max_dim){
mbj_rank <-MBJ(mode=Mode,vcv,Datatree_Out$data,Datatree_In$data,ranks=rank,random_samples=random_samples,cumulative=cumulative) #Mantel by Jaccard!
SRankvsMSig[,rank] <- unlist(mbj_rank[1])
SRankvsPSig[,rank] <- unlist(mbj_rank[2])
}

SRankM <- melt(SRankvsMSig)
SRankP <- melt(SRankvsPSig)
names(SRankM) <- names(SRankP) <- c("Type","Rank","Value")

SRank <- as.data.frame(rbind(SRankM,SRankP))
SRank$Rank <- as.factor(SRank$Rank)
SRank$Type <- as.factor(SRank$Type)

#plot the results
pdf_name <- if(cumulative){paste(Name,"_Cumulative_",Mode,".pdf",sep="")} else {paste(Name,"_SingleRanks_",Mode,".pdf",sep="")}
CairoPDF(pdf_name,8,8)
ggplot(data=SRank, aes(x=Rank, y=Value, colour=Type)) + 
scale_x_discrete() +
geom_jitter(alpha=0.75,position = position_dodge(width = .75)) +
geom_boxplot(outlier.colour=NA) + geom_hline(yintercept=0.01) +
theme_bw()
dev.off()

}

