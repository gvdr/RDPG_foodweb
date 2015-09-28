library(vegan)
library(phytools)
library(geiger)
library(igraph)
library(magrittr)
source("jaccard.R")
source("../Phylosig_functions.R")

name <- "Serengeti"
Rohr_i <- 0.4831
Petchey_i <- 0.5743

###input data
setwd(paste("~/Projects/wd/",name,sep=""))
Links_file <- "interactions.txt"
Links_order <- c(1,2) #from prey to predator
Csv_settings <- list("sep"="\t","head"=FALSE)

#Read food web links
Links_file %>%
read.csv(,sep=Csv_settings$sep,header=Csv_settings$head) %>%
as.matrix %>%
.[,Links_order] -> Links

###build the graph and recover adjacency matrix###
Graph = graph.edgelist(Links, directed=TRUE)
Adj <- get.adjacency(Graph)


### estimate full rank traits
names <- row.names(Adj)
Traits <- get_traits(Adj,names=names)
Traits_In <- Traits$In
Traits_Out <- Traits$Out

### performance
max_rank <- 10
randomizations <- 999
#Method 1, deterministic threshold approximation
Performance_trunk <- array(NA,c(max_rank,3))
mAdj <- as.matrix(Adj)
mTraits_Out <- as.matrix(Traits_Out)
mTraits_In <- as.matrix(Traits_In)
for(rank in 1:max_rank){
  Performance_trunk[rank,] <- F_support(mAdj,trunk_matrix(mTraits_Out,mTraits_In,rank,0.5))
}
colnames(Performance_trunk) = c("Sensitivity","Precision","Accuracy")

#Method 2, randomization
Performance_rand_rep_S <- array(NA,c(randomizations,max_rank,3))
for(rand in 1:randomizations){
for(rank in 1:max_rank){
  Performance_rand_rep_S[rand,rank,] <- F_support(mAdj,rand_matrix(mTraits_Out,mTraits_In,rank))
}
}

setEPS()
postscript(paste(name,"_accuracy_Rohr_Petchey.eps",sep=""))
boxplot(Performance_rand_rep_S[,,3],ylim=c(0,1),xlab="Trait Dimension",ylab="Accuracy")
points(Performance_trunk[,3],col="red",ylim=c(0,1),pch=20)
abline(h=Rohr_i)
abline(h=Petchey_i,lty=2)
dev.off()


setEPS()
postscript(paste(name,"_sensitivity_Rohr_Petchey.eps",sep=""))
boxplot(Performance_rand_rep_S[,,1],ylim=c(0,1),xlab="Trait Dimension",ylab="Sensitivity")
points(Performance_trunk[,1],col="red",ylim=c(0,1),pch=20)
abline(h=Rohr_i)
abline(h=Petchey_i,lty=2)
dev.off()
