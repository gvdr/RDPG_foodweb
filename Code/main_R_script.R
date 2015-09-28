setwd("H:/Projects/Stouffer_RDPFW/Data/Elaborated")
#read the edgelist (the order should be "prey" to "predator")
links = as.matrix(read.csv("H:/Projects/Stouffer_RDPFW/Data/Raw/W_Links.csv"))[,c(2,1)]
ids = as.matrix(read.csv("H:/Projects/Stouffer_RDPFW/Data/Raw/W_ids.csv"))
Name = "Weddell Sea"

#set randomizations parameters (for accuracy estimation)
randomizations = 999
connected_only = F
max_rank = 15


#load the library in this order: there is a masking
library('tcltk') #for the plotting
library('shapes') #for the procrustes' bed
library('igraph') #for the graph stuff (load it after shapes as it masks V)
library('ppls')
library('ggplot2') #for the plotting
library('car') #for the plotting
library('leaps')




#build the graph and recover adjacency matrix
Graph = graph.edgelist(Links, directed=TRUE)

#give dull names to non living things in Weddell
labels <- Ids
labels[489] <- c("STUFF1")
labels[490] <- c("STUFF2")
labels[491] <- c("STUFF3")
labels[492] <- c("STUFF4")
V(Graph)$name <- labels

n=vcount(Graph)
Adj = get.adjacency(Graph)
names = row.names(Adj)

#compute Singular Value Decomposition
SVD = svd(Adj)
U = SVD$u
S = diag(SVD$d)
Ssqrt = structure(vapply(S, sqrt, numeric(1)),dim=dim(S))
V = SVD$v


#get elbows

getElbows(SVD$d[SVD$d > quantile(SVD$d,0.05)])


#compute full rank in- and out- traits
traits_in =  V %*% Ssqrt
traits_out = U %*% Ssqrt
row.names(traits_in) = row.names(traits_out) = names
traits_in_df = data.frame(traits_in)
traits_out_df = data.frame(traits_out)


#Have we done everything right?
if( sum((Adj - round(traits_out %*% t(traits_in)))^2) == 0 ) {print("Everything ok")}
#we have to round the estimated matrix due to numerical limits

mean_ds <- array(NA,max_rank)
var_ds <- array(NA,max_rank)
dists_nonrand <- array(NA,max_rank)

#Compute fitting performance statistics and write them as csv

#Method 1, deterministic threshold approximation
Performance_trunk <- array(NA,c(max_rank,3))
for(rank in 1:max_rank){
  Performance_trunk[rank,] <- F_support(Adj,trunk_matrix(traits_out,traits_in,rank,0.5))
}
colnames(Performance_trunk) = c("Sensitivity","Precision","Accuracy")

write.csv(Performance_trunk, file = paste(Name,"_trunk_performance.csv",sep=""))

#Method 2, randomization
Performance_rand_rep_S <- array(NA,c(randomizations,max_rank,3))
for(rand in 1:randomizations){
for(rank in 1:max_rank){
  Performance_rand_rep_S[rand,rank,] <- F_support(Adj,rand_matrix(traits_out,traits_in,rank))
}
}

write.csv(Performance_rand_rep_S, file = paste(Name,"_rand_rep_performance.csv",sep=""))