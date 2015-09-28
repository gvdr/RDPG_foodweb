library(plyr)
library(igraph)
source("RDPG_functions.R")

base_all <- "../Data/Elaborated/Serengeti/Allesina"
base_raw <- "../Data/Raw/Serengeti"
name <- "Serengeti"
Adj <- get_Adj("Serengeti",base=base_raw)
part_file <- paste(base_all,"partitions.csv",sep="/")
prob_file <- paste(base_all,"groups_ps.csv",sep="/")

partitions <- read.csv(part_file,sep="\t")
probabilities <- read.csv(prob_file,sep="\t")
names_ord <- as.character(partitions[,1])
Adj <- t(Adj[names_ord,names_ord])

names <- row.names(Adj)
Traits <- get_traits(Adj,names=names)
Traits_In <- Traits$In
Traits_Out <- Traits$Out
SerSt <- matrix(0,14,14)

Size <- dim(probabilities)[1] 
for(i in 1:Size){
SerSt[probabilities[i,1],probabilities[i,2]] <- probabilities[i,3]
}

block_sizes <- count(partitions,"group")$freq
N_species <- sum(block_sizes)
Perf <- matrix(0,100,3)
for(i in 1:100){
test1 <- sbm.game(N_species, SerSt, block_sizes, directed=1,loops=1) 
test_mat <- as.matrix(get.adjacency(test1))
Sens <- sum(test_mat * Adj) / sum(Adj)
Spec <- sum((1-test_mat) * (1 - Adj)) / sum(1-Adj)
Accuracy <- (sum(test_mat * Adj) + sum((1-test_mat) * (1 - Adj))) / (dim(Adj)[1]^2)
Perf[i,] <- c(Sens, Spec, Accuracy)
}

rand_matrix <- function(Traits_out,Traits_in,rank){
  Size <- dim(Traits_out)[1]
  Prob <- as.matrix(Traits_out[,1:rank]) %*% as.matrix(t(Traits_in[,1:rank]))
  Adj_out <- t(matrix(ifelse(runif(Size^2) > t(Prob), 0,1),Size,Size))
  return(Adj_out)
}

trunk_matrix <- function(Traits_out,Traits_in,rank){
  Size <- dim(Traits_out)[1]
  Prob <- as.matrix(Traits_out[,1:rank]) %*% as.matrix(t(Traits_in[,1:rank]))
  Adj_out <- t(matrix(ifelse(0.5 > t(Prob), 0,1),Size,Size))
  return(Adj_out)
}

F_support <- function(Obs,Rand){
Size <-  dim(Obs)[1]
Sens <- sum(Rand * Obs) / sum(Obs)
Spec <- sum((1-Rand)*(1-Obs)) / sum(1-Obs)
Acc <- (sum(Rand * Obs) + sum((1-Rand)*(1-Obs))) / (Size^2)
return(c(Sens,Spec,Acc))
}

Traits <- get_traits(Adj)
max_rank <- 20
randomizations <- 999
#Method 2, randomization
Performance_rand_rep_S <- array(NA,c(randomizations,max_rank,3))
for(rand in 1:randomizations){
for(rank in 1:max_rank){
  Performance_rand_rep_S[rand,rank,] <- F_support(Adj,rand_matrix(Traits$Out,Traits$In,rank))
}
}

Performance_trunk <- array(NA,c(max_rank,3))
for(rank in 1:max_rank){
  Performance_trunk[rank,] <- F_support(Adj,trunk_matrix(Traits$Out,Traits$In,rank))
}

pdf(paste("../Manuscript/images/fitting_perf/accuracy/",name,"_zoom.pdf",sep=""),6,6)
boxplot(Performance_rand_rep_S[,,3],ylim=c(0.9,1),xlab="Trait Dimension",ylab="Accuracy")
points(Performance_trunk[,3],col="red",ylim=c(0.9,1),pch=20)
abline(h=min(Perf[,3]))
abline(h=max(Perf[,3]))
dev.off()


pdf(paste("../Manuscript/images/fitting_perf/sensitivity/",name,".pdf",sep=""),6,6)
boxplot(Performance_rand_rep_S[,,1],ylim=c(0,1),xlab="Trait Dimension",ylab="Sensitivity")
points(Performance_trunk[,1],col="red",ylim=c(0,1),pch=20)
abline(h=min(Perf[,1]))
abline(h=max(Perf[,1]))
dev.off()

Adj_prob <- Adj
for(name_from in as.character(partitions[,1])){
	for(name_to in as.character(partitions[,1])){
		Adj_prob[name_from,name_to] <- SerSt[partitions[partitions[,1]==name_from,2],partitions[partitions[,1]==name_to,2]]
	}
}

Connect <- sum(Adj)/(dim(Adj)[1]^2)
Adj_c <- Adj
for(name_from in as.character(partitions[,1])){
	for(name_to in as.character(partitions[,1])){
		Adj_c[name_from,name_to] <- Connect
	}
}


