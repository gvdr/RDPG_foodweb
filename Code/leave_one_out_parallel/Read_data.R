#R
#let's read up all the data
library(igraph)
library(magrittr)
library(ggplot2)
library(reshape2)
#setwd("~/Projects/wd/Weddell/rdpg/Weddell/") #Set the working directory
source("RDPG_functions.R") #load the standard function

max_rank <- 26 #set the makimum rank used in the analysis, TODO infer it from existing folders
base <- "/Weddell_postprob_"
Adj <- get_Adj(name="Weddell",base="/scratch-local/gvd16/Weddell") #get original Adjacency matrix from raw data
Adj_list <- replicate(max_rank, matrix(NA,492,492), simplify=FALSE) #prepare a list of empty matrices, each for every rank
for(rank in 1:max_rank){
	for(i in 1:dim(Adj)[1]){
		for(j in 1:dim(Adj)[2]){
			file_rij <- paste("output/",rank,base,i,"_",j,".csv",sep="")
			if(file.exists(file_rij)){
				Adj_list[[rank]][i,j] <- read.csv(file_rij)$x #fill them with the a posteriori probability of interaction
			}
		}
	}	
}

get_diffs <- function(max_rank,Adj,names){
	diffs <- array(NA,dim=c(length(Adj),2,max_rank))
	vAdj <- as.vector(Adj)
	for(rank in 1:max_rank){
		diffs[,1,rank] <- as.vector(Adj_list[[rank]]) - vAdj
		diffs[,2,rank] <- 1 - is.na(as.vector(Adj_list[[rank]]))
	}
	return(diffs)
}


#compute prediction rates
prediction_ratio_rank <- function(Adj,diffs,rank){
	diffs_ranked <- diffs[,1,rank]
	existings <- as.logical(diffs[,2,rank])
	#Two function to distinguish between zeroes and ones in the adjacency matrix.
	#'vAdj' is the flatten ('as.vector()') of 'Adj'
	vAdj <- as.vector(Adj)
	times_adj <- function(X){X %<>% .[as.logical(vAdj * existings)]}
	times_nadj <- function(X){X %<>% .[as.logical((1-vAdj) * existings)]}
	sum_Adj <- sum(vAdj[as.logical(existings)])
	sum_nAdj <- sum(1-vAdj[as.logical(existings)])
	len_Adj <- length(vAdj[as.logical(existings)])

	rdiffs_Adj <- sum(round(times_adj(diffs_ranked)))
	diffs_Adj <- sum(times_adj(diffs_ranked))
	rdiffs_nAdj <- sum(round(times_nadj(diffs_ranked)))
	diffs_nAdj <- sum(times_nadj(diffs_ranked))

	r_ones <- (sum_Adj  + rdiffs_Adj)/sum_Adj
	ones <- (sum_Adj  + diffs_Adj)/sum_Adj
	r_zeroes <- (sum_nAdj  - rdiffs_nAdj)/sum_nAdj
	zeroes <- (sum_nAdj  - diffs_nAdj)/sum_nAdj
	r_whole <- (len_Adj  - (rdiffs_nAdj + abs(rdiffs_Adj)))/len_Adj
	whole <- (len_Adj - (diffs_nAdj + abs(diffs_Adj)))/len_Adj
	data <- c(r_ones,ones,r_zeroes,zeroes,r_whole,whole)
	return(data)
}

prediction_ratio <- function(Adj,diffs,max_rank){
	void <- rep(NA,max_rank)
	data <- data.frame('r_ones'=void,'ones'=void,'r_zeroes'=void,'zeroes'=void,'r_whole'=void,'whole'=void)
	data$rank <-c(1:max_rank)
	for(rank in 1:max_rank){
		data_rank <- prediction_ratio_rank(Adj,diffs,rank)
		data$r_ones[rank] <- data_rank[1]
		data$ones[rank] <- data_rank[2]
		data$r_zeroes[rank] <- data_rank[3]
		data$zeroes[rank] <- data_rank[4]
		data$r_whole[rank] <- data_rank[5]
		data$whole[rank] <- data_rank[6]
	}
	return(data)
}


W_diffs <- get_diffs(max_rank,Adj,as.character(1:max_rank))

W_P_D <- prediction_ratio(Adj,W_diffs,max_rank)

melted = melt(W_P_D, id.vars="rank")

write.csv(melted, "W_P_D_melted.csv")

pdf("Weddell_leave_one_out_7_April.pdf", 6,6)
ggplot(data=melted, aes(x=rank, y=value, group=variable, colour=variable)) + geom_line()
dev.off()
