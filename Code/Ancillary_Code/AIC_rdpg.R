max_rank <- 25


setwd("./Data/Raw/")

library(vegan)
library(phytools)
library(geiger)
library(igraph)
library(magrittr)

source("../../Code/RDPG_functions.R")


list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
  full.names=FALSE, ignore.case=FALSE) {
  # use full.names=TRUE to pass to file.info
  all <- list.files(path, pattern, all.dirs,
           full.names=TRUE, recursive=FALSE, ignore.case)
  dirs <- all[file.info(all)$isdir]
  # determine whether to return full names or just dir names
  if(isTRUE(full.names))
    return(dirs)
  else
    return(basename(dirs))
}

bern_like <- function(Observed,Probabilities){
  obs <- log(Probabilities)
  non_obs <- log(1 - Probabilities)
  obs_sum <- sum(Observed * obs, na.rm=TRUE)
  non_obs_sum <- sum((1 - Observed) * non_obs, na.rm=TRUE)
  return(obs_sum + non_obs_sum)
}

get_probs <- function(Traits_out,Traits_in,rank){
  Probabilities <- Traits_out[,1:rank] %*% t(Traits_in[,1:rank])
  return(Probabilities)
}

par_pen <- function(Adjacency,rank_length){
	penalty <- 2 * 2 * dim(Adjacency)[1] * rank_length
	return(penalty)
}

par_pen_c <- function(Adjacency,rank_length){
	penalty <- 2 * 2 * dim(Adjacency)[1] * rank_length
	twokkplusone <- 2 * 2 * dim(Adjacency)[1] * rank_length * ( 2 * dim(Adjacency)[1] * rank_length + 1)
	nminuskminusone <- (dim(Adjacency)[1])^2 - 2 * dim(Adjacency)[1] * rank_length - 1
	return(penalty + twokkplusone/nminuskminusone)
}

get_AICc <- function(Adjacency,Traits,ranki){
  Probs <- get_probs(as.matrix(Traits$Out),as.matrix(Traits$In),ranki)
  loglike <- bern_like(Adjacency,Probs)
  pen <- par_pen_c(Adjacency,ranki)
  AIC_score <- pen - 2 * loglike
  return(AIC_score)
}

get_LikeLi <- function(Adjacency,Traits,ranki){
  Probs <- get_probs(as.matrix(Traits$Out),as.matrix(Traits$In),ranki)
  loglike <- bern_like(Adjacency,Probs)
  return(loglike)
}
get_BIC <- function(Adjacency,Traits,ranki){
  Probs <- get_probs(as.matrix(Traits$Out),as.matrix(Traits$In),ranki)
  loglike <- bern_like(Adjacency,Probs)
  pen <- 4 * log(dim(Adjacency)[1]) * dim(Adjacency)[1] * ranki
  AIC_score <- pen - 2 * loglike
  return(AIC_score)
}

get_AIC <- function(Adjacency,Traits,ranki){
  Probs <- get_probs(as.matrix(Traits$Out),as.matrix(Traits$In),ranki)
  loglike <- bern_like(Adjacency,Probs)
  pen <- par_pen(Adjacency,ranki)
  AIC_score <- pen - 2 * loglike
  return(AIC_score)
}

AIC_score <- data.frame(
	Web_name   = character(max_rank*length(list.dirs())),
	Rank_model = integer(max_rank*length(list.dirs())),
	AIC_score  = double(max_rank*length(list.dirs())),
	AICc_score  = double(max_rank*length(list.dirs())),
	BIC_score  = double(max_rank*length(list.dirs())),
	stringsAsFactors=FALSE
	)
i = 1
for(dir_name in list.dirs()){
  Adj <- get_Adj(dir_name,base=dir_name)
  Traits <- get_traits(Adj)
  for(ranki in 1:max_rank){
    AICscore <- get_AIC(Adj,Traits,ranki)
    AICcscore <- get_AICc(Adj,Traits,ranki)
    BICscore <- get_BIC(Adj,Traits,ranki)
    if(AICscore == Inf){print(paste(dir_name,ranki," fucked up"))}
    AIC_score[i,] <- c(dir_name, as.integer(ranki),
      AICscore, AICcscore, BICscore)
    i = i+1
  }
}

Likelihoods <- data.frame(
	Web_name   = character(max_rank*length(list.dirs())),
	Rank_model = integer(max_rank*length(list.dirs())),
	LogLike  = double(max_rank*length(list.dirs())),
	stringsAsFactors=FALSE
	)
i = 1
for(dir_name in list.dirs()){
  Adj <- get_Adj(dir_name,base=dir_name)
  Traits <- get_traits(Adj)
  for(ranki in 1:max_rank){
    LogLikef <- get_LikeLi(Adj,Traits,ranki)
    Likelihoods[i,] <- c(dir_name, as.integer(ranki),LogLikef)
    i = i+1
  }
}

writetable("AIC.csv",AIC_score,sep=",")
