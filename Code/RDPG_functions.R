#library('tcltk') #for the plotting
#library('shapes') #for the procrustes' bed
library('igraph') #for the graph stuff (load it after shapes as it masks V)
#library('ppls')
library('ggplot2') #for the plotting
#library('car') #for the plotting
#library('leaps')
library('RcppEigen')
library('inline')
library('magrittr')

#overall function for leave one out analysis
do_all <- function(name,max_rank){
  Adj <- get_Adj(name,base=name)
  connectivity <- sum(Adj) / (dim(Adj)[1]^2) 
  
  names <- 1:dim(Adj)[1]

  diffs <- get_diffs(max_rank,Adj,names)

  preds <- prediction_ratio(Adj,diffs)
  write.csv(file=paste(name,"preds.csv",sep="_"),preds)
  
  pdf(paste(name,"r_preds.pdf",sep="_"),8,8)
  plot(preds$r_ones,ylim=c(0,1),pch="1",ylab="Correct prediction rate", xlab="coordinates")
  points(preds$r_zeroes,ylim=c(0,1),pch="0")
  points(preds$r_whole,ylim=c(0,1),pch=16)
  abline(h=connectivity)
  dev.off()

  return(preds)
}

#function for rank by rank leave one out analysis
do_all_rank <- function(name,rank){
  Adj <- get_Adj(name)
  names <- 1:dim(Adj)[1]

  diffs <- get_diffs_rank(rank,Adj,names)
  write.csv(file=paste(name,"_diffs_",rank,".csv",sep=""),diffs)

  return(diffs)
}

#ancillary function to list all directories
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


#Get and Adjacency matrix for the Food Web in folder 'name'
#Read readme for data format
get_Adj <- function(name,sep="\t",head=FALSE,
  Links_file='interactions.txt',Links_order=c(2,1),base=''){
    Links_file <- paste(base,Links_file,sep="/")
    Links_order <- Links_order #from prey to predator
    Csv_settings <- list(sep=sep,"head"=head)
    #Read food web links
    Links_file %>%
    read.csv(,sep=Csv_settings$sep,header=Csv_settings$head) %>%
    as.matrix %>%
    .[,Links_order] -> Links
    ###build the graph and recover adjacency matrix###
    Graph <- graph.edgelist(Links, directed=TRUE)
    Adj <- get.adjacency(Graph)
    return(as.matrix(Adj))
}

test_traits <- function(Adj, Traits_out, Traits_in){
  if( sum((Adj - round(as.matrix(Traits_out) %*% t(as.matrix(Traits_in))))^2) == 0 ){
    print("Everything ok")
  }
}


#RcppEigen code to compute SVD decomposition of a matrix
#It doesn't use LAPACK routines, so it may avoid problems with inverses
codeEigen='
  const Eigen::Map<Eigen::MatrixXd> m (as<Eigen::Map<Eigen::MatrixXd> >(m_ ));

  Eigen::JacobiSVD <Eigen::MatrixXd>svd(m,
                   Eigen::ComputeThinU|Eigen::ComputeThinV);
  return List::create( Rcpp::Named("u")=svd.matrixU(),
                       Rcpp::Named("d")=svd.singularValues(),
                       Rcpp::Named("v")=svd.matrixV() );
'
svdEigen <- cxxfunction(signature(m_="matrix"), codeEigen, plugin="RcppEigen")

#'get_traits'
#computes the functional traits for matrix 'Adj'
#'names' define the species names
#it returns a list of dataframes if 'dataframe=TRUE'
#or a list of matrices if 'dataframe=FALSE'
get_traits <- function(Adj, dataframe=TRUE, names=FALSE){
  SVD <- svdEigen(Adj) #svd decomposition fo 'Adj' in u*d*v
  U <- SVD$u
  V <- SVD$v
  S <- diag(SVD$d)
  Ssqrt <- structure(vapply(S, sqrt, numeric(1)),dim=dim(S))
  traits_in <-  V %*% Ssqrt #scales the foraging traits
  traits_out <- U %*% Ssqrt #scales the vulnerability traits
  if(dataframe){
    traits_in <- data.frame(traits_in)
    traits_out <- data.frame(traits_out)
  }
  if(names){
    row.names(traits_in) <- row.names(traits_out) <- names #set species names
    colnames(traits_in) <- colnames(traits_out) <- names #set species names
  }
  traits <- list("In"=traits_in, "Out"=traits_out)
  return(traits)
}

#'get_traits_ij'
#computes the functional traits for matrix 'Adj' with
#entry 'Adj[i,j]' unknown
#'i' and 'j' can be index-values or species names as defined in 'names'
#it returns a list of dataframes if 'df=TRUE'
#or a list of matrices if 'df=FALSE'
get_traits_ij <- function(Adj,i,j,connectivity,df=FALSE,names){
  Adj_ij <- Adj #the original matrix
  Adj_ij[i,j] <- connectivity #remove knowledge about the interaction from i to j
  Traits_ij <- get_traits(Adj_ij,dataframe=df,names=names) #compute traits
  return(Traits_ij)
}

#'get_traits_ij'
#computes the functional traits for matrix 'Adj' with
#entry 'Adj[i,j]' unknown
#'i' and 'j' can be index-values or species names as defined in 'names'
#it returns a list of dataframes if 'df=TRUE'
#or a list of matrices if 'df=FALSE'
Leave_one_out_ij <- function(Adj,i,j,connectivity,rank,names){
  len <- dim(Adj)[1]
  Traits <- get_traits_ij(as.matrix(Adj),i,j,connectivity=connectivity,df=FALSE,names)
  ij <- (Traits$Out[,1:rank] %*% t(Traits$In[,1:rank]))[i,j]
  return(ij)
}


#Produce a matrix where each entry is the leave_one_out fitted probability
#obtained by inputing that entry as a prior in Adj
Leave_one_out_mat <- function(Adj,rank,names){
  len <- dim(Adj)[1]
  connectivity <- sum(Adj) / (len^2) 
  Adj_unknown <- matrix(0,len,len)
  pb <- txtProgressBar(min = 0, max = len^2, style = 3)
  k <- 0
  for(i in 1:len){
    for(j in 1:len){
        Traits <- get_traits_ij(as.matrix(Adj),i,j,connectivity=connectivity,df=FALSE,names)
        Adj_unknown[i,j] <- (Traits$Out[,1:rank] %*% t(Traits$In[,1:rank]))[i,j]
    k <- k+1
    setTxtProgressBar(pb, k)
    }
  }
  close(pb)
  return(Adj_unknown)
}

#compute a vector of differences (signed) between the original adjacency
#and the leave one out mat, for a rank approximation up to max_rank
get_diffs <- function(max_rank,Adj,names){
  diffs <- matrix(NA,length(Adj),max_rank)
  vAdj <- as.vector(Adj)
  for(rank in 1:max_rank){
    diffs[,rank] <- as.vector(Leave_one_out_mat(Adj,rank,names)) - vAdj
  }
  return(diffs)
}

#compute a vector of differences (signed) between the original adjacency
#and the leave one out mat, for a rank approximation equal to 'rank'
get_diffs_rank <- function(rank,Adj,names){
  diffs <- rep(NA,length(Adj))
  vAdj <- as.vector(Adj)
  diffs <- as.vector(Leave_one_out_mat(Adj,rank,names)) - vAdj
  return(diffs)
}


#compute prediction rates
prediction_ratio <- function(Adj,diffs){

  #Two function to distinguish between zeroes and ones in the adjacency matrix.
  #'vAdj' is the flatten ('as.vector()') of 'Adj'
  vAdj <- as.vector(Adj)
  times_adj <- function(X){X %<>% .[as.logical(vAdj)]}
  times_nadj <- function(X){X %<>% .[as.logical(1-vAdj)]}
  sum_Adj <- sum(vAdj)
  sum_nAdj <- sum(1-vAdj)
  len_Adj <- length(vAdj)

  rdiffs_Adj <- colSums(apply(round(diffs), 2, times_adj))
  diffs_Adj <- colSums(apply(diffs, 2, times_adj))
  rdiffs_nAdj <- colSums(apply(round(diffs), 2, times_nadj))
  diffs_nAdj <- colSums(apply(diffs, 2, times_nadj))

  r_ones <- (sum_Adj  + rdiffs_Adj)/sum_Adj
  ones <- (sum_Adj  + diffs_Adj)/sum_Adj
  r_zeroes <- (sum_nAdj  - rdiffs_nAdj)/sum_nAdj
  zeroes <- (sum_nAdj  - diffs_nAdj)/sum_nAdj
  r_whole <- (len_Adj  - (rdiffs_nAdj + abs(rdiffs_Adj)))/len_Adj
  whole <- (len_Adj - (diffs_nAdj + abs(diffs_Adj)))/len_Adj
  data <- data.frame('r_ones'=r_ones,'ones'=ones,'r_zeroes'=r_zeroes,'zeroes'=zeroes,'r_whole'=r_whole,'whole'=whole)
  return(data)
}


#Mantel test of vcv.phylo and Jaccard similarity
MBJ <- function(VCV,Traits_out,Traits_in,ranks=size,random_samples=FALSE, mode="all",cumulative=TRUE){
  size <- dim(Traits_in)[1]
  if(cumulative){Raw_k_fw <- Traits_out[,1:ranks] %*% t(Traits_in[,1:ranks])}
  else {Raw_k_fw <- Traits_out[,ranks] %*% t(Traits_in[,ranks])}
  if(random_samples){
    Mantel_stats <- rep(NA,random_samples)
    Mantel_signs <- rep(NA,random_samples)
    for(samples in 1:random_samples){
      sample_fw <- matrix(ifelse(runif(size*size) > c(Raw_k_fw), 0, 1),size,size)
      sample_fw <- graph.adjacency(sample_fw)
      Jaccard_sims <- similarity.jaccard(sample_fw, mode=mode)
      m_k <- mantel(Jaccard_sims,VCV,na.rm=TRUE)
      Mantel_stats[samples] <- m_k$statistic
      Mantel_signs[samples] <- m_k$signif
    }
  } else {
    sample_fw <- round(Raw_k_fw)
    sample_fw <- graph.adjacency(sample_fw)
    Jaccard_sims <- similarity.jaccard(sample_fw, mode=mode)
    m_k <- mantel(Jaccard_sims,VCV,na.rm=TRUE)
    Mantel_stats <- m_k$statistic
    Mantel_signs <- m_k$signif
  }
  return(list("Mantel statistic" = Mantel_stats, "p-value" = Mantel_signs))
}
