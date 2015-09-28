#!Kmult! compute the observed Blomberg's K
#for an observed trait matrix !x!
#on an observed phylogeny !phy!
#it takes as input the ancillary variables
# !C!, !Ccol!, !Css!, !D.mat!, !K.denom! so to avoid repeated computations
# as the only variable changing in the function !Test.Kmult! is the row order of !x!

Kmult = function(x, phy, C, Ccol, Css, D.mat, K.denom){  
  #compute the observed variance covariance structure
  a.obs = Ccol%*%x/Css  #as for APE's evol.vcv() implementation
  
  #compute the distance between the observed and phylogenetic aware vcv
  distmat = as.matrix(dist(rbind(as.matrix(x),a.obs)))
  
  #sum distances root vs. tips
  MSEobs.d = sum(distmat[(1:N),(N+1)]^2)
  
  dist.adj = as.matrix(dist(rbind((D.mat%*%(x-(ones%*%a.obs))),0)))
  MSE.d = sum(dist.adj[(1:N),(N+1)]^2) #sum distances for transformed data)
  
  #compute Blomberg's K
  K.stat = (MSEobs.d/MSE.d)/K.denom
  
  return(K.stat)
}

#!Test.Kmult! compute Blomberg's k for a set of trait data randomization
#so to return the statical significance of the observed K of !x! on !phy!
#Moreover, it computes the ancillary input of !Kmult!

Test.Kmult <- function(x, phy, iter=999, visual=FALSE){
  library(ape) #necessary to call !evol.vcv()!
  
  
  #impose matrix form to the trait data
  x = as.matrix(x)
  
  N<-length(phy$tip.label) #dimension of the problem
  ones<-array(1,N) #an N-vector of 1s
  
  
  C<-vcv.phylo(phy) #the expected vcv structure under bbm
  
  #compute the ancillary C transforms
  eigC <- eigen(C)
  D.mat <-solve(eigC$vectors %*% diag(sqrt(eigC$values)) %*% t(eigC$vectors))
  Ccol <- colSums(solve(C))
  Css <- sum(solve(C))
  
  #compute the denominator in Blomberg's K formula (doesn't involve !x!)
  K.denom<-(sum(diag(C))- N*solve(t(ones)%*%solve(C)%*%ones)) / (N-1)

  #compute Blomberg's K for the observed !x!
  K.obs<-Kmult(x, phy, C, Ccol, Css, D.mat, K.denom)
  
  #compute the statical significance of !K.obs!
  P.val <- 1
  K.val <- rep(0, iter)
  for (i in 1:iter){
    #shuffle the rows of x.r, randomising the observation
    x.r <- as.matrix(x[sample(nrow(x)),])
    rownames(x.r)<-rownames(x) #to avoid it to be reordered back
    
    #compute K for x.r
    K.rand = Kmult(x.r, phy, C, Ccol, Css, D.mat, K.denom)
    
    #compare the random and observed K
    P.val = ifelse(K.rand >= K.obs, P.val+1, P.val)
    K.val[i] = K.rand
  }
  P.val = P.val/(iter + 1)
  K.val[iter + 1] = K.obs
  
  #if !visual = TRUE! plot the K.random vs K.obs histogram
  if(visual){
    hist(K.val, 30, freq = TRUE, col = "gray", xlab = "Phylogenetic Signal")
    arrows(K.obs, 50, K.obs, 5, length = 0.1, lwd = 2)
  }
  
  return(list(phy.signal = K.obs, pvalue = P.val))

}
