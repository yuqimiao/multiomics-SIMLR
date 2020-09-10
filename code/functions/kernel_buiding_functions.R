
library(SIMLR)
library(Matrix)
library("parallel")
library(gplots)
library(igraph)
library("clValid")

# # [self] kernel_calculating funciton:
# using the kernels from different type of data toget the fused matrix directly
# input:
#   x -- m features* n subjects  to calculate kernels 
  # allk_fun -- k neighbors need to be considered
  # sigma_fun -- hyperparameter in kernel building function
  # is.square -- whether square the distance
  # is.normal -- whether use the normal kernel, if not using the RBF kernel
  # is.SP -- whether use the SIMLR proceed
# output:
  # Dkernels -- kernel for the current data type

kernel_calculating = function(x_fun,allk_fun = 30,sigma_fun = 2,is.square = T,is.normal = T,is.SP = T) {
  if(allk_fun<(nrow(x_fun))) {
    # distance calc and sort
    if(is.square){
      Diff_fun = dist2(t(x_fun))^2
    }else{
      Diff_fun = dist2(t(x_fun))
    }

    Diff_sort_fun = t(apply(Diff_fun,MARGIN=2,FUN=sort)) 
    ## sort every column and transpose => ith row is the sorted distance btw ith subject and others
    
    # calc symmetric mu_{ij} for every subject
    TT = apply(Diff_sort_fun[,2:(allk_fun+1)],MARGIN=1,FUN=mean) + .Machine$double.eps
    ## calc mean distancefor the first k neighbors of every subjects(row), length = k
    TT = matrix(data = TT, nrow = length(TT), ncol = 1) ## dim k*1
    Sig = apply(array(0,c(nrow(TT),ncol(TT))),MARGIN=1,FUN=function(x) {x=TT[,1]})
    Sig = Sig + t(Sig)
    Sig = Sig / 2
    Sig_valid = array(0,c(nrow(Sig),ncol(Sig)))
    Sig_valid[which(Sig > .Machine$double.eps,arr.ind=TRUE)] = 1
    Sig = Sig * Sig_valid + .Machine$double.eps
    ## make the symmetric mu_{ij}
    
    if(is.normal){
      W = dnorm(Diff_fun,0,sigma_fun*Sig)
    }else{
      W = dnorm(Diff_fun,0,sigma_fun*Sig)*sigma_fun*Sig
    }
    
    # Output the symmetric kernel by matrix
    D_Kernels = Matrix((W + t(W)) / 2, sparse=TRUE, doDiag=FALSE)
    
    if(is.SP){
      # SIMLR proceed
      #update from SIMLR function (??)
      K = D_Kernels
      k = 1/sqrt(diag(K)+1)
      # diagnal of the matrix is the highest similarity
      G = K * (k %*% t(k))
      ## with diag(K) all ~ 0, k is just n 1s, thus G is K
      G1 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
      G2 = t(G1)
      # Use the average difference btw 2 self similarities and pairwise similarity as kenels
      D_Kernels_tmp = (G1 + G2 - 2*G)/2
      D_Kernels_tmp = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
      D_Kernels = Matrix(D_Kernels_tmp, sparse=TRUE, doDiag=FALSE)
    }
    return(D_Kernels)
  }
}


## dist2 function:
# Goal: calculating the Eucledian square distance btw subjects in matrices
# Input:
#   x -- n1*d1 matrix, every row is a subject
# c -- n2*d2 matrix, every row is a subject
# Procedure:
#   1. If no c, calculating pair-wise subjects difference in x
# 2. The matrix calculation separate 3 calculation in $d_{ij} = \sum_{p=1}^d(x_{ip}-x_{jp})^2$:
#   $\sum_{p=1}^dx_{ip}^2$
#   $\sum_{p=1}^dx_{jp}^2$
#   $\sum_{p=1}^dx_{ip}x_{jp}$
#   Output:
#   distance matric: n1*n2

dist2 = function( x, c = NA ) {
  
  # set the parameters for x
  if(is.na(c)) {
    c = x
  }
  
  # compute the dimension
  n1 = nrow(x)
  d1 = ncol(x)
  n2 = nrow(c)
  d2 = ncol(c)
  if(d1!=d2) {
    stop("Data dimension does not match dimension of centres.")
  }
  
  # compute the distance
  dist = t(rep(1,n2) %*% t(apply(t(x^2),MARGIN=2,FUN=sum))) + 
    (rep(1,n1) %*% t(apply(t(c^2),MARGIN=2,FUN=sum))) - 
    2 * (x%*%t(c))
  
  return(dist)
  
}