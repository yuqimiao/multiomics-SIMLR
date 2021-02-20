dist_kernel = function( Diff, sigma = 1, k_neighb = 30, is.dist = T ) {


  Diff_sort = t(apply(Diff,MARGIN=2,FUN=sort))

  # compute the combined kernels
  m = dim(Diff)[1]
  n = dim(Diff)[2]
  allk = k_neighb
  if(allk<m){
        TT = apply(Diff_sort[,2:(allk+1)],MARGIN=1,FUN=mean) + .Machine$double.eps
        TT = matrix(data = TT, nrow = length(TT), ncol = 1)
        Sig = apply(array(0,c(nrow(TT),ncol(TT))),MARGIN=1,FUN=function(x) {x=TT[,1]})
        Sig = Sig + t(Sig)
        Sig = Sig / 2
        Sig_valid = array(0,c(nrow(Sig),ncol(Sig)))
        Sig_valid[which(Sig > .Machine$double.eps,arr.ind=TRUE)] = 1
        Sig = Sig * Sig_valid + .Machine$double.eps
        W = dnorm(Diff,0,sigma*Sig)
        D_Kernels = Matrix((W + t(W)) / 2, sparse=TRUE, doDiag=FALSE)
  }

  # compute D_Kernels

  if(is.dist){
    K = D_Kernels
    k = 1/sqrt(diag(K)+1)
    G = K * (k %*% t(k))
    G1 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
    G2 = t(G1)
    D_Kernels_tmp = (G1 + G2 - 2*G)/2
    D_Kernels_tmp = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
    D_Kernels = Matrix(D_Kernels_tmp, sparse=TRUE, doDiag=FALSE)
  }

  return(D_Kernels)

}

dist_kernel_RKHS = function(Diff, sigma = 1, k_neighb = 30 ) {


  Diff_sort = t(apply(Diff,MARGIN=2,FUN=sort))

  # compute the combined kernels
  m = dim(Diff)[1]
  n = dim(Diff)[2]
  allk = k_neighb
  if(allk<m){
    TT = apply(Diff_sort[,2:(allk+1)],MARGIN=1,FUN=mean) + .Machine$double.eps
    TT = matrix(data = TT, nrow = length(TT), ncol = 1)
    Sig = apply(array(0,c(nrow(TT),ncol(TT))),MARGIN=1,FUN=function(x) {x=TT[,1]})
    Sig = Sig + t(Sig)
    Sig = Sig / 2
    Sig_valid = array(0,c(nrow(Sig),ncol(Sig)))
    Sig_valid[which(Sig > .Machine$double.eps,arr.ind=TRUE)] = 1
    Sig = Sig * Sig_valid + .Machine$double.eps
    W = dnorm(Diff,0,sigma*Sig)*sigma*Sig
    D_Kernels = Matrix((W + t(W)) / 2, sparse=TRUE, doDiag=FALSE)
  }

  # compute D_Kernels

  K = D_Kernels
  k = 1/sqrt(diag(K)+1)
  G = K * (k %*% t(k))
  G1 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
  G2 = t(G1)
  D_Kernels_tmp = (G1 + G2 - 2*G)/2
  D_Kernels_tmp = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
  D_Kernels = Matrix(D_Kernels_tmp, sparse=TRUE, doDiag=FALSE)


  return(D_Kernels)

}
