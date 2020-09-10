
library(SIMLR)
library(Matrix)
library("parallel")
library(gplots)
library(igraph)
library("clValid")

source("./code/R/SIMLR.R")
source("./code/R/compute.multiple.kernel.R")
source("./code/R/network.diffusion.R")
source("./code/R/utils.simlr.R")
source("./code/R/tsne.R")

dyn.load("./code/R/projsplx_R.so")
# normalize2:
# normalize according to columns
#   input:
#     x -- n subjects * m features to calculate kernels

normalize2 = function(x) {
  x = as.matrix(x);
  mean = apply(x, 2, mean)
  sd = apply(x, 2, sd)
  sd[sd==0] = 1
  xNorm = t((t(x) - mean) / sd)
  return(abs(xNorm))
}

# kernel_calculating funciton:
# using the kernels from different type of data toget the fused matrix directly
#   input:
#     x -- n subjects * m features to calculate kernels
#     allk_fun -- k neighbors need to be considered
#     sigma_fun -- hyperparameter in kernel building function
#   output:
#       Dkernels -- kernel for the current data type

# kernel_calculating = function(x_fun,allk_fun = 30,sigma_fun = 2) {
#   if(allk_fun<(nrow(x_fun)-1)) {
#     # distance calc and sort
#     Diff_fun = dist2(t(x_fun))^2
#     Diff_sort_fun = t(apply(Diff_fun,MARGIN=2,FUN=sort))
#     ## sort every column and transpose => ith row is the sorted distance btw ith subject and others
#
#     # calc symmetric mu_{ij} for every subject
#     TT = apply(Diff_sort_fun[,2:(allk_fun+1)],MARGIN=1,FUN=mean) + .Machine$double.eps
#     ## calc mean distancefor the first k neighbors of every subjects(row), length = k
#     TT = matrix(data = TT, nrow = length(TT), ncol = 1) ## dim k*1
#     Sig = apply(array(0,c(nrow(TT),ncol(TT))),MARGIN=1,FUN=function(x) {x=TT[,1]})
#     Sig = Sig + t(Sig)
#     Sig = Sig / 2
#     Sig_valid = array(0,c(nrow(Sig),ncol(Sig)))
#     Sig_valid[which(Sig > .Machine$double.eps,arr.ind=TRUE)] = 1
#     Sig = Sig * Sig_valid + .Machine$double.eps
#     ## make the symmetric mu_{ij}
#
#     # calc the kernel using normal diensity
#     W = dnorm(Diff_fun,0,sigma_fun*Sig)*sigma_fun*Sig
#
#     # Output the symmetric kernel by matrix
#     D_Kernels = Matrix((W + t(W)) / 2, sparse=TRUE, doDiag=FALSE)
#
#     # update from SIMLR function (??)
#     K = D_Kernels
#     k = 1/sqrt(diag(K)+1)
#     # diagnal of the matrix is the highest similarity
#     G = K * (k %*% t(k))
#     ## with diag(K) all ~ 0, k is just n 1s, thus G is K
#     G1 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
#     G2 = t(G1)
#     # ??
#     D_Kernels_tmp = (G1 + G2 - 2*G)/2
#     D_Kernels_tmp = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
#     D_Kernels = Matrix(D_Kernels_tmp, sparse=TRUE, doDiag=FALSE)
#     return(D_Kernels)
#   }
# }


"SIMLR_multi" = function(D_Kernels, c, no.dim = NA, k = 10, if.impute = FALSE, normalize = FALSE, cores.ratio = 1 ) {

  #set any required parameter to the defaults
  if(is.na(no.dim)) {
    no.dim = c
  }

  # # check the if.impute parameter
  # if(if.impute == TRUE) {
  #   X = t(X)
  #   X_zeros = which(X==0,arr.ind=TRUE)
  #   if(length(X_zeros)>0) {
  #     R_zeros = as.vector(X_zeros[,"row"])
  #     C_zeros = as.vector(X_zeros[,"col"])
  #     ind = (C_zeros - 1) * nrow(X) + R_zeros
  #     X[ind] = as.vector(colMeans(X))[C_zeros]
  #   }
  #   X = t(X)
  # }
  #
  # # check the normalize parameter
  # if(normalize == TRUE) {
  #   X = t(X)
  #   X = X - min(as.vector(X))
  #   X = X / max(as.vector(X))
  #   C_mean = as.vector(colMeans(X))
  #   X = apply(X,MARGIN=1,FUN=function(x) return(x-C_mean))
  # }
  #
  # # start the clock to measure the execution time
  ptm = proc.time()
  #
  # set some parameters
  NITER = 30
  num = ncol(D_Kernels[[1]])
  r = -1
  beta = 0.8
  #
  # cat("Computing the multiple Kernels.\n")
  #
  # # compute the kernels
  # D_Kernels = multiple.kernel(t(X),cores.ratio)
  #
  # set up some parameters
  alphaK = 1 / rep(length(D_Kernels),length(D_Kernels))
  distX = array(0,c(dim(D_Kernels[[1]])[1],dim(D_Kernels[[1]])[2]))
  for (i in 1:length(D_Kernels)) {
    distX = distX + D_Kernels[[i]]
  }
  distX = distX / length(D_Kernels)

  # sort distX for rows
  res = apply(distX,MARGIN=1,FUN=function(x) return(sort(x,index.return = TRUE)))
  distX1 = array(0,c(nrow(distX),ncol(distX)))
  idx = array(0,c(nrow(distX),ncol(distX)))
  for(i in 1:nrow(distX)) {
    distX1[i,] = res[[i]]$x
    idx[i,] = res[[i]]$ix
  }

  A = array(0,c(num,num))
  di = distX1[,2:(k+2)]
  rr = 0.5 * (k * di[,k+1] - apply(di[,1:k],MARGIN=1,FUN=sum))
  id = idx[,2:(k+2)]

  numerator = (apply(array(0,c(length(di[,k+1]),dim(di)[2])),MARGIN=2,FUN=function(x) {x=di[,k+1]}) - di)
  temp = (k*di[,k+1] - apply(di[,1:k],MARGIN=1,FUN=sum) + .Machine$double.eps)
  denominator = apply(array(0,c(length(temp),dim(di)[2])),MARGIN=2,FUN=function(x) {x=temp})
  temp = numerator / denominator
  a = apply(array(0,c(length(t(1:num)),dim(di)[2])),MARGIN=2,FUN=function(x) {x=1:num})
  A[cbind(as.vector(a),as.vector(id))] = as.vector(temp)
  if(r<=0) {
    r = mean(rr)
  }
  lambda = max(mean(rr),0)
  A[is.nan(A)] = 0
  A0 = (A + t(A)) / 2
  S0 = max(max(distX)) - distX

  # print(cat("Performing network diffiusion.\n"))
  #
  # # perform network diffiusion
  # S0 = network.diffusion(S0,k)

  # compute dn
  S0 = dn(S0,'ave')
  S = S0
  D0 = diag(apply(S,MARGIN=2,FUN=sum))
  L0 = D0 - S

  eig1_res = eig1(L0,c,0)
  F_eig1 = eig1_res$eigvec
  temp_eig1 = eig1_res$eigval
  evs_eig1 = eig1_res$eigval_full

  # perform the iterative procedure NITER times
  converge = vector()
  for(iter in 1:NITER) {

    cat("Iteration: ",iter,"\n")

    distf = L2_distance_1(t(F_eig1),t(F_eig1))
    A = array(0,c(num,num))
    b = idx[,2:dim(idx)[2]]
    a = apply(array(0,c(num,ncol(b))),MARGIN=2,FUN=function(x){ x = 1:num })
    inda = cbind(as.vector(a),as.vector(b))
    ad = (distX[inda]+lambda*distf[inda])/2/r
    dim(ad) = c(num,ncol(b))

    # call the c function for the optimization
    c_input = -t(ad)
    c_output = t(ad)
    ad = t(.Call("projsplx_R",c_input,c_output))

    A[inda] = as.vector(ad)
    A[is.nan(A)] = 0
    A = (A + t(A)) / 2
    S = (1 - beta) * S + beta * A
    S = as.matrix(S)
    #S = network.diffusion(S,k)
    D = diag(apply(S,MARGIN=2,FUN=sum))
    L = D - S
    F_old = F_eig1
    eig1_res = eig1(L,c,0)
    F_eig1 = eig1_res$eigvec
    temp_eig1 = eig1_res$eigval
    ev_eig1 = eig1_res$eigval_full
    evs_eig1 = cbind(evs_eig1,ev_eig1)
    DD = vector()
    for (i in 1:length(D_Kernels)) {
      temp = (.Machine$double.eps+D_Kernels[[i]]) * (S+.Machine$double.eps)
      DD[i] = mean(apply(temp,MARGIN=2,FUN=sum))
    }
    alphaK0 = umkl(DD)
    alphaK0 = alphaK0 / sum(alphaK0)
    alphaK = (1-beta) * alphaK + beta * alphaK0
    alphaK = alphaK / sum(alphaK)
    fn1 = sum(ev_eig1[1:c])
    fn2 = sum(ev_eig1[1:(c+1)])
    converge[iter] = fn2 - fn1
    if (iter<10) {
      if (ev_eig1[length(ev_eig1)] > 0.000001) {
        lambda = 1.5 * lambda
        r = r / 1.01
      }
    }
    else {
      if(converge[iter]>converge[iter-1]) {
        S = S_old
        if(converge[iter-1] > 0.2) {
          warning('Maybe you should set a larger value of c.')
        }
        break
      }
    }
    S_old = S

    # compute Kbeta
    distX = D_Kernels[[1]] * alphaK[1]
    for (i in 2:length(D_Kernels)) {
      distX = distX + as.matrix(D_Kernels[[i]]) * alphaK[i]
    }

    # sort distX for rows
    res = apply(distX,MARGIN=1,FUN=function(x) return(sort(x,index.return = TRUE)))
    distX1 = array(0,c(nrow(distX),ncol(distX)))
    idx = array(0,c(nrow(distX),ncol(distX)))
    for(i in 1:nrow(distX)) {
      distX1[i,] = res[[i]]$x
      idx[i,] = res[[i]]$ix
    }

  }
  LF = F_eig1
  D = diag(apply(S,MARGIN=2,FUN=sum))
  L = D - S

  # compute the eigenvalues and eigenvectors of P
  eigen_L = eigen(L)
  U = eigen_L$vectors
  D = eigen_L$values

  if (length(no.dim)==1) {
    print(ncol(U))
    U_index = seq(ncol(U),(ncol(U)-no.dim+1))
    F_last = tsne(S,k=no.dim,initial_config=U[,U_index])
  }
  else {
    F_last = list()
    for (i in 1:length(no.dim)) {
      U_index = seq(ncol(U),(ncol(U)-no.dim[i]+1))
      F_last[i] = tsne(S,k=no.dim[i],initial_config=U[,U_index])
    }
  }

  # compute the execution time
  execution.time = proc.time() - ptm

  cat("Performing Kmeans.\n")
  y = kmeans(F_last,c,nstart=200)

  ydata = tsne(S)

  # create the structure with the results
  results = list()
  results[["y"]] = y
  results[["S"]] = S
  results[["F"]] = F_last
  results[["ydata"]] = ydata
  results[["alphaK"]] = alphaK
  results[["execution.time"]] = execution.time
  results[["converge"]] = converge
  results[["LF"]] = LF

  return(results)

}
# using SNF as initial
"SIMLR_multi_SNF" = function(D_Kernels, c, SNF_init,no.dim = NA, k = 10, if.impute = FALSE, normalize = FALSE, cores.ratio = 1 ) {

  #set any required parameter to the defaults
  if(is.na(no.dim)) {
    no.dim = c
  }
  ptm = proc.time()
  #
  # set some parameters
  NITER = 30
  num = ncol(D_Kernels[[1]])
  r = -1
  beta = 0.8

  ## initial weight and weighted kernels
  alphaK = 1 / rep(length(D_Kernels),length(D_Kernels))
  distX = array(0,c(dim(D_Kernels[[1]])[1],dim(D_Kernels[[1]])[2]))
  for (i in 1:length(D_Kernels)) {
    distX = distX + D_Kernels[[i]]
  }
  distX = distX / length(D_Kernels)

  # sort distX for rows
  res = apply(distX,MARGIN=1,FUN=function(x) return(sort(x,index.return = TRUE)))
  distX1 = array(0,c(nrow(distX),ncol(distX)))
  idx = array(0,c(nrow(distX),ncol(distX)))
  for(i in 1:nrow(distX)) {
    distX1[i,] = res[[i]]$x
    idx[i,] = res[[i]]$ix
  }

  A = array(0,c(num,num))
  di = distX1[,2:(k+2)]
  rr = 0.5 * (k * di[,k+1] - apply(di[,1:k],MARGIN=1,FUN=sum))
  id = idx[,2:(k+2)]

  numerator = (apply(array(0,c(length(di[,k+1]),dim(di)[2])),MARGIN=2,FUN=function(x) {x=di[,k+1]}) - di)
  temp = (k*di[,k+1] - apply(di[,1:k],MARGIN=1,FUN=sum) + .Machine$double.eps)
  denominator = apply(array(0,c(length(temp),dim(di)[2])),MARGIN=2,FUN=function(x) {x=temp})
  temp = numerator / denominator
  a = apply(array(0,c(length(t(1:num)),dim(di)[2])),MARGIN=2,FUN=function(x) {x=1:num})
  A[cbind(as.vector(a),as.vector(id))] = as.vector(temp)
  if(r<=0) {
    r = mean(rr)
  }
  lambda = max(mean(rr),0)
  A[is.nan(A)] = 0
  A0 = (A + t(A)) / 2
  S0 = SNF_init

  # print(cat("Performing network diffiusion.\n"))
  #
  # # perform network diffiusion
  # S0 = network.diffusion(S0,k)

  # compute dn
  S0 = dn(S0,'ave')
  S = S0
  D0 = diag(apply(S,MARGIN=2,FUN=sum))
  L0 = D0 - S

  eig1_res = eig1(L0,c,0)
  F_eig1 = eig1_res$eigvec
  temp_eig1 = eig1_res$eigval
  evs_eig1 = eig1_res$eigval_full

  # perform the iterative procedure NITER times
  converge = vector()
  for(iter in 1:NITER) {

    cat("Iteration: ",iter,"\n")

    distf = L2_distance_1(t(F_eig1),t(F_eig1))
    A = array(0,c(num,num))
    b = idx[,2:dim(idx)[2]]
    a = apply(array(0,c(num,ncol(b))),MARGIN=2,FUN=function(x){ x = 1:num })
    inda = cbind(as.vector(a),as.vector(b))
    ad = (distX[inda]+lambda*distf[inda])/2/r
    dim(ad) = c(num,ncol(b))

    # call the c function for the optimization
    c_input = -t(ad)
    c_output = t(ad)
    ad = t(.Call("projsplx_R",c_input,c_output))

    A[inda] = as.vector(ad)
    A[is.nan(A)] = 0
    A = (A + t(A)) / 2
    S = (1 - beta) * S + beta * A
    S = as.matrix(S)
    #S = network.diffusion(S,k)
    D = diag(apply(S,MARGIN=2,FUN=sum))
    L = D - S
    F_old = F_eig1
    eig1_res = eig1(L,c,0)
    F_eig1 = eig1_res$eigvec
    temp_eig1 = eig1_res$eigval
    ev_eig1 = eig1_res$eigval_full
    evs_eig1 = cbind(evs_eig1,ev_eig1)
    DD = vector()
    for (i in 1:length(D_Kernels)) {
      temp = (.Machine$double.eps+D_Kernels[[i]]) * (S+.Machine$double.eps)
      DD[i] = mean(apply(temp,MARGIN=2,FUN=sum))
    }
    alphaK0 = umkl(DD)
    alphaK0 = alphaK0 / sum(alphaK0)
    alphaK = (1-beta) * alphaK + beta * alphaK0
    alphaK = alphaK / sum(alphaK)
    fn1 = sum(ev_eig1[1:c])
    fn2 = sum(ev_eig1[1:(c+1)])
    converge[iter] = fn2 - fn1
    if (iter<10) {
      if (ev_eig1[length(ev_eig1)] > 0.000001) {
        lambda = 1.5 * lambda
        r = r / 1.01
      }
    }
    else {
      if(converge[iter]>converge[iter-1]) {
        S = S_old
        if(converge[iter-1] > 0.2) {
          warning('Maybe you should set a larger value of c.')
        }
        break
      }
    }
    S_old = S

    # compute Kbeta
    distX = D_Kernels[[1]] * alphaK[1]
    for (i in 2:length(D_Kernels)) {
      distX = distX + as.matrix(D_Kernels[[i]]) * alphaK[i]
    }

    # sort distX for rows
    res = apply(distX,MARGIN=1,FUN=function(x) return(sort(x,index.return = TRUE)))
    distX1 = array(0,c(nrow(distX),ncol(distX)))
    idx = array(0,c(nrow(distX),ncol(distX)))
    for(i in 1:nrow(distX)) {
      distX1[i,] = res[[i]]$x
      idx[i,] = res[[i]]$ix
    }

  }
  LF = F_eig1
  D = diag(apply(S,MARGIN=2,FUN=sum))
  L = D - S

  # compute the eigenvalues and eigenvectors of P
  eigen_L = eigen(L)
  U = eigen_L$vectors
  D = eigen_L$values

  if (length(no.dim)==1) {
    print(ncol(U))
    U_index = seq(ncol(U),(ncol(U)-no.dim+1))
    F_last = tsne(S,k=no.dim,initial_config=U[,U_index])
  }
  else {
    F_last = list()
    for (i in 1:length(no.dim)) {
      U_index = seq(ncol(U),(ncol(U)-no.dim[i]+1))
      F_last[i] = tsne(S,k=no.dim[i],initial_config=U[,U_index])
    }
  }

  # compute the execution time
  execution.time = proc.time() - ptm

  cat("Performing Kmeans.\n")
  y = kmeans(F_last,c,nstart=200)

  ydata = tsne(S)

  # create the structure with the results
  results = list()
  results[["y"]] = y
  results[["S"]] = S
  results[["F"]] = F_last
  results[["ydata"]] = ydata
  results[["alphaK"]] = alphaK
  results[["execution.time"]] = execution.time
  results[["converge"]] = converge
  results[["LF"]] = LF

  return(results)

}

# cw for constraint weights
"SIMLR_multi_SNF_cw" = function(D_Kernels, c, w1=0.5,w2 = 0.5,SNF_init,no.dim = NA, k = 10, if.impute = FALSE, normalize = FALSE, cores.ratio = 1 ) {

  #set any required parameter to the defaults
  if(is.na(no.dim)) {
    no.dim = c
  }
  ptm = proc.time()
  #
  # set some parameters
  NITER = 30
  num = ncol(D_Kernels[[1]])
  r = -1
  beta = 0.8

  ## initial weight and weighted kernels
  alphaK = 1 / rep(length(D_Kernels),length(D_Kernels))
  distX = array(0,c(dim(D_Kernels[[1]])[1],dim(D_Kernels[[1]])[2]))
  for (i in 1:length(D_Kernels)) {
    distX = distX + D_Kernels[[i]]
  }
  distX = distX / length(D_Kernels)

  # sort distX for rows
  res = apply(distX,MARGIN=1,FUN=function(x) return(sort(x,index.return = TRUE)))
  distX1 = array(0,c(nrow(distX),ncol(distX)))
  idx = array(0,c(nrow(distX),ncol(distX)))
  for(i in 1:nrow(distX)) {
    distX1[i,] = res[[i]]$x
    idx[i,] = res[[i]]$ix
  }

  A = array(0,c(num,num))
  di = distX1[,2:(k+2)]
  rr = 0.5 * (k * di[,k+1] - apply(di[,1:k],MARGIN=1,FUN=sum))
  id = idx[,2:(k+2)]

  numerator = (apply(array(0,c(length(di[,k+1]),dim(di)[2])),MARGIN=2,FUN=function(x) {x=di[,k+1]}) - di)
  temp = (k*di[,k+1] - apply(di[,1:k],MARGIN=1,FUN=sum) + .Machine$double.eps)
  denominator = apply(array(0,c(length(temp),dim(di)[2])),MARGIN=2,FUN=function(x) {x=temp})
  temp = numerator / denominator
  a = apply(array(0,c(length(t(1:num)),dim(di)[2])),MARGIN=2,FUN=function(x) {x=1:num})
  A[cbind(as.vector(a),as.vector(id))] = as.vector(temp)
  if(r<=0) {
    r = mean(rr)
  }
  lambda = max(mean(rr),0)
  A[is.nan(A)] = 0
  A0 = (A + t(A)) / 2
  S0 = SNF_init

  # print(cat("Performing network diffiusion.\n"))
  #
  # # perform network diffiusion
  # S0 = network.diffusion(S0,k)

  # compute dn
  S0 = dn(S0,'ave')
  S = S0
  D0 = diag(apply(S,MARGIN=2,FUN=sum))
  L0 = D0 - S

  eig1_res = eig1(L0,c,0)
  F_eig1 = eig1_res$eigvec
  temp_eig1 = eig1_res$eigval
  evs_eig1 = eig1_res$eigval_full

  # perform the iterative procedure NITER times
  converge = vector()
  for(iter in 1:NITER) {

    cat("Iteration: ",iter,"\n")

    distf = L2_distance_1(t(F_eig1),t(F_eig1))
    A = array(0,c(num,num))
    b = idx[,2:dim(idx)[2]]
    a = apply(array(0,c(num,ncol(b))),MARGIN=2,FUN=function(x){ x = 1:num })
    inda = cbind(as.vector(a),as.vector(b))
    ad = (distX[inda]+lambda*distf[inda])/2/r
    dim(ad) = c(num,ncol(b))

    # call the c function for the optimization
    c_input = -t(ad)
    c_output = t(ad)
    ad = t(.Call("projsplx_R",c_input,c_output))

    A[inda] = as.vector(ad)
    A[is.nan(A)] = 0
    A = (A + t(A)) / 2
    S = (1 - beta) * S + beta * A
    S = as.matrix(S)
    #S = network.diffusion(S,k)
    D = diag(apply(S,MARGIN=2,FUN=sum))
    L = D - S
    F_old = F_eig1
    eig1_res = eig1(L,c,0)
    F_eig1 = eig1_res$eigvec
    temp_eig1 = eig1_res$eigval
    ev_eig1 = eig1_res$eigval_full
    evs_eig1 = cbind(evs_eig1,ev_eig1)
    DD = vector()
    for (i in 1:length(D_Kernels)) {
      temp = (.Machine$double.eps+D_Kernels[[i]]) * (S+.Machine$double.eps)
      DD[i] = mean(apply(temp,MARGIN=2,FUN=sum))
    }
    alphaK0 = umkl(DD)
    alphaK0 = alphaK0 / sum(alphaK0)
    alphaK = (1-beta) * alphaK + beta * alphaK0
    ############################################################
    if(is.na(w1)){
      alphaK = alphaK / sum(alphaK)
    }else{
      alphaK[1:55] = (alphaK[1:55]/sum(alphaK[1:55]))*w1
      alphaK[56:110] = (alphaK[56:110]/sum(alphaK[56:110]))*w2
    }
    #############################################################
    fn1 = sum(ev_eig1[1:c])
    fn2 = sum(ev_eig1[1:(c+1)])
    converge[iter] = fn2 - fn1
    if (iter<10) {
      if (ev_eig1[length(ev_eig1)] > 0.000001) {
        lambda = 1.5 * lambda
        r = r / 1.01
      }
    }
    else {
      if(converge[iter]>converge[iter-1]) {
        S = S_old
        if(converge[iter-1] > 0.2) {
          warning('Maybe you should set a larger value of c.')
        }
        break
      }
    }
    S_old = S

    # compute Kbeta
    distX = D_Kernels[[1]] * alphaK[1]
    for (i in 2:length(D_Kernels)) {
      distX = distX + as.matrix(D_Kernels[[i]]) * alphaK[i]
    }

    # sort distX for rows
    res = apply(distX,MARGIN=1,FUN=function(x) return(sort(x,index.return = TRUE)))
    distX1 = array(0,c(nrow(distX),ncol(distX)))
    idx = array(0,c(nrow(distX),ncol(distX)))
    for(i in 1:nrow(distX)) {
      distX1[i,] = res[[i]]$x
      idx[i,] = res[[i]]$ix
    }

  }
  LF = F_eig1
  D = diag(apply(S,MARGIN=2,FUN=sum))
  L = D - S

  # compute the eigenvalues and eigenvectors of P
  eigen_L = eigen(L)
  U = eigen_L$vectors
  D = eigen_L$values

  if (length(no.dim)==1) {
    print(ncol(U))
    U_index = seq(ncol(U),(ncol(U)-no.dim+1))
    F_last = tsne(S,k=no.dim,initial_config=U[,U_index])
  }
  else {
    F_last = list()
    for (i in 1:length(no.dim)) {
      U_index = seq(ncol(U),(ncol(U)-no.dim[i]+1))
      F_last[i] = tsne(S,k=no.dim[i],initial_config=U[,U_index])
    }
  }

  # compute the execution time
  execution.time = proc.time() - ptm

  cat("Performing Kmeans.\n")
  y = kmeans(F_last,c,nstart=200)

  ydata = tsne(S)

  # create the structure with the results
  results = list()
  results[["y"]] = y
  results[["S"]] = S
  results[["F"]] = F_last
  results[["ydata"]] = ydata
  results[["alphaK"]] = alphaK
  results[["execution.time"]] = execution.time
  results[["converge"]] = converge
  results[["LF"]] = LF

  return(results)

}

##icw stands for inner constraint weight_sum
"SIMLR_multi_SNF_icw" = function(D_Kernels, c, w1 = 0.5,w2 = 0.5, u = 20,SNF_init,no.dim = NA, k = 10, if.impute = FALSE, normalize = FALSE, cores.ratio = 1 ) {

  #set any required parameter to the defaults
  if(is.na(no.dim)) {
    no.dim = c
  }
  ptm = proc.time()
  #
  # set some parameters
  NITER = 30
  num = ncol(D_Kernels[[1]])
  r = -1
  beta = 0.8

  ## initial weight and weighted kernels
  alphaK = 1 / rep(length(D_Kernels),length(D_Kernels))
  distX = array(0,c(dim(D_Kernels[[1]])[1],dim(D_Kernels[[1]])[2]))
  for (i in 1:length(D_Kernels)) {
    distX = distX + D_Kernels[[i]]
  }
  distX = distX / length(D_Kernels)

  # sort distX for rows
  res = apply(distX,MARGIN=1,FUN=function(x) return(sort(x,index.return = TRUE)))
  distX1 = array(0,c(nrow(distX),ncol(distX)))
  idx = array(0,c(nrow(distX),ncol(distX)))
  for(i in 1:nrow(distX)) {
    distX1[i,] = res[[i]]$x
    idx[i,] = res[[i]]$ix
  }

  A = array(0,c(num,num))
  di = distX1[,2:(k+2)]
  rr = 0.5 * (k * di[,k+1] - apply(di[,1:k],MARGIN=1,FUN=sum))
  id = idx[,2:(k+2)]

  numerator = (apply(array(0,c(length(di[,k+1]),dim(di)[2])),MARGIN=2,FUN=function(x) {x=di[,k+1]}) - di)
  temp = (k*di[,k+1] - apply(di[,1:k],MARGIN=1,FUN=sum) + .Machine$double.eps)
  denominator = apply(array(0,c(length(temp),dim(di)[2])),MARGIN=2,FUN=function(x) {x=temp})
  temp = numerator / denominator
  a = apply(array(0,c(length(t(1:num)),dim(di)[2])),MARGIN=2,FUN=function(x) {x=1:num})
  A[cbind(as.vector(a),as.vector(id))] = as.vector(temp)
  if(r<=0) {
    r = mean(rr)
  }
  lambda = max(mean(rr),0)
  A[is.nan(A)] = 0
  A0 = (A + t(A)) / 2
  S0 = SNF_init

  # print(cat("Performing network diffiusion.\n"))
  #
  # # perform network diffiusion
  # S0 = network.diffusion(S0,k)

  # compute dn
  S0 = dn(S0,'ave')
  S = S0
  D0 = diag(apply(S,MARGIN=2,FUN=sum))
  L0 = D0 - S

  eig1_res = eig1(L0,c,0)
  F_eig1 = eig1_res$eigvec
  temp_eig1 = eig1_res$eigval
  evs_eig1 = eig1_res$eigval_full

  # perform the iterative procedure NITER times
  converge = vector()
  for(iter in 1:NITER) {

    cat("Iteration: ",iter,"\n")

    distf = L2_distance_1(t(F_eig1),t(F_eig1))
    A = array(0,c(num,num))
    b = idx[,2:dim(idx)[2]]
    a = apply(array(0,c(num,ncol(b))),MARGIN=2,FUN=function(x){ x = 1:num })
    inda = cbind(as.vector(a),as.vector(b))
    ad = (distX[inda]+lambda*distf[inda])/2/r
    dim(ad) = c(num,ncol(b))

    # call the c function for the optimization
    c_input = -t(ad)
    c_output = t(ad)
    ad = t(.Call("projsplx_R",c_input,c_output))

    A[inda] = as.vector(ad)
    A[is.nan(A)] = 0
    A = (A + t(A)) / 2
    S = (1 - beta) * S + beta * A
    S = as.matrix(S)
    #S = network.diffusion(S,k)
    D = diag(apply(S,MARGIN=2,FUN=sum))
    L = D - S
    F_old = F_eig1
    eig1_res = eig1(L,c,0)
    F_eig1 = eig1_res$eigvec
    temp_eig1 = eig1_res$eigval
    ev_eig1 = eig1_res$eigval_full
    evs_eig1 = cbind(evs_eig1,ev_eig1)
    DD = vector()
    for (i in 1:length(D_Kernels)) {
      temp = (.Machine$double.eps+D_Kernels[[i]]) * (S+.Machine$double.eps)
      DD[i] = mean(apply(temp,MARGIN=2,FUN=sum))
    }
    alphaK0 = umkl_test(DD,w1,w2,u)
    alphaK0 = alphaK0 / sum(alphaK0)
    alphaK = (1-beta) * alphaK + beta * alphaK0
    alphaK = alphaK / sum(alphaK)
    fn1 = sum(ev_eig1[1:c])
    fn2 = sum(ev_eig1[1:(c+1)])
    converge[iter] = fn2 - fn1
    if (iter<10) {
      if (ev_eig1[length(ev_eig1)] > 0.000001) {
        lambda = 1.5 * lambda
        r = r / 1.01
      }
    }
    else {
      if(converge[iter]>converge[iter-1]) {
        S = S_old
        if(converge[iter-1] > 0.2) {
          warning('Maybe you should set a larger value of c.')
        }
        break
      }
    }
    S_old = S

    # compute Kbeta
    distX = D_Kernels[[1]] * alphaK[1]
    for (i in 2:length(D_Kernels)) {
      distX = distX + as.matrix(D_Kernels[[i]]) * alphaK[i]
    }

    # sort distX for rows
    res = apply(distX,MARGIN=1,FUN=function(x) return(sort(x,index.return = TRUE)))
    distX1 = array(0,c(nrow(distX),ncol(distX)))
    idx = array(0,c(nrow(distX),ncol(distX)))
    for(i in 1:nrow(distX)) {
      distX1[i,] = res[[i]]$x
      idx[i,] = res[[i]]$ix
    }

  }
  LF = F_eig1
  D = diag(apply(S,MARGIN=2,FUN=sum))
  L = D - S

  # compute the eigenvalues and eigenvectors of P
  eigen_L = eigen(L)
  U = eigen_L$vectors
  D = eigen_L$values

  if (length(no.dim)==1) {
    print(ncol(U))
    U_index = seq(ncol(U),(ncol(U)-no.dim+1))
    F_last = tsne(S,k=no.dim,initial_config=U[,U_index])
  }
  else {
    F_last = list()
    for (i in 1:length(no.dim)) {
      U_index = seq(ncol(U),(ncol(U)-no.dim[i]+1))
      F_last[i] = tsne(S,k=no.dim[i],initial_config=U[,U_index])
    }
  }

  # compute the execution time
  execution.time = proc.time() - ptm

  cat("Performing Kmeans.\n")
  y = kmeans(F_last,c,nstart=200)

  ydata = tsne(S)

  # create the structure with the results
  results = list()
  results[["y"]] = y
  results[["S"]] = S
  results[["F"]] = F_last
  results[["ydata"]] = ydata
  results[["alphaK"]] = alphaK
  results[["execution.time"]] = execution.time
  results[["converge"]] = converge
  results[["LF"]] = LF

  return(results)

}

umkl_test = function( D, beta = NA , w1, w2, u = 20) {

  # set some parameters
  if(is.na(beta)) {
    beta = 1 / length(D)
  }
  tol = 1e-4
  u = u
  logU = log(u)

  # compute Hbeta
  res_hbeta = Hbeta_test(D, beta,w1,w2)
  H = res_hbeta$H
  thisP = res_hbeta$P

  betamin = -Inf
  betamax = Inf
  # evaluate whether the perplexity is within tolerance
  Hdiff = H - logU
  tries = 0
  while (abs(Hdiff) > tol && tries < 30) {
    #if not, increase or decrease precision
    if (Hdiff > 0) {
      betamin = beta
      if(abs(betamax)==Inf) {
        beta = beta * 2
      }
      else {
        beta = (beta + betamax) / 2
      }
    }
    else {
      betamax = beta
      if(abs(betamin)==Inf) {
        beta = beta * 2
      }
      else {
        beta = (beta + betamin) / 2
      }
    }
    # compute the new values
    res_hbeta = Hbeta(D, beta)
    H = res_hbeta$H
    thisP = res_hbeta$P
    Hdiff = H - logU
    tries = tries + 1
  }

  return(thisP)

}
Hbeta_test = function( D, beta,w1,w2 ) {

  D = (D - min(D)) / (max(D) - min(D) + .Machine$double.eps)
  P = exp(-D * beta)
  sumP = sum(P)
  H = log(sumP) + beta * sum(D * P) / sumP
  P[1:55] =( P / sum(P[1:55]))*w1
  P[56:110] = (P / sum(P[56:110]))*w2

  return(list(H=H,P=P))

}
