setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")

library(dplyr)
library(parallel)
library(clValid)
library(SNFtool)
library(igraph)
library(Matrix)


##########
"dist_2" = function( x, c = NA ) {

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
standardNormalization = function (x)
{
  x = as.matrix(x)
  mean = apply(x, 2, mean)
  sd = apply(x, 2, sd)
  sd[sd == 0] = 1
  xNorm = t((t(x) - mean)/sd)
  return(xNorm)
}

kernel_calculating = function(x_fun,allk_fun = 30,sigma_fun = 2,is.square = T,is.normal = T,is.SP = T) {
  x_fun = t(standardNormalization(t(x_fun)))
  if(allk_fun<(nrow(x_fun))) {
    # distance calc and sort
    if(is.square){
      Diff_fun = dist_2(t(x_fun))^2
    }else{
      Diff_fun = dist_2(t(x_fun))
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
    # if is.SP, using dij = kii+kjj-2kijto calculate distance based on kernel value
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

normalize_kernel <- function(X) {
  ## This function normalize the kernel functions
  row.sum.mdiag <- rowSums(X) - diag(X)
  row.sum.mdiag[row.sum.mdiag == 0] <- 1
  X <- X/(2 * (row.sum.mdiag))
  diag(X) <- 0.5
  return(X)
}

multiple.kernel2 = function( x, cores.ratio = 1 ) {

  # set the parameters
  kernel.type = list()
  kernel.type[1] = list("poly")
  kernel.params = list()
  kernel.params[1] = list(0)

  # compute some parameters from the kernels
  N = dim(x)[1]
  KK = 0
  sigma = seq(2,1,-0.5)

  # compute and sort Diff
  Diff = dist2(x)^2
  Diff_sort = t(apply(Diff,MARGIN=2,FUN=sort))

  # compute the combined kernels
  m = dim(Diff)[1]
  n = dim(Diff)[2]
  allk = 30

  ## calculate the mean k-neignbor distance for every subject
  if(allk>(nrow(x)-1)){
    allk = nrow(x)-2
  }

  TT = apply(Diff_sort[,2:(allk+1)],MARGIN=1,FUN=mean) + .Machine$double.eps
  TT = matrix(data = TT, nrow = length(TT), ncol = 1)
  Sig = apply(array(0,c(nrow(TT),ncol(TT))),MARGIN=1,FUN=function(x) {x=TT[,1]})
  Sig = Sig + t(Sig)
  Sig = Sig / 2
  Sig_valid = array(0,c(nrow(Sig),ncol(Sig)))
  Sig_valid[which(Sig > .Machine$double.eps,arr.ind=TRUE)] = 1
  Sig = Sig * Sig_valid + .Machine$double.eps

  # setup a parallelized estimation of the kernels
  cores = as.integer(cores.ratio * (detectCores() - 1))
  if (cores < 1 || is.na(cores) || is.null(cores)) {
    cores = 1
  }

  cl = makeCluster(cores)

  clusterEvalQ(cl, {library(Matrix)})

  D_Kernels = list()
  D_Kernels = unlist(parLapply(cl,1:length(sigma),fun=function(l,x_fun=x,Diff_sort_fun=Diff_sort,allk_fun=allk,
                                                               Diff_fun=Diff,sigma_fun=sigma,KK_fun=KK, sig_fun = Sig) {
    W = dnorm(Diff_fun,0,sigma_fun[l]*sig_fun)*Sig*sigma_fun
    D_Kernels[[KK_fun+l]] = Matrix((W + t(W)) / 2, sparse=TRUE, doDiag=FALSE)
    return(D_Kernels)

  }))

  stopCluster(cl)

  # compute D_Kernels
  for (i in 1:length(D_Kernels)) {
    K = D_Kernels[[i]]
    k = 1/sqrt(diag(K)+1)
    G = K * (k %*% t(k))
    G1 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
    G2 = t(G1)
    D_Kernels_tmp = (G1 + G2 - 2*G)/2
    D_Kernels_tmp = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
    Kernels_tmp = matrix(1, nrow = nrow(D_Kernels_tmp), ncol = ncol(D_Kernels_tmp)) - D_Kernels_tmp
    D_Kernels[[i]] = Matrix(Kernels_tmp, sparse=TRUE, doDiag=FALSE)
  }

  return(D_Kernels)

}

dominateset <- function(xx,KK=20) {
  ### This function outputs the top KK neighbors.
  zero <- function(x) {
    s = sort(x, index.return=TRUE)
    x[s$ix[1:(length(x)-KK)]] = 0
    return(x)
  }
  normalize <- function(X) X / rowSums(X)
  A = matrix(0,nrow(xx),ncol(xx));
  for(i in 1:nrow(xx)){
    A[i,] = zero(xx[i,]);
  }

  return(A)
}

dominateset2 <- function(xx,KK=20) {
  ### This function outputs the top KK neighbors.
  zero <- function(x) {
    s = sort(x, index.return=TRUE)
    x[s$ix[1:(length(x)-KK)]] = 0
    return(x)
  }

  return(A)
}

"eig1" <- function( A, c = NA, isMax = NA, isSym = NA ) {
  ## Calculate eigen vector and eigen value given dimension
  # set the needed parameters
  if(is.na(c)) {
    c = dim(A)[1]
  }
  if(c>dim(A)[1]) {
    c = dim(A)[1]
  }
  if(is.na(isMax)) {
    isMax = 1
  }
  if(is.na(isSym)) {
    isSym = 1
  }

  # compute the eigenvalues and eigenvectors of A
  if(isSym==1) {
    eigen_A = eigen(A,symmetric=TRUE)
  }
  else {
    eigen_A = eigen(A)
  }
  v = eigen_A$vectors
  d = eigen_A$values

  # sort the eigenvectors
  if(isMax == 0) {
    eigen_A_sorted = sort(d,index.return=TRUE)
  }
  else {
    eigen_A_sorted = sort(d,decreasing=TRUE,index.return=TRUE)
  }
  d1 = eigen_A_sorted$x
  idx = eigen_A_sorted$ix
  idx1 = idx[1:c]

  # compute the results
  eigval = d[idx1]
  eigvec = Re(v[,idx1])
  eigval_full = d[idx]

  return(list(eigval=eigval,eigvec=eigvec,eigval_full=eigval_full))

}

"umkl2" = function( D, beta = NA,u = 2 ) {

  # set some parameters
  if(is.na(beta)) {
    beta = 1 / length(D)
  }
  tol = 1e-4
  u = u
  logU = log(u)

  # compute Hbeta
  res_hbeta = Hbeta(D, beta)
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

"Hbeta" = function( D, beta ) {

  D = (D - min(D)) / (max(D) - min(D) + .Machine$double.eps)
  P = exp(-D * beta)
  sumP = sum(P)
  H = log(sumP) + beta * sum(D * P) / sumP
  P = P / sumP

  return(list(H=H,P=P))

}

silhouette_similarity = function (group, similarity_matrix) {
  similarity_matrix = as.matrix(similarity_matrix)
  similarity_matrix <- (similarity_matrix + t(similarity_matrix))/2
  diag(similarity_matrix) = 0
  normalize <- function(X) X/rowSums(X)
  similarity_matrix <- normalize(similarity_matrix)
  n <- length(group)
  if (!all(group == round(group)))
    stop("'group' must only have integer codes")
  cluster_id <- sort(unique(group <- as.integer(group)))
  k <- length(cluster_id)
  if (k <= 1 || k >= n)
    return(NA)
  doRecode <- (any(cluster_id < 1) || any(cluster_id > k))
  if (doRecode)
    group <- as.integer(fgroup <- factor(group))
  cluster_id <- sort(unique(group))
  wds <- matrix(NA, n, 3, dimnames = list(names(group), c("cluster",
                                                          "neighbor", "sil_width")))
  for (j in 1:k) {
    index <- (group == cluster_id[j])
    Nj <- sum(index)
    wds[index, "cluster"] <- cluster_id[j]
    dindex <- rbind(apply(similarity_matrix[!index, index,
                                            drop = FALSE], 2, function(r) tapply(r, group[!index],
                                                                                 mean)))
    maxC <- apply(dindex, 2, which.max)
    wds[index, "neighbor"] <- cluster_id[-j][maxC]
    s.i <- if (Nj > 1) {
      a.i <- colSums(similarity_matrix[index, index])/(Nj -
                                                         1)
      b.i <- dindex[cbind(maxC, seq(along = maxC))]
      ifelse(a.i != b.i, (a.i - b.i)/pmax(b.i, a.i), 0)
    }
    else 0
    wds[index, "sil_width"] <- s.i
  }
  attr(wds, "Ordered") <- FALSE
  class(wds) <- "silhouette"
  wds
}


#########
# GL_SIMLR function: (global and local simlarity learning)
#
# input:
#   data_list: a list of the data types to be integrated, the data columns are samples, rows are features, all data types should have same number of columns
#   NITER: max interation times
#   c: number of clusters, estimated before
#   B: number of nearest neighbors,
#   rho: tunning parameter for weight entrophy
#   beta: the hyperparameter for local importance
#   gamma: hyper parameter for weight term
# output:
#   w_list: the contribution weight of every kernel
#   k: the final weighted and normalized kernel
#   cluster_result: spectral Clustering result on the weighted kernel matrix

#######
#Single kernel
## test version, output iteration steps, remove for final function
GL_SIMLR_1 = function(data_list,dk_list = NA, NITER = 50, c = 4, B = 20,rho = 0.1,beta = 1, gamma = 1, kernel_type = 1,is.normal = F){

  ## Calculate distance matrix for all data ttpes
  kernel_list = list()
  if(is.na(dk_list)){
    for(i in 1:length(data_list)){
      mk = kernel_calculating(data_list[[i]], is.normal = is.normal)
      if(kernel_type == 1){
        kernel_list[[i]] = (max(mk)-mk)
      } else if (kernel_type == 2){
        kernel_list[[i]] = (2-mk)
      }

      kernel_list[[i]] = (kernel_list[[i]]+t(kernel_list[[i]]))/2

    }
  }else{
    for(i in 1:length(dk_list)){
      if(kernel_type == 1){
        kernel_list[[i]] = (max(dk_list[[i]])-dk_list[[i]])
      } else if (kernel_type == 2){
        kernel_list[[i]] = (2-dk_list[[i]])
      }
    }
  }

  ### For local matrix
  local_list = NULL
  for(i in 1:length(kernel_list)){
    local_list[[i]] = dominateset(kernel_list[[i]],KK = B)
  }


  ## Normalization
  for(i in 1:length(kernel_list)){
    kernel_list[[i]] = normalize_kernel(kernel_list[[i]])
  }

  ## Calculate kernel terms
  global_F = rep(0,length(kernel_list))
  local_F = rep(0,length(kernel_list))
  for (i in 1:length(kernel_list)){
    for (j in 1:length(kernel_list)){
      global_F[i] = global_F[i] + sum(kernel_list[[i]]*(kernel_list[[j]])/sum(kernel_list[[j]]^2))
      local_F[i] = local_F[i] + sum(local_list[[i]]*(local_list[[j]])/sum(local_list[[j]]^2))
    }
  }

  ## initialization

  ### For kernels

  k0 = matrix(0, nrow(kernel_list[[1]]),nrow(kernel_list[[1]]))
  for (i in 1:length(kernel_list)){
    k0 = k0+kernel_list[[i]]
  }

  k0 = k0/length(kernel_list)
  w_list = rep(1/length(kernel_list),length(kernel_list))
  D0 = diag(apply(k0,MARGIN=2,FUN=sum))
  L0 = D0 - k0

  eig1_res = eig1(L0,c,0)
  H = eig1_res$eigvec
  temp_eig1 = eig1_res$eigval
  evs_eig1 = eig1_res$eigval_full
  k = k0

  ## start iteration
  converge = vector()
  w_list_all = rep(1/length(kernel_list),length(kernel_list))
  cluster_all = kmeans(as.matrix(H),c,nstart = 200)$cluster
  for(iter in 1:NITER){

    ## given H update w
    for(i in 1:length(w_list)){
      w_list[i] = exp((global_F[i] + beta*local_F[i] + gamma*sum(diag(t(H)%*%kernel_list[[i]]%*%H)))/rho)
      # w_list[i] = exp((global_F[i]+local_F[i]+sum(diag(t(H)%*%normalize_kernel(kernel_list[[i]])%*%H)))/rho)
    }
    w_list = w_list/sum(w_list)
    w_list_all = rbind(w_list_all, w_list)

    ##  update new k_old
    k_old = k
    k = matrix(0,nrow(kernel_list[[1]]),nrow(kernel_list[[1]]))
    for (i in 1:length(kernel_list)){
      k = k + w_list[[i]]*kernel_list[[i]]
    }

    ### set k to be symmetric and normalized
    # k = normalize_kernel(k)
    # k = (k+t(k))/2

    ## given w update H
    D = diag(apply(k, MARGIN = 2, FUN = sum))
    L = D - k
    H_old = H

    eig1_res = eig1(L,c,0)
    H = eig1_res$eigvec
    temp_eig1 = eig1_res$eigval
    evs_eig1 = eig1_res$eigval_full
    cluster_all = rbind(cluster_all, kmeans(as.matrix(H),c,nstart = 200)$cluster)
    ## iteration judge
    converge[iter] = sum(evs_eig1[1:(c+1)]) - sum(evs_eig1[1:c])
    # converge[iter] = evs_eig1[c] - evs_eig1[c-1]
    if(iter>5){
      if(converge[iter] > converge[iter-1]){
        k = k_old
        if(converge[iter]>0.2){
          warning('Maybe you should set a larger value of c.')
        }
        break
      }
    }
  }
  cluster = kmeans(as.matrix(H_old),c,nstart = 200)$cluster
  return(list(cluster = cluster, k = k, w_list = w_list,converge = converge,w_list_all = w_list_all, cluster_all = cluster_all))
}
# data_list
# D_list = NA
# NITER = 50
# c = 4
# B = 20
# rho = 0.1
# gamma = NA
# beta = NA
# tol = 1e-10
# normal_type = 1
# F_norm = F
GL_SIMLR_4 = function(data_list, D_list = NA, NITER = 50, c = 4, B = 20,rho = 0.1,gamma = NA,beta = NA,tol = 1e-10,normal_type = 1, F_norm = F,standardize_term = T){

  ## Calculate kernel matrix and distance matrix for all data types
  kernel_list = list() ## the kerenl matrix, using the density of normal distribution
  if(is.na(D_list)){
    D_list = list() ## using kernel to calculate distance here: Dij = kii+kjj-2kij
    for(i in 1:length(data_list)){
      kernel_list[[i]] = kernel_calculating(data_list[[i]], is.SP = F)
      D_list[[i]] = kernel_calculating(data_list[[i]], is.SP = T)
      # D_list[[i]] = dist2(data_list[[i]])
      kernel_list[[i]] = (kernel_list[[i]]+t(kernel_list[[i]]))/2
    }
  }else{
    for(i in 1:length(D_list)){
      kernel_list[[i]] = dist_kernel(D_list[[i]], is.dist = F)
      D_list[[i]] = dist_kernel(D_list[[i]], is.dist = T)
    }
  }

  ## Normalization
  for(i in 1:length(kernel_list)){
    if(normal_type == 1){
      kernel_list[[i]] = normalize_kernel(kernel_list[[i]])
    }else{
      D = diag(rowSums(kernel_list[[i]]))
      kernel_list[[i]] = solve(D^(1/2)) %*% kernel_list[[i]] %*% solve(D^(1/2))
    }

  }

  ### For local matrix
  local_list = NULL
  for(i in 1:length(kernel_list)){
    local_list[[i]] = dominateset(kernel_list[[i]],KK = B)
  }


  ## Calculate kernel terms
  global_F = rep(0,length(kernel_list))
  local_F = rep(0,length(kernel_list))
  for (i in 1:length(kernel_list)){
    for (j in 1:length(kernel_list)){
      if(F_norm){
        global_F[i] = global_F[i] + sum((kernel_list[[i]]/sum(kernel_list[[i]]^2))*(kernel_list[[j]]/sum(kernel_list[[j]]^2)))
        local_F[i] = local_F[i] + sum((local_list[[i]]/sum(local_list[[i]]^2))*(local_list[[j]]/sum(local_list[[j]]^2)))
      }else{
        global_F[i] = global_F[i] + sum(kernel_list[[i]]*kernel_list[[j]])
        local_F[i] = local_F[i] + sum(local_list[[i]]*local_list[[j]])
      }
    }
  }
  global_local_term = tibble(global = global_F , local = local_F)



  ## initialization

  ### For kernels
  k0 = matrix(0, nrow(kernel_list[[1]]),nrow(kernel_list[[1]]))
  for (i in 1:length(kernel_list)){
    k0 = k0+kernel_list[[i]]
  }

  k0 = k0/length(kernel_list)
  w_list = rep(1/length(kernel_list),length(kernel_list))
  d0 = apply(k0,MARGIN=1,FUN=sum)
  d0[d0 == 0] = .Machine$double.eps
  D0 = diag(d0)
  L0 = D0-k0
  L0 =  diag(1/sqrt(d0))%*% L0 %*%  diag(1/sqrt(d0))
  # L0 = matrix(1, dim(D0)[1],dim(D0)[2]) - solve(D0^(1/2))%*%k0%*%solve(D0^(1/2))

  eig1_res = eig1(L0,c,0)
  H = eig1_res$eigvec
  H_c = t(apply(as.matrix(H), MARGIN = 1, FUN = function(x){x/sqrt(sum(x^2))}))
  temp_eig1 = eig1_res$eigval
  evs_eig1 = eig1_res$eigval_full
  k = k0

  ##  calculate hyper parameter initialization
  ### Calculate weighted distance
  estimate_gamma = ifelse(is.na(gamma),T,F)
  estimate_beta = ifelse(is.na(beta),T,F)
  D_w = matrix(0, nrow(D_list[[1]]),ncol(D_list[[1]]))
  for (i in 1:length(D_list)){
    D_w = D_w+D_list[[i]]*w_list[[i]]
  }

  D_w_sort = t(apply(D_w, MARGIN = 1, FUN = sort))
  par_t = sum(apply(D_w_sort, MARGIN = 1, FUN = function(x){sum(x[(B+2)]-x[2:(B+1)])}))/(2*nrow(D_list[[1]]))
  gamma = ifelse(estimate_gamma, par_t, gamma)
  beta = ifelse(estimate_beta, par_t, beta)



  ## start iteration
  converge = vector()
  w_list_all = rep(1/length(kernel_list),length(kernel_list))
  cluster_all = kmeans(H_c,c,nstart = 200)$cluster
  GL_trace_all = NULL
  H_sum = NULL
  L = L0
  delta_eig = 0
  for(iter in 1:NITER){
    ## given H update w
    GL_trace3 = NULL
    for(i in 1:length(w_list)){
      L_i = diag(rowSums(kernel_list[[i]]))-kernel_list[[i]]
      GL_trace = sum(diag(t(H)%*%L_i%*%H))
      w_list[i] = exp((beta*(global_F[i] + local_F[i]) - gamma*GL_trace)/rho)
      GL_trace3 = c(GL_trace3, GL_trace)
    }
    if(standardize_term){
      w_list = exp((beta*(global_F-mean(global_F) + local_F-mean(local_F)) - gamma*(GL_trace3-mean(GL_trace3)))/rho)
    }
    H_sum = c(H_sum, sum(H))
    w_list = w_list/sum(w_list)
    w_list_all = rbind(w_list_all, w_list)
    GL_trace_all = rbind(GL_trace_all,GL_trace3)

    ##  update new k_old
    k_old = k
    k = matrix(0,nrow(kernel_list[[1]]),nrow(kernel_list[[1]]))
    s = matrix(0,nrow(kernel_list[[1]]),nrow(kernel_list[[1]]))
    for (i in 1:length(kernel_list)){
      k = k + w_list[[i]]*kernel_list[[i]]
      s = s + w_list[[i]]*local_list[[i]]
    }

    ### set k to be symmetric and normalized
    # k = normalize_kernel(k)
    # k = (k+t(k))/2

    ## given w update H
    d = apply(k, MARGIN = 1, FUN = sum)
    d[d==0] = .Machine$double.eps
    D = diag(d)
    L = D-k
    L =  diag(1/sqrt(d))%*% L %*%  diag(1/sqrt(d))
    H_old = H

    eig1_res = eig1(L,c,0)
    H = eig1_res$eigvec
    H_c = t(apply(as.matrix(H), MARGIN = 1, FUN = function(x){x/sqrt(sum(x^2))}))
    temp_eig1 = eig1_res$eigval
    evs_eig1 = eig1_res$eigval_full
    cluster_all = rbind(cluster_all, kmeans(H_c,c,nstart = 200)$cluster)



    ## iteration judge
    obj_1 = 0
    for(i in 1:length(w_list)){
      obj_1 = obj_1 +
        beta*(sum((as.matrix(k)/sum(as.matrix(k)^2))*(kernel_list[[i]]/sum(kernel_list[[i]]^2)))+
                sum((as.matrix(s)/sum(as.matrix(s)^2))*(local_list[[i]]/sum(local_list[[i]]^2))))
    }
    obj_2 = gamma*sum(diag(t(H)%*%L%*%H))
    obj_3 = rho*sum(w_list*log(w_list))
    converge[iter] = obj_1+obj_2+obj_3

    if(iter<=10){
      delta_eig = ifelse(evs_eig1[length(evs_eig1)] > 1e-6, 1,0)
    }else{
      if(converge[iter] - converge[iter-1] < tol){
        k = k_old
        # if(converge[iter]>0.2){
        #   warning('Maybe you should set a larger value of c.')
        # }
        break
      }
    }
    ## update hyperparameter
    D_w = matrix(0, nrow(D_list[[1]]),ncol(D_list[[1]]))
    for (i in 1:length(D_list)){
      D_w = D_w+D_list[[i]]*w_list[[i]]
    }

    D_w_sort = t(apply(D_w, MARGIN = 1, FUN = sort))
    par_t = sum(apply(D_w_sort, MARGIN = 1, FUN = function(x){sum(x[(B+2)]-x[2:(B+1)])}))/(2*nrow(D_list[[1]]))
    gamma = ifelse(estimate_gamma, par_t*(1+0.5*delta_eig), gamma)
    beta = ifelse(estimate_beta, par_t, beta)
  }


  ## clustering
  cluster = kmeans(as.matrix(H),c,nstart = 200)$cluster
  return(list(cluster = cluster,
              k = k,
              w_list = tibble(w_list),
              GL_trace = tibble(GL_trace3),
              GL_trace_all = GL_trace_all,
              global_local_term = global_local_term,
              converge = converge,
              w_list_all = w_list_all,
              cluster_all = cluster_all,
              H_sum = H_sum,
              beta = beta,
              gamma = gamma))
}



GL_SIMLR_tunning =  function(data_list = NA, D_list = NA, c, gamma_grid = NA, beta_grid = NA,rho_grid = 0.1,true_label = NA){
  ## In this setting, gamma and beta must have the same length, and set grid in parameter tunning pool
  ## parameter grid
  m = 1
  par_list = NULL
  if(!is.na(gamma_grid)){
    for(i in 1:length(beta_grid)){
      for(j in 1:length(gamma_grid)){
        for(k in 1:length(rho_grid)){
          par_list[[m]] = list(gamma = gamma_grid[j],
                               beta = beta_grid[i],
                               rho = rho_grid[k])
          m = m+1
        }
      }
    }
  }else{
    for(k in 1:length(rho_grid)){
      par_list[[m]] = list(rho = rho_grid[k])
      m = m+1
    }
  }

  # start parameter selection
  ## Using silouhette as criteria (need further justification of the R package)
  tunning_crit = NULL
  for(j in 1:length(par_list)){
    if(j %% 10 == 0){ print(j)}
    if(!is.na(gamma_grid)) {
      gamma = par_list[[j]]$gamma
      beta = par_list[[j]]$beta
    }else{
      gamma = NA
      beta = NA
    }
    rho = par_list[[j]]$rho
    tryCatch({
      gl_res = GL_SIMLR_4(data_list = data_list, D_list = D_list, c = c, gamma = gamma,beta = beta, rho = rho)
      sil_ind = summary(silhouette_similarity(gl_res$cluster, gl_res$k))$avg.width
      w1 = gl_res$w_list[1,1][[1]]
      w2 = gl_res$w_list[2,1][[1]]
      if(!is.na(true_label)){
        nmi_ind = igraph::compare(gl_res$cluster, true_label, method = "nmi")
        rand_ind = igraph::compare(gl_res$cluster, true_label, method = "rand")
        crit = c(gamma = gl_res$gamma,beta = gl_res$beta, rho = rho, silourtte = sil_ind,nmi = nmi_ind, rand = rand_ind, w1 = w1,w2 = w2)
      }else{
        crit = c(gamma = gl_res$gamma, beta = gl_res$beta,rho = rho, silourtte = sil_ind, w1 = w1,w2 = w2)
      }


      ## collection
      tunning_crit = rbind(tunning_crit, crit)
    },
    error = function(e){cat("!!!!!!!!!!!!!!!!!!!",conditionMessage(e), rho)}
    )
  }

  tunning_crit = as_tibble(tunning_crit)
  rho_opt = tunning_crit %>% filter(silourtte == max(silourtte)) %>% pull(rho)
  if(!is.na(gamma_grid)){
    gamma_opt = tunning_crit %>% filter(silourtte == max(silourtte)) %>% pull(gamma)
    beta_opt = tunning_crit %>% filter(silourtte == max(silourtte)) %>% pull(beta)
  }else(
    gamma_opt = NA
  )
  res_opt = GL_SIMLR_4(data_list = data_list, D_list = D_list, c = c, gamma = gamma_opt[1],beta = beta_opt[1], rho = rho_opt)

  ## best performance
  return(list(tunning_crit = tunning_crit,
              res_opt = res_opt,
              rho_opt = rho_opt,
              gamma_opt = res_opt$gamma,
              beta_opt = beta_opt,
              gl_opt = sum(res_opt$global_local_term),
              GL_trace_opt = sum(res_opt$GL_trace)))
}


#######
# multiple kernel
GL_SIMLR_m = function(data_list, NITER = 50, c = 4, B = 20,rho = 0.1,gamma = 1,beta = 1,tol = 1e-10,kernel_type = 1,is.normal = F){


  ## Calculate normalize kernels for all data types
  kernel_list = list()
  for(i in 1:length(data_list)){
    mk = multiple.kernel2(data_list[[i]])
    kernel_list = c(kernel_list, mk)
    # kernel_list[[i]] =lapply(mk, FUN = normalize_kernel)
  }

  ### For local matrix
  local_list = NULL
  for(i in 1:length(kernel_list)){
    local_list[[i]] = dominateset(kernel_list[[i]],KK = B)
  }


  ## Normalization
  for(i in 1:length(kernel_list)){
    kernel_list[[i]] = normalize_kernel(kernel_list[[i]])
  }

  ## Calculate kernel terms
  global_F = rep(0,length(kernel_list))
  local_F = rep(0,length(kernel_list))
  for (i in 1:length(kernel_list)){
    for (j in 1:length(kernel_list)){
      global_F[i] = global_F[i] + sum(kernel_list[[i]]*(kernel_list[[j]])/sum(kernel_list[[j]]^2))
      local_F[i] = local_F[i] + sum(local_list[[i]]*(local_list[[j]])/sum(local_list[[j]]^2))
    }
  }

  ## initialization

  ### For kernels
  k0 = matrix(0, nrow(kernel_list[[1]]),nrow(kernel_list[[1]]))
  for (i in 1:length(kernel_list)){
    k0 = k0+kernel_list[[i]]
  }

  k0 = k0/length(kernel_list)
  w_list = rep(1/length(kernel_list),length(kernel_list))
  D0 = diag(apply(k0,MARGIN=2,FUN=sum))
  L0 = D0 - k0

  eig1_res = eig1(L0,c,0)
  H = eig1_res$eigvec
  temp_eig1 = eig1_res$eigval
  evs_eig1 = eig1_res$eigval_full
  k = k0


  ## start iteration
  converge = vector()
  w_list_all = rep(1/length(kernel_list),length(kernel_list))
  cluster_all = kmeans(as.matrix(H),c,nstart = 200)$cluster
  for(iter in 1:NITER){

    ## given H update w
    for(i in 1:length(w_list)){
      w_list[i] = exp((global_F[i] + beta*local_F[i] + gamma*sum(diag(t(H)%*%kernel_list[[i]]%*%H)))/rho)
      # w_list[i] = exp((global_F[i]+local_F[i]+sum(diag(t(H)%*%normalize_kernel(kernel_list[[i]])%*%H)))/rho)
    }
    w_list = w_list/sum(w_list)
    w_list_all = rbind(w_list_all, w_list)

    ##  update new k_old
    k_old = k
    k = matrix(0,nrow(kernel_list[[1]]),nrow(kernel_list[[1]]))
    s = matrix(0,nrow(kernel_list[[1]]),nrow(kernel_list[[1]]))
    for (i in 1:length(kernel_list)){
      k = k + w_list[[i]]*kernel_list[[i]]
      s = s + w_list[[i]]*local_list[[i]]
    }

    ## set k to be symmetric and normalized
    k = normalize_kernel(k)
    k = (k+t(k))/2

    ## given w update H
    D = diag(apply(k, MARGIN = 2, FUN = sum))
    L = D - k
    H_old = H

    eig1_res = eig1(L,c,0)
    H = eig1_res$eigvec
    temp_eig1 = eig1_res$eigval
    evs_eig1 = eig1_res$eigval_full
    cluster_all = rbind(cluster_all, kmeans(as.matrix(H),c,nstart = 200)$cluster)

    ## iteration judge
    obj_1 = 0
    for(i in 1:length(w_list)){
      obj_1 = obj_1 + sum(k*(kernel_list[[i]]/sum(kernel_list[[i]]^2)))+beta* sum(s*(local_list[[i]]/sum(local_list[[i]]^2)))
    }
    obj_2 = sum(diag(t(H)%*%k%*%H))
    obj_3 = sum(w_list*log(w_list))
    converge[iter] = obj_1+obj_2+obj_3
    # converge[iter] = evs_eig1[c] - evs_eig1[c-1]
    if(iter>5){
      if(converge[iter] - converge[iter-1] < tol){
        k = k_old
        # if(converge[iter]>0.2){
        #   warning('Maybe you should set a larger value of c.')
        # }
        break
      }
    }
  }
  cluster = kmeans(as.matrix(H_old),c,nstart = 200)$cluster
  return(list(cluster = cluster, k = k, w_list = w_list,converge = converge,w_list_all = w_list_all, cluster_all = cluster_all))
}

