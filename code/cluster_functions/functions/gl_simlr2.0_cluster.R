# function: gl_simlr
# goal: using kernel alignment and trace term to find weight between kernels to construct similarity matrix
library(parallel)
library(clValid)
library(dplyr)
library(SNFtool)
library(igraph)
library(Matrix)
library(quadprog)
gl_simlr = function(data_list = NA, dist_list = NA, kernel_list = NA, B = 30, c = 4, max = 30, tol = 1e-10,
                    rho = 0.1,gamma = 0.5,beta = 0.5,delta = 0,kernel_type = "Dk",sigma_dk = 2,
                    mu_aff = 0.5, isdk = F, standardization = 2){

  # kernel calculation
  if(is.na(kernel_list)){
    if(is.na(dist_list)){
      if(kernel_type == "Dk"){
        dist_kernel_list = lapply(data_list, FUN = function(x) D_kernels(x_fun = x,allk_fun = B,sigma_fun = sigma_dk,standardization = standardization) )
        if(isdk){
          kernel_list = lapply(dist_kernel_list, FUN = function(x){x[[2]]})
        }else{
          kernel_list = lapply(dist_kernel_list, FUN = function(x){x[[1]]})
        }

      }else if(kernel_type == "affinity"){
        if(standardization == 1){
          data_list = lapply(data_list, standardNormalization)
        }else if(standardization == 2){
          data_list = lapply(data_list, function(x) {
            apply(x, 2, function(y) {
              (y-min(y))/(max(y)-min(y))
            })
          })
        }
        kernel_list = lapply(data_list, FUN = function(x) affinityMatrix(dist2(x)^(1/2),sigma = mu_aff))
      }
    }else{ ## didn't test whether applicable to dist kernel
      if(kernel_type == "Dk"){
        dist_kernel_list = lapply(data_list, FUN = function(x) D_kernels(Diff_fun = x,allk_fun = B,sigma_fun = sigma_dk,standardization = standardization))
        if(isdk){
          kernel_list = lapply(dist_kernel_list, FUN = function(x){x[[2]]})
        }else{
          kernel_list = lapply(dist_kernel_list, FUN = function(x){x[[1]]})
        }
      }else if(kernel_type == "affinity"){
        kernel_list = lapply(dist_list, FUN = function(x) affinityMatrix(as.matrix(x)^(1/2)))
      }
    }
  }

  # kernel_list = lapply(dist_kernel_list, FUN = function(x){x[[1]]})
  local_kernel_list = lapply(kernel_list, FUN = function(x) dominateset(x, KK = B))
  # D_list = lapply(dist_kernel_list, FUN = function(x){x[[2]]})

  # Kernel normalization
  norm_kernel_list = lapply(kernel_list, FUN = function(x){
    D = diag(rowSums(x))
    solve(D^(1/2))%*%as.matrix(x)%*%solve(D^(1/2))
  }
  )

  norm_local_kernel_list = lapply(local_kernel_list, FUN = function(x){
    D = diag(rowSums(x))
    solve(D^(1/2))%*%as.matrix(x)%*%solve(D^(1/2))
  }
  )
  # Global local term calculation
  global_term = kernel_similarity(norm_kernel_list)
  local_term = kernel_similarity(norm_local_kernel_list)

  # w0, H0 initialization
  M = length(kernel_list)
  n = dim(kernel_list[[1]])[1]
  w = rep(1/M,M)
  K = matrix(0,n,n)


  for(i in 1:M){
    K = K+norm_kernel_list[[i]]*w[i]
  }

  D = diag(rowSums(K))
  L = D-K
  L = solve(D^(1/2))%*%L%*%solve(D^(1/2))
  eig_res = eig1(L, c, isMax = 0)
  H = eig_res$eigvec

  # iteration
  criteria = NULL
  w_all = NULL
  trace_term = NULL
  for(t in 1:max){
    ## iteration of new w given H
    trace = NULL
    for(i in 1:M){
      L_i = diag(rowSums(norm_kernel_list[[i]])) - norm_kernel_list[[i]]
      trace = c(trace, sum(diag(t(H)%*%L_i%*%H)))
      w[i]=exp((delta*(global_term[i]+beta*local_term[i])-gamma*sum(diag(t(H)%*%L_i%*%H)))/rho)
    }
    w = w/sum(w)
    w_all = rbind(w_all, w)
    trace_term = rbind(trace_term, trace)

    ## iteration of H given w
    K = matrix(0,n,n)
    for(i in 1:M){
      K = K+norm_kernel_list[[i]]*w[i]
    }
    D = diag(rowSums(K))
    L = D-K
    L = solve(D^(1/2))%*%L%*%solve(D^(1/2))
    eig_res = eig1(L, c, isMax = 0)
    H = eig_res$eigvec

    ## stop  critieria
    obj = -delta*(sum(w*global_term+w*beta*local_term))+gamma*sum(diag(t(H)%*%L%*%H))+rho*sum(w*log(w))
    criteria = c(criteria,obj)
    if(t > 10){
      if ((criteria[t]-criteria[t-1])<tol){
        cluster = kmeans(H,c,nstart = 200)$cluster
        break
      }
    }

  }
  return(list(w = w,
              K = K,
              w_all = w_all,
              trace_term = trace_term,
              global_term = global_term,
              local_term = local_term,
              cluster = cluster,
              iteration = t))
}


Crazy_5term_clust = function(data_list = NA, # a list of data matrices to integrate
                             kernel_list = NA,
                             k = 30, # #neighbors
                             c = 3, # #clusters
                             kernel_fun = "D_kernels", # the definition of kernels
                             rho = 0.1, # tuning parameter for entropy penalty
                             alpha0,# tuning parameter for <S, K>
                             alpha, # tuning parameter for <S, W>
                             beta = NA, # tuning parameter for ||S||
                             gamma = NA, # tuning parameter for laplacian
                             stopping = 10^-3, # stopping rule
                             n_ite = 50, # max number of iterations
                             normalize = 2,
                             print_details = F, # whether to print stepwise details
                             self_weight = F,
                             umkl = F
)
{
  if(is.na(kernel_list)){
    if(kernel_fun == "D_kernels"){
      Kmat_ls = lapply(data_list, function(x) as.matrix(D_kernels(x_fun = x)$ori_Kernels))
    }else if(kernel_fun == "affinity"){
      Kmat_ls = lapply(data_list, function(x) affinityMatrix(dist2(x)^(1/2)))
    }
  }else{
    Kmat_ls = kernel_list
  }


  # normalization
  s_Kmat_ls = Kmat_ls # scaled kernel matrices
  s_Qmat_ls = lapply(s_Kmat_ls, FUN = function(x) dominateset(x, KK = k))
  if(normalize == 1){
    for (ll in 1:length(s_Kmat_ls)) {
      s_Kmat_ls[[ll]] = s_Kmat_ls[[ll]] / norm(s_Kmat_ls[[ll]], type = 'F')
      s_Qmat_ls[[ll]] = s_Qmat_ls[[ll]] / norm(s_Qmat_ls[[ll]], type = 'F')
    }
  }else if(normalize == 2){
    s_Kmat_ls = lapply(s_Kmat_ls, normalize_kernel)
    s_Qmat_ls = lapply(s_Qmat_ls, normalize_kernel)
  }

  if(umkl){
    Dk_ls = lapply(s_Kmat_ls, D_kernels_calc)
    Dq_ls = lapply(s_Qmat_ls, D_kernels_calc)
  }

  n_sam = nrow(s_Kmat_ls[[1]]) # number of samples
  n_ker = length(s_Kmat_ls) # number of kernels
  I_n = diag(n_sam) # diagonal matrix
  old_w = rep(1/n_ker, n_ker) # initialize weights
  old_S = matrix(0, nrow = n_sam, ncol = n_sam) # initialize S
  for (ll in 1:n_ker) {old_S = old_S + old_w[ll] * s_Kmat_ls[[ll]]}
  # old_L = matrix(0, nrow = n_sam, ncol = c) # initialize L
  old_eig = eigen(old_S)
  old_L = old_eig$vectors[,order(old_eig$values)][,1:c]
  # old_L[which(response == 1), 1] = old_L[which(response == 0), 2] = 1
  # old_L = t(t(old_L)/sqrt(colSums(old_L))) # scale old_L
  # A = case_control_mat(response) # case control status matrix
  # s_A = A / norm(A, type = 'F') # scaled A


  ### Start iteration
  for (iter in 1:n_ite) {
    w_s_Kmat = w_s_Qmat = matrix(0, n_sam, n_sam) # weighted average of scaled kernel matrices
    for (ll in 1:n_ker) {
      w_s_Kmat = w_s_Kmat + old_w[ll] * s_Kmat_ls[[ll]]
      w_s_Qmat = w_s_Qmat + old_w[ll] * s_Qmat_ls[[ll]]
    }

    # parameter estimation
    if(is.na(beta)&is.na(gamma)){
      D_w_Kmat = dist_kernels(w_s_Kmat)
      D_w_Qmat = dist_kernels(w_s_Qmat)
      par_est = (1/(2*n_sam))*sum(apply(alpha0*D_w_Kmat+alpha*D_w_Qmat, MARGIN = 1, function(x){
        s_x = sort(x)
        rep(s_x[k+2],k) - s_x[2:(k+1)]
      }))
      print(par_est)
      gamma = beta = par_est
    }
    ### Initialization
    new_w = old_w
    new_L = old_L
    new_S = old_S
    linear_terms_S = alpha0 * w_s_Kmat + alpha * w_s_Qmat + gamma * (old_L %*% t(old_L))

    ### Update S
    for (i in 1:n_sam) {
      QP_results = solve.QP(Dmat = I_n * beta * 2, # quadratic programming
                            dvec = linear_terms_S[i,],
                            Amat = t(rbind(rep(1,n_sam), I_n)),
                            bvec = c(1, rep(0,n_sam)),
                            meq = 1)
      new_S[,i] = QP_results$solution
    }
    new_S = (new_S+t(new_S))/2 # make sure S is symmetric

    ### Update L
    Laplacian = I_n - new_S
    eg_results = eigen(Laplacian) # eigen-decompositions
    new_L = eg_results$vectors[,order(eg_results$values)][,1:c]# extract the two eigenvectors

    ### Update w
    if(umkl){
        DD = vector()
        for (i in 1:length(Dk_ls)) {
          temp = (.Machine$double.eps + alpha0*Dk_ls[[i]]+alpha*Dq_ls[[i]]) *
            (new_S + .Machine$double.eps)
          DD[i] = mean(apply(temp, MARGIN = 2, FUN = sum))
        }
        new_w = umkl.cimlr.gl(DD)
        new_w = new_w/sum(new_w)
    }else{
      first_term = vector() # 1st term in optimization
      for (ll in 1:n_ker) {
        first_term[ll] = sum((alpha0*s_Kmat_ls[[ll]]+alpha*s_Qmat_ls[[ll]]) * new_S)
        if(self_weight){
          new_w[ll] = 1/(2*sqrt(first_term[ll])) # self_weight estimation
        }else{
          new_w[ll] = exp(first_term[ll] / rho) # the new weights for kernels
        }

      }
      new_w = new_w/sum(new_w)
    }
     # scale the kernels

    ### Print details
    if (print_details) {
      cat(paste0('Iteration ',iter, ':\n  Optimal weights: '))
      cat(new_w)
      cat('\n')
      opt_value = - sum(first_term) +
        rho * sum(new_w * log(new_w)) +
        #-alpha * sum(new_S * s_A) +
        beta * norm(new_S, type = 'F')^2 +
        gamma * sum(diag(t(new_L) %*% (I_n - new_S) %*% new_L))
      cat(paste0('  Criterion: ', opt_value, '\n'))
    }

    ### Whether to stop
    if (iter>=3 & max(abs(new_w-old_w))<= stopping) {break} else {
      old_w = new_w
      old_S = new_S
      old_L = new_L
    }
  }

  names(new_w) = names(s_Kmat_ls)
  return(list(w = new_w,
              S = new_S,
              L = new_L))
}

Crazy_5term_clust_dist = function(data_list = NA, # a list of data matrices to integrate
                                   kernel_list = NA,
                                   k = 30, # #neighbors
                                   c = 3, # #clusters
                                   kernel_fun = "D_kernels", # the definition of kernels
                                   rho = 0.1, # tuning parameter for entropy penalty
                                   alpha0,# tuning parameter for <S, K>
                                   alpha, # tuning parameter for <S, W>
                                   beta = NA, # tuning parameter for ||S||
                                   gamma = NA, # tuning parameter for laplacian
                                   stopping = 10^-3, # stopping rule
                                   n_ite = 50, # max number of iterations
                                   normalize = 2,
                                   print_details = F, # whether to print stepwise details
                                   self_weight = F,
                                   umkl = F
)
{
  if(is.na(kernel_list)){
    if(kernel_fun == "D_kernels"){
      Kmat_ls = lapply(data_list, function(x) as.matrix(D_kernels(x_fun = x)$ori_Kernels))
    }else if(kernel_fun == "affinity"){
      Kmat_ls = lapply(data_list, function(x) affinityMatrix(dist2(x)^(1/2)))
    }
  }else{
    Kmat_ls = kernel_list
  }


  # normalization
  s_Kmat_ls = Kmat_ls # scaled kernel matrices
  s_Qmat_ls = lapply(s_Kmat_ls, FUN = function(x) dominateset(x, KK = k))
  if(normalize == 1){
    for (ll in 1:length(s_Kmat_ls)) {
      s_Kmat_ls[[ll]] = s_Kmat_ls[[ll]] / norm(s_Kmat_ls[[ll]], type = 'F')
      s_Qmat_ls[[ll]] = s_Qmat_ls[[ll]] / norm(s_Qmat_ls[[ll]], type = 'F')
    }
  }else if(normalize == 2){
    s_Kmat_ls = lapply(s_Kmat_ls, normalize_kernel)
    s_Qmat_ls = lapply(s_Qmat_ls, normalize_kernel)
  }


  Dk_ls = lapply(s_Kmat_ls, D_kernels_calc)
  Dq_ls = lapply(s_Qmat_ls, D_kernels_calc)


  n_sam = nrow(s_Kmat_ls[[1]]) # number of samples
  n_ker = length(s_Kmat_ls) # number of kernels
  I_n = diag(n_sam) # diagonal matrix
  old_w = rep(1/n_ker, n_ker) # initialize weights
  old_S = matrix(0, nrow = n_sam, ncol = n_sam) # initialize S
  for (ll in 1:n_ker) {old_S = old_S + old_w[ll] * s_Kmat_ls[[ll]]}
  # old_L = matrix(0, nrow = n_sam, ncol = c) # initialize L
  old_eig = eigen(old_S)
  old_L = old_eig$vectors[,order(old_eig$values)][,1:c]
  # old_L[which(response == 1), 1] = old_L[which(response == 0), 2] = 1
  # old_L = t(t(old_L)/sqrt(colSums(old_L))) # scale old_L
  # A = case_control_mat(response) # case control status matrix
  # s_A = A / norm(A, type = 'F') # scaled A


  ### Start iteration
  for (iter in 1:n_ite) {
    w_s_Kmat = w_s_Qmat = matrix(0, n_sam, n_sam) # weighted average of scaled kernel matrices
    for (ll in 1:n_ker) {
      w_s_Kmat = w_s_Kmat + old_w[ll] * s_Kmat_ls[[ll]]
      w_s_Qmat = w_s_Qmat + old_w[ll] * s_Qmat_ls[[ll]]
    }

    # parameter estimation
    if(is.na(beta)&is.na(gamma)){
      D_w_Kmat = dist_kernels(w_s_Kmat)
      D_w_Qmat = dist_kernels(w_s_Qmat)
      par_est = (1/(2*n_sam))*sum(apply(alpha0*D_w_Kmat+alpha*D_w_Qmat, MARGIN = 1, function(x){
        s_x = sort(x)
        rep(s_x[k+2],k) - s_x[2:(k+1)]
      }))
      print(par_est)
      gamma = beta = par_est
    }
    ### Initialization
    new_w = old_w
    new_L = old_L
    new_S = old_S
    linear_terms_S = alpha0 * w_s_Kmat + alpha * w_s_Qmat + gamma * (old_L %*% t(old_L))

    ### Update S
    for (i in 1:n_sam) {
      QP_results = solve.QP(Dmat = I_n * beta * 2, # quadratic programming
                            dvec = linear_terms_S[i,],
                            Amat = t(rbind(rep(1,n_sam), I_n)),
                            bvec = c(1, rep(0,n_sam)),
                            meq = 1)
      new_S[,i] = QP_results$solution
    }
    new_S = (new_S+t(new_S))/2 # make sure S is symmetric

    ### Update L
    Laplacian = I_n - new_S
    eg_results = eigen(Laplacian) # eigen-decompositions
    new_L = eg_results$vectors[,order(eg_results$values)][,1:c]# extract the two eigenvectors

    ### Update w
    if(umkl){
      DD = vector()
      for (i in 1:length(Dk_ls)) {
        temp = (.Machine$double.eps + alpha0*Dk_ls[[i]]+alpha*Dq_ls[[i]]) *
          (new_S + .Machine$double.eps)
        DD[i] = mean(apply(temp, MARGIN = 2, FUN = sum))
      }
      new_w = umkl.cimlr.gl(DD)
      new_w = new_w/sum(new_w)
    }else{
      first_term = vector() # 1st term in optimization
      for (ll in 1:n_ker) {
        first_term[ll] = sum((alpha0*s_Kmat_ls[[ll]]+alpha*s_Qmat_ls[[ll]]) * new_S)
        if(self_weight){
          new_w[ll] = 1/(2*sqrt(first_term[ll])) # self_weight estimation
        }else{
          new_w[ll] = exp(first_term[ll] / rho) # the new weights for kernels
        }

      }
      new_w = new_w/sum(new_w)
    }
    # scale the kernels

    ### Print details
    if (print_details) {
      cat(paste0('Iteration ',iter, ':\n  Optimal weights: '))
      cat(new_w)
      cat('\n')
      opt_value = - sum(first_term) +
        rho * sum(new_w * log(new_w)) +
        #-alpha * sum(new_S * s_A) +
        beta * norm(new_S, type = 'F')^2 +
        gamma * sum(diag(t(new_L) %*% (I_n - new_S) %*% new_L))
      cat(paste0('  Criterion: ', opt_value, '\n'))
    }

    ### Whether to stop
    if (iter>=3 & max(abs(new_w-old_w))<= stopping) {break} else {
      old_w = new_w
      old_S = new_S
      old_L = new_L
    }
  }

  names(new_w) = names(s_Kmat_ls)
  return(list(w = new_w,
              S = new_S,
              L = new_L))
}

Crazy_5term_clust_tune = function(kernel_list,true_label, c = 4, beta_list = c(0.01,0.1), gamma_list = c(0.01,0.1)){
  par = list()
  m = 1
  for(i in 1:length(beta_list)){
    for(j in 1:length(gamma_list)){
      par[[m]] = list(beta = beta_list[i],
                      gamma = gamma_list[j])
      m = m+1
    }
  }
  # print(par)
  res_list = map(par, function(x){tryCatch({Crazy_5term_clust(kernel_list = kernel_list,
                                                              c = c,
                                                              alpha = 1,
                                                              alpha0 = 3,
                                                              rho = 0.1,
                                                              normalize = 0,
                                                              beta = x$beta,
                                                              gamma = x$gamma
  )},
  error = function(e){cat("!!!!!!!!!!!!!!!!!!!",conditionMessage(e))})})
  compare = map_dfr(map(res_list, "S"),function(x) {
    cluster = spectralClustering(as.matrix(x),c)
    sil = summary(silhouette_similarity(cluster,as.matrix(x)))$avg.width
    nmi = igraph::compare(cluster, true_label,method = "nmi")
    rand = igraph::compare(cluster, true_label,method = "rand")
    return(list(sil = sil, nmi = nmi,rand = rand))
  })

  return(list(compare = compare, bst_res = res_list[[which.max(sil_compare)]]))
}



#############################################
# Additional functions
#############################################
"dist2" = function( x, c = NA ) {

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


standardNormalization = function (x){
  x = as.matrix(x)
  mean = apply(x, 2, mean)
  sd = apply(x, 2, sd)
  sd[sd == 0] = 1
  xNorm = t((t(x) - mean)/sd)
  return(xNorm)
}

# using the std of data: (y-min(y))/(max(y)-min(y)) make sure entries in (0,1)
# using the defination as formula in SIMLR
D_kernels = function(x_fun = NA,Diff_fun = NA,allk_fun = 30,sigma_fun = 2,is.square = F,is.normal = T,is.SP = T,standardization = 1) {
  if(is.na(Diff_fun)){
    if(standardization == 1){
      x_fun = standardNormalization(x_fun)
    }else if(standardization == 2){
      x_fun = apply(x_fun, 2, function(y) {(y-min(y))/(max(y)-min(y))})
    }


    if(allk_fun<(nrow(x_fun))) {
      # distance calc and sort
      if(is.square){
        Diff_fun = dist2(x_fun,x_fun)^2
      }else{
        Diff_fun = dist2(x_fun,x_fun)^(1/2)
      }
    }
  }
  diag(Diff_fun) = 0
  Diff_sort_fun = t(apply(Diff_fun,MARGIN=2,FUN=sort))
  ## sort every column and transpose => ith row is the sorted distance btw ith subject and others

  # calc symmetric mu_{ij} for every subject
  TT = apply(Diff_sort_fun[,2:(allk_fun+1)],MARGIN=1,FUN=mean) + .Machine$double.eps
  ## calc mean distance for the first k neighbors of every subjects(row), length = k
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
  W_ori = D_Kernels = Matrix((W + t(W)) / 2, sparse=TRUE, doDiag=FALSE)

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
  Kernels = max(as.matrix(D_Kernels))-D_Kernels
  return(list(ori_Kernels = W_ori, Kernels = Kernels, D_Kernels = D_Kernels))
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

kernel_similarity = function(kernel_list){
  term = rep(0, length(kernel_list))
  for(i in 1:length(kernel_list)){
    for(j in 1:length(kernel_list)){
      term[i] = term[i]+sum(kernel_list[[i]]*kernel_list[[j]])
    }
  }
  return(term)
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

dist_kernels = function(kernel){
  K = kernel
  k = 1/sqrt(diag(K)+1)
  ## diagnal of the matrix is the highest similarity
  G = K * (k %*% t(k))
  ## with diag(K) all ~ 0, k is just n 1s, thus G is K
  G1 = apply(array(0,c(length(diag(G)),length(diag(G)))),MARGIN=2,FUN=function(x) {x=diag(G)})
  G2 = t(G1)
  ## Use the average difference btw 2 self similarities and pairwise similarity as kenels
  D_Kernels_tmp = (G1 + G2 - 2*G)/2
  D_Kernels = D_Kernels_tmp - diag(diag(D_Kernels_tmp))
  return(D_Kernels)
}

multi_kernels_gl = function (x, cores.ratio = 0.5,is.dist = F) {
  # data _normalization
  # x = apply(x, 2, function(y) {
  #   (y-min(y))/(max(y)-min(y))
  # })

  N = dim(x)[1]
  KK = 0
  sigma = seq(2, 1, -0.1)
  mu = seq(0.1,1, 0.1)
  allk = c(10,30)
  cores = as.integer(cores.ratio * (detectCores() - 1))
  if (cores < 1 || is.na(cores) || is.null(cores)) {
    cores = 1
  }
  cl = makeCluster(cores)
  clusterEvalQ(cl, {
    library(Matrix)
    D_kernels = function(x_fun = NA,Diff_fun = NA,allk_fun = 30,sigma_fun = 2,is.square = F,is.normal = T,is.SP = T,standardization = 1) {
      if(is.na(Diff_fun)){
        if(standardization == 1){
          x_fun = standardNormalization(x_fun)
        }else if(standardization == 2){
          x_fun = apply(x_fun, 2, function(y) {(y-min(y))/(max(y)-min(y))})
        }


        if(allk_fun<(nrow(x_fun))) {
          # distance calc and sort
          if(is.square){
            Diff_fun = dist2(x_fun,x_fun)^2
          }else{
            Diff_fun = dist2(x_fun,x_fun)^(1/2)
          }
        }
      }
      diag(Diff_fun) = 0
      Diff_sort_fun = t(apply(Diff_fun,MARGIN=2,FUN=sort))
      ## sort every column and transpose => ith row is the sorted distance btw ith subject and others

      # calc symmetric mu_{ij} for every subject
      TT = apply(Diff_sort_fun[,2:(allk_fun+1)],MARGIN=1,FUN=mean) + .Machine$double.eps
      ## calc mean distance for the first k neighbors of every subjects(row), length = k
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
      W_ori = D_Kernels = Matrix((W + t(W)) / 2, sparse=TRUE, doDiag=FALSE)

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
      Kernels = max(as.matrix(D_Kernels))-D_Kernels
      return(list(ori_Kernels = W_ori, Kernels = Kernels, D_Kernels = D_Kernels))
    }
    library(SNFtool)
  })

  D_Kernels = list()
  D_Kernels = unlist(parLapply(cl, 1:length(allk), fun = function(l, x_fun = x, B = allk,
                                                                  sig_df_list = sigma, mu_aff_list = mu, KK_fun = KK) {
    if (B[l] < (nrow(x_fun) - 1)) {
      for (j in 1:length(sig_df_list)){
        name = paste("dk_B",B[l], "_sig",sig_df_list[[j]],sep = "" )
        D_Kernels[[name]] = Matrix(D_kernels(x_fun = x,allk_fun = B[l],sigma_fun = sig_df_list[[j]],standardization = 0)[[1]])
      }

      # for(t in 1:length(mu_aff_list)){
      #   name = paste("aff_B",B[l], "_mu",mu_aff_list[[t]],sep = "" )
      #   D_Kernels[[name]] = Matrix(affinityMatrix(dist2(x,x)^(1/2),K = B[l], sigma = mu_aff_list[[t]]))
      # }

    }

    return(D_Kernels)
  }
  ))
  stopCluster(cl)
  return(D_Kernels)
}


# Distance matrix calculated from kernels
D_kernels_calc = function(kernel){
  K = kernel
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
  return(D_Kernels)
}

# umkl function
"umkl.cimlr.gl" = function( D, beta = NA ) {

  # set some parameters
  if(is.na(beta)) {
    beta = 1 / length(D)
  }
  tol = 1e-4
  u = 50
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
        beta = beta / 2
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
