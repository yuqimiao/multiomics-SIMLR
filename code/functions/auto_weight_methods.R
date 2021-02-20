# Auto-weighted multiple kernel graph learning (AMGL, in ref30)
# input:
#   1. kernel_list
#   2. c: estimated number of clusters
# output:
#   1. H: the first c eigen vectors of L_final
#   2. clusters: k-means result to rows of H

AMGL = function(kernel_list, c = 4, max = 30, tol = 1e-10){

  # Kernel normalization
  norm_kernel_list = lapply(kernel_list, FUN = function(x){
    D = diag(rowSums(x))
    solve(D^(1/2))%*%as.matrix(x)%*%solve(D^(1/2))
  }
  )


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
      w[i]=1/2*(sqrt(trace))
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
    obj = sum(sqrt(diag(t(H)%*%L%*%H)))
    criteria = c(criteria,obj)
    if(t > 10){
      if ((criteria[t]-criteria[t-1])<tol){
        cluster = kmeans(H,c,nstart = 200)$cluster
        break
      }
    }

  }
  names(w) = names(kernel_list)
  return(list(w = w,
              K = K,
              w_all = w_all,
              trace_term = trace_term,
              cluster = cluster,
              iteration = t))
}

# MVCSK: multiview clustering with single kernel
# input:
#   1. kernel_list
#   2. c: number of clusters
#   3. hyper-parameter: lambda -- for trace term,
#                       mu -- for penalty term
# output:
#   1. S: simlarity matrix
#   2. H: the first c eig vec from S
#   3. w: final weight

library(quadprog)
MVCSK = function(kernel_list = NA,
                 c = 3, # #clusters
                 lambda = 1e-3,
                 mu = 0.1,
                 stopping = 10^-3,
                 print_details = F, # whether to print stepwise details
                 n_ite = 50, # max number of iterations
                 normalize = 2){
  # normalization
  s_Kmat_ls = kernel_list # scaled kernel matrices
  if(normalize == 1){
    for (ll in 1:length(s_Kmat_ls)) {
      s_Kmat_ls[[ll]] = s_Kmat_ls[[ll]] / norm(s_Kmat_ls[[ll]], type = 'F')
    }
  }else if(normalize == 2){
    s_Kmat_ls = lapply(s_Kmat_ls, normalize_kernel)
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
    w_s_Kmat = matrix(0, n_sam, n_sam) # weighted average of scaled kernel matrices
    for (ll in 1:n_ker) {
      w_s_Kmat = w_s_Kmat + old_w[ll] * s_Kmat_ls[[ll]]
    }

    ### Initialization
    new_w = old_w
    new_L = old_L
    new_S = old_S
    d_mat = dist2(new_L)
    diag(d_mat) = 0
    linear_terms_S = (lambda/2)*d_mat - 2*w_s_Kmat

    ### Update S
    for (i in 1:n_sam) {
      QP_results = solve.QP(Dmat = 2*(I_n * mu + w_s_Kmat), # quadratic programming
                            dvec = -linear_terms_S[i,],
                            Amat =  I_n,
                            bvec =  rep(0,n_sam),
                            meq = 0)
                            # Amat = t(rbind(rep(1,n_sam), I_n)),
                            # bvec = c(1, rep(0,n_sam)),
                            # meq = 1)
      new_S[,i] = QP_results$solution
    }
    new_S = (new_S+t(new_S))/2 # make sure S is symmetric

    ### Update L
    Laplacian = I_n - new_S
    eg_results = eigen(Laplacian) # eigen-decompositions
    new_L = eg_results$vectors[,order(eg_results$values)][,1:c]# extract the two eigenvectors

    ### Update w
    sim_term = map_dbl(kernel_list,.f = function(x){sum(diag(x - 2*x%*%new_S + t(new_S)%*%x%*%new_S))})
    new_w = 1/(2*sqrt(sim_term))
    new_w = new_w/sum(new_w) # scale the kernels

    ### Print details
    if (print_details) {
      cat(paste0('Iteration ',iter, ':\n  Optimal weights: '))
      cat(new_w)
      cat('\n')
      opt_value = - sum(sim_term) +
        #-alpha * sum(new_S * s_A) +
        mu * norm(new_S, type = 'F')^2 +
        lambda * sum(diag(t(new_L) %*% (I_n - new_S) %*% new_L))
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
