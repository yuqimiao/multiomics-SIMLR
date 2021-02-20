library(quadprog)

case_control_mat = function(response) {
  # response: a vector of binary outcomes (0s and 1s)
  # Construct a matrix to show the case/control status between responses: 1 (same) -1 (different)
  return(((1-as.matrix(stats::dist(response, method = 'manhattan')))-0.5)*2)
}

Crazy_5term_opt = function(Kmat_ls, # a list of kernel matrices for all samples
                           response, # binary outcomes
                           rho, # tuning parameter for entropy penalty
                           alpha, # tuning parameter for <S, A/||A||>
                           beta, # tuning parameter for ||S||
                           gamma, # tuning parameter for laplacian
                           stopping = 10^-3, # stopping rule
                           n_ite = 50, # max number of iterations
                           print_details = F # whether to print stepwise details
)
{
  s_Kmat_ls = Kmat_ls # scaled kernel matrices
  for (ll in 1:length(s_Kmat_ls)) {
    s_Kmat_ls[[ll]] = s_Kmat_ls[[ll]] / norm(s_Kmat_ls[[ll]], type = 'F')
  }
  n_sam = nrow(s_Kmat_ls[[1]]) # number of samples
  n_ker = length(s_Kmat_ls) # number of kernels
  I_n = diag(n_sam) # diagonal matrix
  old_w = rep(1/n_ker, n_ker) # initialize weights
  old_S = matrix(NA, nrow = n_sam, ncol = n_sam) # initialize S
  old_L = matrix(0, nrow = n_sam, ncol = 2) # initialize L
  old_L[which(response == 1), 1] = old_L[which(response == 0), 2] = 1
  old_L = t(t(old_L)/sqrt(colSums(old_L))) # scale old_L
  A = case_control_mat(response) # case control status matrix
  s_A = A / norm(A, type = 'F') # scaled A

  ### Start iteration
  for (k in 1:n_ite) {
    w_s_Kmat = matrix(0, n_sam, n_sam) # weighted average of scaled kernel matrices
    for (ll in 1:n_ker) {w_s_Kmat = w_s_Kmat + old_w[ll] * s_Kmat_ls[[ll]]}

    ### Initialization
    new_w = old_w
    new_L = old_L
    new_S = old_S
    linear_terms_S = w_s_Kmat + alpha * s_A + gamma * (old_L %*% t(old_L))

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
    new_L = eg_results$vectors[, c(n_sam-1, n_sam)] # extract the two eigenvectors

    ### Update w
    first_term = vector() # 1st term in optimization
    for (ll in 1:n_ker) {
      first_term[ll] = sum(s_Kmat_ls[[ll]] * new_S)
      new_w[ll] = exp(first_term[ll] / rho) # the new weights for kernels
    }
    new_w = new_w/sum(new_w) # scale the kernels

    ### Print details
    if (print_details) {
      cat(paste0('Iteration ',k, ':\n  Optimal weights: '))
      cat(new_w)
      cat('\n')
      opt_value = - sum(first_term) +
        rho * sum(new_w * log(new_w)) -
        alpha * sum(new_S * s_A) +
        beta * norm(new_S, type = 'F')^2 +
        gamma * sum(diag(t(new_L) %*% (I_n - new_S) %*% new_L))
      cat(paste0('  Criterion: ', opt_value, '\n'))
    }

    ### Whether to stop
    if (k>=3 & max(abs(new_w-old_w))<= stopping) {break} else {
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

## changing to cluster (unsupervised)
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

Crazy_5term_clust = function(data_list, # a list of data matrices to integrate
                           k = 30, # #neighbors
                           c = 3, # #clusters
                           kernel_fun = "D_kernels", # the definition of kernels
                           rho, # tuning parameter for entropy penalty
                           alpha0,# tuning parameter for <S, K>
                           alpha, # tuning parameter for <S, W>
                           beta, # tuning parameter for ||S||
                           gamma, # tuning parameter for laplacian
                           stopping = 10^-3, # stopping rule
                           n_ite = 50, # max number of iterations
                           normalize = T,
                           print_details = F # whether to print stepwise details
)
{
  if(kernel_fun == "D_kernels"){
    Kmat_ls = lapply(data_list, function(x) as.matrix(D_kernels(x_fun = x)$Kernels))
  }else if(kernel_fun == "affinity"){
    Kmat_ls = lapply(data_list, function(x) affinityMatrix(dist2(x)^(1/2)))
  }

  s_Kmat_ls = Kmat_ls # scaled kernel matrices
  s_Qmat_ls = lapply(s_Kmat_ls, FUN = function(x) dominateset(x, KK = k))
  if(normalize == 1){
    for (ll in 1:length(s_Kmat_ls)) {
      s_Kmat_ls[[ll]] = s_Kmat_ls[[ll]] / norm(s_Kmat_ls[[ll]], type = 'F')
      s_Qmat_ls[[ll]] = s_Qmat_ls[[ll]] / norm(s_Qmat_ls[[ll]], type = 'F')
    }
  }else{
    s_Kmat_ls = lapply(s_Kmat_ls, normalize_kernel)
    s_Qmat_ls = lapply(s_Qmat_ls, normalize_kernel)
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
  for (k in 1:n_ite) {
    w_s_Kmat = w_s_Qmat = matrix(0, n_sam, n_sam) # weighted average of scaled kernel matrices
    for (ll in 1:n_ker) {
      w_s_Kmat = w_s_Kmat + old_w[ll] * s_Kmat_ls[[ll]]
      w_s_Qmat = w_s_Qmat + old_w[ll] * s_Qmat_ls[[ll]]
      }

    ### Initialization
    new_w = old_w
    new_L = old_L
    new_S = old_S
    linear_terms_S = alpha0*w_s_Kmat + alpha * w_s_Qmat + gamma * (old_L %*% t(old_L))

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
    first_term = vector() # 1st term in optimization
    for (ll in 1:n_ker) {
      first_term[ll] = sum((s_Kmat_ls[[ll]]+alpha*s_Qmat_ls[[ll]]) * new_S)
      new_w[ll] = exp(first_term[ll] / rho) # the new weights for kernels
    }
    new_w = new_w/sum(new_w) # scale the kernels

    ### Print details
    if (print_details) {
      cat(paste0('Iteration ',k, ':\n  Optimal weights: '))
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
    if (k>=3 & max(abs(new_w-old_w))<= stopping) {break} else {
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



rho = 10^-2
alpha = 1
beta = 10^0
gamma = 10^0
