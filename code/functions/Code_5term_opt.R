library(quadprog)

case_control_mat = function(response) {
  # response: a vector of binary outcomes (0s and 1s)
  # Construct a matrix to show the case/control status between responses: 1 (same) -1 (different)
  return(2*as.matrix(stats::dist(response, method = 'manhattan'))-1)
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



rho = 10^-2
alpha = 1
beta = 10^0
gamma = 10^0
