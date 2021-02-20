## See how normalization influence the global-local term
# data_list = list(t(data$data1),t(data$data2),t(data$data3))
# D_list = NA

kernel_validation = function(data_list, B = 30, D_list = NA, is.dnorm = T, normal_type = 2){
  ## calculate global kernel
  kernel_list = list() ## the kerenl matrix, using the density of normal distribution
  if(is.na(D_list)){
    D_list = list() ## using kernel to calculate distance here: Dij = kii+kjj-2kij
    for(i in 1:length(data_list)){
      kernel_list[[i]] = kernel_calculating(data_list[[i]], is.SP = F, is.normal = is.dnorm)
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
  ## calculate local kernel
  local_list = NULL
  for(i in 1:length(kernel_list)){
    local_list[[i]] = dominateset(kernel_list[[i]],KK = B)
  }


  ## 3 different normalization
  for(i in 1:length(kernel_list)){
    if(normal_type == 0){
      ## no normalization
      kernel_list = kernel_list
      local_list = local_list
    }else if(normal_type == 1){
      kernel_list[[i]] = normalize_kernel(kernel_list[[i]])
      local_list[[i]] = normalize_kernel(local_list[[i]])
    }else if(normal_type == 2){
      D = diag(rowSums(kernel_list[[i]]))
      kernel_list[[i]] = solve(D^(1/2)) %*% kernel_list[[i]] %*% solve(D^(1/2))
      Dl = diag(rowSums(local_list[[i]]))
      local_list[[i]] = solve(Dl^(1/2)) %*% local_list[[i]] %*% solve(Dl^(1/2))
    }else{
      kernel_list[[i]] = kernel_list[[i]]/sum(kernel_list[[i]]^2)
      local_list[[i]] = local_list[[i]]/sum(local_list[[i]]^2)
    }
  }



  ## Calculate kernel terms
  global_F = rep(0,length(kernel_list))
  local_F = rep(0,length(kernel_list))
  global_F_ori = rep(0,length(kernel_list))
  local_F_ori = rep(0,length(kernel_list))
  for (i in 1:length(kernel_list)){
    for (j in 1:length(kernel_list)){
      global_F[i] = global_F[i] + sum(kernel_list[[i]]*kernel_list[[j]])
      local_F[i] = local_F[i] + sum(local_list[[i]]*local_list[[j]])
    }
  }

  global_local_term = tibble(global = global_F , local = local_F)
  return(global_local_term = global_local_term)
}




# Initial conclusion: the kernel normalization will change













