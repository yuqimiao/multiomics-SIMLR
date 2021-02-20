simulation_2 = function(x_size = 150,eff_size = c(5,5), sub_ratio = rep(1/3,3),dist = c("norm", "norm"),sigma = 2,n1_r = 0.6,n2_r = 0.8,r1 = 1,r2 = 3,uninfo_r1 = 0,uninfo_r2 = 1){
  ##GE
  if(dist[1] == "norm"){
    GE = normal_sim(x_size, sub_ratio =c(sub_ratio[1],sub_ratio[2]+sub_ratio[3]),eff_size[1],sigma = sigma)
    #c(sub_ratio[1]+sub_ratio[2], sub_ratio[3])
  }else if(dist[1] == "logit"){
    GE = logit_sim(x_size, sub_ratio =c(sub_ratio[1],sub_ratio[2]+sub_ratio[3]),eff_size[1],sigma = sigma)
  }
  ##MI
  if(dist[2] == "norm"){
    MI = normal_sim(x_size, sub_ratio =c(sub_ratio[1]+sub_ratio[2], sub_ratio[3]),eff_size[2],sigma = sigma)

  }else if(dist[2] == "logit"){
    MI = logit_sim(x_size, sub_ratio =c(sub_ratio[1]+sub_ratio[2], sub_ratio[3]),eff_size[2],sigma = sigma)
  }

  ## feature weight for GE
  eff_index = c(50,50)
  ge_weight=abs(runif(1000,0,1))
  s=sample(1:eff_index[1],n1_r*eff_index[1])
  ge_weight[s]=runif(n1_r*eff_index[1],r1,r2)
  ss=sample((eff_index[1]+1):1000,n2_r*eff_index[1])
  ge_weight[ss]=runif(n2_r*eff_index[1],uninfo_r1,uninfo_r2)

  mi_weight=abs(runif(1000,0,1))
  s=sample(1:eff_index[2],n1_r*eff_index[2])
  mi_weight[s]=runif(n1_r*eff_index[2],r1,r2)
  ss=sample((eff_index[2]+1):1000,n2_r*eff_index[2])
  mi_weight[ss]=runif(n2_r*eff_index[2],uninfo_r1,uninfo_r2)
  a = sub_ratio[1]*x_size
  b = (sub_ratio[1]+sub_ratio[2])*x_size
  true_label = c(rep(1,a),rep(2,(b-a)),rep(3,x_size-b))

  return(list(GE = GE, MI = MI,ge_weight = ge_weight, mi_weight = mi_weight,true_label = true_label))
}

normal_sim = function(x_size = 150, sub_ratio = c(2/3,1/3), eff_size = 5,sigma = 2){
  data1 = matrix(0, x_size, 1000)
  a = sub_ratio[1]*x_size
  sub1_index = 1:a
  sub2_index = (a+1):x_size
  data1[sub1_index,1:50] = rnorm(a*50, -eff_size, sqrt(sigma))
  data1[sub2_index,1:50] = rnorm((x_size-a)*50, eff_size, sqrt(sigma))
  data1[,51:1000] = rnorm(x_size*950, 0, sqrt(sigma))
  return(data1)
}

beta_sim = function(x_size = 150, sub_ratio = c(2/3,1/3), par = list(c(1,5),c(5,1)),sigma = 2){
  data1 = matrix(0, x_size, 1000)
  a = sub_ratio[1]*x_size
  sub1_index = 1:a
  sub2_index = (a+1):x_size
  data1[sub1_index,1:50] = rbeta(a*50, shape1 = par[[1]][1], shape2 = par[[1]][2]) ## hypo
  data1[sub2_index,1:50] = rbeta((x_size-a)*50, shape1 = par[[2]][1], shape2 = par[[2]][2])
  data1[,51:1000] = rbeta(x_size*950, 5,5, sqrt(sigma))
  return(data1)
}

logit_sim = function(x_size = 150, sub_ratio = c(2/3,1/3), eff_size = 5,sigma = 2){
  data1 = matrix(0, x_size, 1000)
  a = sub_ratio[1]*x_size
  sub1_index = 1:a
  sub2_index = (a+1):x_size
  data1[sub1_index,1:50] = sapply(rnorm(a*50, -eff_size, sqrt(sigma)), FUN = function(x) exp(x)/(1+exp(x)))
  data1[sub2_index,1:50] = sapply(rnorm((x_size-a)*50, eff_size, sqrt(sigma)),FUN = function(x) exp(x)/(1+exp(x)))
  data1[,51:1000] = sapply(rnorm(x_size*950, 0, sqrt(sigma)),FUN = function(x) exp(x)/(1+exp(x)))
  return(data1)
}

simulation_verify = function(K,sim,truelabel,name,normalize = FALSE,cores.ratio = 0.5,save = FALSE){

  Data1 = sim$GE
  Data2 = sim$MI
  if(normalize == TRUE){
    Data1 = standardNormalization(Data1)
    Data2 = standardNormalization(Data2)
  }
  weight1 = sim$ge_weight
  weight2 = sim$mi_weight

  ## normalization
  # Data1 = standardNormalization(Data1)
  # Data2 = standardNormalization(Data2)


  ### use the weighted kernel with SP to construct SNF

  mk1 = multiple.kernel(Data1,cores.ratio = cores.ratio)
  d1 = matrix(0,150,150)
  for(i in 1:length(mk1)){
    d1 = d1+as.matrix(mk1[[i]])
  }
  d1 = d1/length(mk1)
  k1 = max(max(d1))-d1
  gap1 = estimateNumberOfClustersGivenGraph(k1, 2:5)$`Eigen-gap best`

  mk2 = multiple.kernel(Data2,cores.ratio = cores.ratio)
  d2 = matrix(0,150,150)
  for (i in 1:length(mk2)){
    d2 = d2+as.matrix(mk2[[i]])
  }
  d2 = d2/length(mk2)
  k2 = max(max(d2)) - d2
  gap2 = estimateNumberOfClustersGivenGraph(k2, 2:5)$`Eigen-gap best`

  ## calculate SIMLR for 2 data types separately
  simlr1 = SIMLR(t(Data1), c = K,cores.ratio = cores.ratio)
  simlr2= SIMLR(t(Data2), c = K,cores.ratio = cores.ratio)

  nmi1 = compare(simlr1$y$cluster,truelabel,method = "nmi")
  nmi2 = compare(simlr2$y$cluster,truelabel,method = "nmi")
  rand1 = compare(simlr1$y$cluster,truelabel,method = "rand")
  rand2 = compare(simlr2$y$cluster,truelabel,method = "rand")

  ## optimal kernels
  # optimal_mk1 = as.matrix(mk1[[which.max(simlr1$alphaK)]])
  # optimal_mk2 = as.matrix(mk2[[which.max(simlr2$alphaK)]])
  D_kernels = list()
  for(i in 1:(length(mk1)+length(mk2))){
    if(i<=55){
      D_kernels[[i]] = mk1[[i]]
    }else{
      D_kernels[[i]] = mk2[[i-length(mk1)]]
    }
  }
  # weighted affinity matrix
  W1 = affinityMatrix(dist2_w(as.matrix(Data1),as.matrix(Data1),weight = weight1)^(1/2))
  W2 = affinityMatrix(dist2_w(as.matrix(Data2),as.matrix(Data2),weight = weight2)^(1/2))
  ## SNF

  # W = SNF(list(optimal_mk1,optimal_mk2), 20 ,10)
  W = SNF(list(W1,W2),20,20)
  gap_w = estimateNumberOfClustersGivenGraph(W,2:5)$`Eigen-gap best`
  cluster_SNF = spectralClustering(W, K)
  nmi_SNF = compare(cluster_SNF, truelabel, method = "nmi")
  rand_SNF = compare(cluster_SNF, truelabel, method = "rand")

  ## SIMLR_MS(multi_SNF)
  #SIMLR_MS = SIMLR_multi_SNF(list(optimal_mk1,optimal_mk2),c = K, SNF_init = W)
  SIMLR_MS = SIMLR_multi_SNF(D_kernels,c = K, SNF_init = W,cores.ratio = cores.ratio)
  gap_sim = estimateNumberOfClustersGivenGraph(SIMLR_MS$S, 2:5)$`Eigen-gap best`
  nmi_SIMLR = compare(SIMLR_MS$y$cluster, truelabel, method = "nmi")
  rand_SIMLR  = compare(SIMLR_MS$y$cluster, truelabel, method = "rand")

  evaluation = tibble(cluster = K, gap1,gap2,gap_w,gap_sim,nmi1,nmi2,nmi_SNF,nmi_SIMLR,rand1,rand2,rand_SNF,rand_SIMLR,GE_weightsum = sum(SIMLR_MS$alphaK[1:55]))
  kernel_list = list(mk1 = mk1,mk2 = mk2)
  sim_res = list(S1 = simlr1, S2 = simlr2,S_MS = SIMLR_MS)
  if(save){
    save(W,file = paste(name,"_",as.character(K),"_SNF.Rdata",sep = ""))
    save(D_kernels,file =paste(name,"_",as.character(K),"_kernel_list.Rdata",sep = ""))
    save(sim_res,file =paste(name,"_",as.character(K),"_Simlr_result.Rdata",sep = ""))
  }
  return(evaluation)
}
