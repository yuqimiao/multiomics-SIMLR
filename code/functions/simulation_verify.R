library(SNFtool)
library(SIMLR)
library(abSNF)

source("functions/Estimate_Number_of_Clusters.R")
source("functions/SIMLR_multi.R")

## simulation function
# input:
#   scenario parameter -- including subtype ratio:(1/3,1/3,1/3), (0.1,0.3,0.6),(0.1,0.1,0.8), 
#                                   effective size: (5,5), (1,10)
#   K -- number of clusters to tunn
# output:
#   simulation table for 1,2,SNF, SIMLR:
#   gap -- the eigen gap for corresponding similarity matrix
#   NMI -- normalized mutual information for spectual clustering results from corresponding similarity matrix
#   rand -- rand index for spectual clustering results from corresponding similarity matrix


### use the affinity matrix def in SNF to construct SNF
# set.seed(123)
# sim = simulation_1() ## balanced case
# K = 3
# truelabel = sim$true_label


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

library(SNFtool)
library(SIMLR)
library(abSNF)

source("functions/Estimate_Number_of_Clusters.R")
source("functions/SIMLR_multi.R")

## simulation function( 3 type)
# input:
#   scenario parameter -- including subtype ratio:(1/3,1/3,1/3), (0.1,0.3,0.6),(0.1,0.1,0.8), 
#                                   effective size: (5,5), (1,10)
#   K -- number of clusters to tunn
# output:
#   simulation table for 1,2,SNF, SIMLR:
#   gap -- the eigen gap for corresponding similarity matrix
#   NMI -- normalized mutual information for spectual clustering results from corresponding similarity matrix
#   rand -- rand index for spectual clustering results from corresponding similarity matrix


### use the affinity matrix def in SNF to construct SNF


simulation_verify_3data = function(K,sim,truelabel,normalize = FALSE,cores.ratio = 0.5){
  data1_info = get_info(sim$data1,sim$weight1,truelabel, K = K, cores.ratio = cores.ratio,normalize = normalize)
  data2_info = get_info(sim$data2,sim$weight2,truelabel, K = K, cores.ratio = cores.ratio,normalize = normalize)
  data3_info = get_info(sim$data3,sim$weight3,truelabel, K = K, cores.ratio = cores.ratio,normalize = normalize)
  
  ## for SNF method
  W = SNF(list(data1_info$af,data2_info$af,data3_info$af),20,20)
  gap_w = estimateNumberOfClustersGivenGraph(W,2:5)$`Eigen-gap best`
  cluster_SNF = spectralClustering(W, K)
  nmi_SNF = compare(cluster_SNF, truelabel, method = "nmi")
  rand_SNF = compare(cluster_SNF, truelabel, method = "rand")
  
  ## for multi-SIMLR method
  D_kernels = list()
  for(i in 1:165){
    if(i<=55){
      D_kernels[[i]] = data1_info$mk[[i]]
    }else if (i<=110){
      D_kernels[[i]] = data2_info$mk[[i-55]]
    }else{
      D_kernels[[i]] = data3_info$mk[[i-110]]
    }
  }
  SIMLR_MS = SIMLR_multi_SNF(D_kernels,c = K, SNF_init = W,cores.ratio = cores.ratio)
  gap_sim = estimateNumberOfClustersGivenGraph(SIMLR_MS$S, 2:5)$`Eigen-gap best`
  nmi_SIMLR = compare(SIMLR_MS$y$cluster, truelabel, method = "nmi")
  rand_SIMLR  = compare(SIMLR_MS$y$cluster, truelabel, method = "rand")
  
  ## For multi-SIMLR-weight method
  ##TO ADD
  
  ## get result
  evaluation = tibble(
    cluster = K, 
    gap1 = data1_info$res_data$gap,
    gap2 = data2_info$res_data$gap,
    gap3 = data3_info$res_data$gap,
    gap_w,
    gap_sim,
    nmi1 = data1_info$res_data$nmi,
    nmi2 = data2_info$res_data$nmi,
    nmi3 = data3_info$res_data$nmi,
    nmi_SNF,
    nmi_SIMLR,
    rand1 = data1_info$res_data$rand,
    rand2 = data2_info$res_data$rand,
    rand3 = data3_info$res_data$rand,
    rand_SNF,
    rand_SIMLR,
    weight1_sum = sum(SIMLR_MS$alphaK[1:55]),
    weight2_sum = sum(SIMLR_MS$alphaK[56:110])
    )
  return(evaluation)
}

get_info = function(data1,weight1,truelabel,K = 2, cores.ratio = 0.5,normalize = T){
  if(normalize){
    data1 = standardNormalization(data1)
    }
  mk = multiple.kernel(data1)
  af = affinityMatrix(dist2_w(as.matrix(data1),as.matrix(data1), weight = weight1))
  gap = estimateNumberOfClustersGivenGraph(af,2:5)$`Eigen-gap best`
  simlr1 = SIMLR(t(data1), K,cores.ratio = cores.ratio)
  nmi = compare(simlr1$y$cluster,truelabel, method = "nmi")
  rand = compare(simlr1$y$cluster,truelabel, method = "rand")
  res_data = tibble(gap = gap, nmi = nmi, rand = rand)
  return(list(res_data = res_data, mk = mk, af = af))
}





















