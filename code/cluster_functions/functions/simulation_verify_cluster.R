
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
library(SNFtool)
library(SIMLR)
library(abSNF)

source("./functions/Estimate_Number_of_Clusters.R")
source("./functions/SIMLR_multi_cluster.R")

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
## simulation_verify_2 add the weight changing simlr methods
simulation_verify_2data = function(K,sim,truelabel,w1 = NA,normalize = FALSE,cores.ratio = 0.5){
  data1_info = get_info(sim$data1,sim$weight1,truelabel, K = K, cores.ratio = cores.ratio,normalize = normalize)
  data2_info = get_info(sim$data2,sim$weight2,truelabel, K = K, cores.ratio = cores.ratio,normalize = normalize)

  ## for SNF method
  W = SNF(list(data1_info$af,data2_info$af),20,20)
  gap_w = estimateNumberOfClustersGivenGraph(W,2:5)$`Eigen-gap best`
  cluster_SNF = spectralClustering(W, K)
  nmi_SNF = compare(cluster_SNF, truelabel, method = "nmi")
  rand_SNF = compare(cluster_SNF, truelabel, method = "rand")

  ## for multi-SIMLR method
  D_kernels = list()
  for(i in 1:110){
    if(i<=55){
      D_kernels[[i]] = data1_info$mk[[i]]
    }else if (i<=110){
      D_kernels[[i]] = data2_info$mk[[i-55]]
    }
  }
  # for simlr_multi_SNF method
  SIMLR_MS = SIMLR_multi_SNF(D_kernels,c = K, SNF_init = W,cores.ratio = cores.ratio)
  gap_sim = estimateNumberOfClustersGivenGraph(SIMLR_MS$S, 2:5)$`Eigen-gap best`
  nmi_SIMLR = compare(SIMLR_MS$y$cluster, truelabel, method = "nmi")
  rand_SIMLR  = compare(SIMLR_MS$y$cluster, truelabel, method = "rand")
  ## for chaning weight methods
  if(!is.na(w1)){
    SIMLR_MS_cw = SIMLR_multi_SNF_cw(D_kernels,w1 = w1, w2 = 1-w1, c = K, SNF_init = W,cores.ratio = cores.ratio)
    nmi_SIMLR_cw = compare(SIMLR_MS_cw$y$cluster, truelabel, method = "nmi")
    rand_SIMLR_cw  = compare(SIMLR_MS_cw$y$cluster, truelabel, method = "rand")
  }

  ## For multi-SIMLR-weight method
  ##TO ADD

  ## get result
  if(!is.na(w1)){
    evaluation = tibble(
      cluster = K,
      gap1 = data1_info$res_data$gap,
      gap2 = data2_info$res_data$gap,
      gap_w,
      gap_sim,
      nmi1 = data1_info$res_data$nmi,
      nmi2 = data2_info$res_data$nmi,
      nmi_SNF,
      nmi_SIMLR,
      nmi_SIMLR_cw,
      rand1 = data1_info$res_data$rand,
      rand2 = data2_info$res_data$rand,
      rand_SNF,
      rand_SIMLR,
      rand_SIMLR_cw,
      weight1_sum = sum(SIMLR_MS$alphaK[1:55]),
      weight1_sum_cw = sum(SIMLR_MS_cw$alphaK[1:55]))
  }else{
    evaluation = tibble(
      cluster = K,
      gap1 = data1_info$res_data$gap,
      gap2 = data2_info$res_data$gap,
      gap_w,
      gap_sim,
      nmi1 = data1_info$res_data$nmi,
      nmi2 = data2_info$res_data$nmi,
      nmi_SNF,
      nmi_SIMLR,
      rand1 = data1_info$res_data$rand,
      rand2 = data2_info$res_data$rand,
      rand_SNF,
      rand_SIMLR,
      weight1_sum = sum(SIMLR_MS$alphaK[1:55]))
  }
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






















