
### use the affinity matrix def in SNF to construct SNF
library(SNFtool)
library(SIMLR)
library(abSNF)

source("R/SIMLR_Estimate_Number_of_Clusters.R")

# set.seed(1)
# simulation1 = simulation_1()
# K = 3 # total # clusters, if only 1 type, using K-1
# ## TODO need to tuning by eigen gap
# alpha = 0.5
# iter = 10
# Data1 = simulation1$GE # 150 samples * 1000 features 
# Data2 = simulation1$MI
# truelabel = simulation1$true_label
# k  = seq(10,30,2)
# sig = seq(2,1,-0.25)
# kernel = list()
# for(i in 1:length(k)){
#   for(j in 1:length(sig)){
#     kernel[[(i-1)*length(sig)+j]] = c(k[i],sig[j] )
#   }
# }
SNF_initial_generater = function(K, alpha,iter,Data1,Data2,truelabel){
  ## normalization
  Data1 = standardNormalization(Data1)
  Data2 = standardNormalization(Data2)
  
  ## calc distance
  Dist1 = dist2(as.matrix(Data1),as.matrix(Data1))
  Dist2 = dist2(as.matrix(Data2),as.matrix(Data2))
  
  ## calc similarity matrix
  W1 = affinityMatrix(Dist1, K, alpha)
  W2 = affinityMatrix(Dist2, K, alpha)
  
  ## SNF construction
  W = SNF(list(W1,W2), K, iter)
  
  ## spectual clustering based on W
  cluster_SNF = spectralClustering(W,K,type = 2)
  nmi_SNF = compare(cluster_SNF, truelabel,method = "nmi")
  rand_SNF = compare(cluster_SNF,truelabel, method = "rand")
  
  
  ### use the weighted kernel with SP to construct SNF
  
  simlr1 = SIMLR(t(Data1), c = K-1)
  simlr2= SIMLR(t(Data2), c = K-1)
  
  
  ## all 55 multiple kernels
  mk1 = multiple.kernel.verify(t(Data1),allk = seq(10,30,2),sigma = seq(2,1,-0.25))
  mk2 = multiple.kernel.verify(t(Data2),allk = seq(10,30,2),sigma = seq(2,1,-0.25))
  ## why can't assign outside?
  
  ## optimal kernels
  optimal_mk1 = as.matrix(mk1[[which.max(simlr1$alphaK)]])
  optimal_mk2 = as.matrix(mk1[[which.max(simlr2$alphaK)]])
  
  ## weighted kernels
  weighted_mk1 = 0
  weighted_mk2 = 0
  for(i in 1:length(mk1)){
    weighted_mk1 = weighted_mk1 + simlr1$alphaK[i]*mk1[[i]]
    weighted_mk2 = weighted_mk2 + simlr2$alphaK[i]*mk2[[i]]
  }
  
  weighted_mk1 = as.matrix(weighted_mk1)
  weighted_mk2 = as.matrix(weighted_mk2)
  
  W_opt = SNF(list(optimal_mk1,optimal_mk2), K,iter)
  W_weight = SNF(list(weighted_mk1,weighted_mk2), K,iter)
  
  cluster_SNF_opt = spectralClustering(W_opt, K)
  cluster_SNF_weight = spectralClustering(W_weight, K)
  
  nmi_SNF_opt = compare(cluster_SNF_opt, truelabel, method = "nmi")
  nmi_SNF_weight = compare(cluster_SNF_weight, truelabel, method = "nmi")
  
  rand_SNF_opt = compare(cluster_SNF_opt, truelabel, method = "rand")
  rand_SNF_weight = compare(cluster_SNF_weight, truelabel, method = "rand") 
  
  evaluation = 
    tibble(nmi = c(nmi_SNF,nmi_SNF_opt,nmi_SNF_weight),
           rand = c(rand_SNF,rand_SNF_opt,rand_SNF_weight)) %>% 
    mutate(type = c("SNF","SNF_opt","SNF_weight"))
  SNF_list = list(W = W, W_opt = W_opt, W_weight = W_weight)
  kernel_list = list(mk = list(mk1 = mk1,mk2 = mk2),
                     aff = list(W1,W2),
                     opt_mk = list(optimal_mk1 = optimal_mk1,optimal_mk2 = optimal_mk2), 
                     weight_mk = list(weighted_mk1 = weighted_mk1,weighted_mk2 = weighted_mk2))
  sim_res = list(S1 = simlr1, S2 = simlr2)
  
  return(list(evaluation = evaluation,SNF_list = SNF_list,kernel_list = kernel_list,sim_res = sim_res))
}


## kernel matrix choosen
#Goal: Use simulation to see which similarity definition works better

# kernel_choosen = function(n = 10){
#   
# }










