## get the repeat ID from cluster
setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
task_ID = as.integer(Sys.getenv("SGE_TASK_ID"))
N= 1:4
n = N[task_ID]

## dependencies
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
library("grDevices")
library("Matrix")
library("dplyr")
library("stringr")
library("SNFtool")
library("SIMLR")
library("abSNF")
library("parallel")
library(gplots)
library(igraph)
library("clValid")

source("./functions/simulation_function.R")
source("./functions/SIMLR_Changing_rho.R")

source("./R/SIMLR.R")
source("./R/compute.multiple.kernel.R")
source("./R/network.diffusion.R")
source("./R/utils.simlr.R")
source("./R/tsne.R")

dyn.load("./R/projsplx_R.so")

## set the simulation scenarios and parameters
eff_size = list(c(5,5),c(3,7))
sub_ratio = list(rep(1/3,3),c(0.3,0.4,0.3))
w_data = mapply(c, seq(0.1,0.9,by = 0.2), seq(0.9,0.1, by = -0.2),SIMPLIFY = F)

par = list()
m = 1
for (j in 1:length(sub_ratio)){
  for (i in 1:length(eff_size)){
    par[[m]] = list("sub_ratio" = sub_ratio[[j]],
                    "eff_size" = eff_size[[i]]
    )
    m = m+1
  }
}

scenario_name = c(paste("sub111_",c("eff55","eff37"),sep = ""),
                  paste("sub343_",c("eff55","eff37"),sep = ""))

repeat_time = 20

## function to get the multiple_kernel, affinity matrix from data
get_info = function(data1,weight1,cores.ratio = 0.5,normalize = T){
  if(normalize){
    data1 = standardNormalization(data1)
  }
  mk = multiple.kernel(data1,cores.ratio = cores.ratio)
  simlr1 = SIMLR(t(data1),2,cores.ratio = cores.ratio)
  af = affinityMatrix(dist2_w(as.matrix(data1),as.matrix(data1), weight = weight1))
  opt_w = simlr1$alphaK
  return(list(opt_w = opt_w, mk = mk, af = af))
}

## repeat simulation
nmi_all = NULL
cluster_all = NULL
for (j in 1:repeat_time){
  sim = simulation_2(sub_ratio = par[[n]]$sub_ratio,
                     eff_size = par[[n]]$eff_size,
                     sigma = (2*(par[[n]]$eff_size)^2))
  truelabel = sim$truelabel
  ## get data info
  data1_info = get_info(sim$data1,sim$weight1)
  data2_info = get_info(sim$data2,sim$weight2)

  ## calculate SNF
  W = SNF(list(data1_info$af,data2_info$af),20,20)
  cluster_SNF = spectralClustering(W,3)
  ## calculate multiple_kernels
  D_kernels = list()
  for(i in 1:110){
    if(i<=55){
      D_kernels[[i]] = data1_info$mk[[i]]
    }else if (i<=110){
      D_kernels[[i]] = data2_info$mk[[i-55]]
    }
  }

  ## get the weight list
  weight_list = list(opt_w1 = data1_info$opt_w, opt_w2 = data2_info$opt_w)

  ## fit for every weight combination
  nmi_fixed = NULL
  cluster = NULL
  for(i in 1:length(w_data)){
    res_f = SIMLR_multi_fix_weight(D_Kernels = D_kernels,c = 3, weight_list = weight_list, w_data = w_data[[i]])
    cluster = rbind(cluster,res_f$y$cluster)
    nmi_fixed = c(nmi_fixed, compare(res_f$y$cluster, truelabel, method = "nmi"))
  }

  names(nmi_fixed) = paste("w1_",seq(1,9,by = 2),sep = "")
  colnames(cluster) = paste("C",1:150, sep = "")
  nmi_fixed = c(nmi_fixed,nmi_SNF = compare(cluster_SNF,truelabel, method = "nmi"))
  ## write in the final table
  nmi_all = rbind(nmi_all, nmi_fixed)
  cluster_all = rbind(cluster_all,as_tibble(cluster))
}


nmi_all = as_tibble(nmi_all) %>%
  mutate(scenario = scenario_name[n]) %>%
  select(scenario,everything())

cluster_all = as_tibble(cluster_all) %>%
  mutate(scenario = scenario_name[n]) %>%
  select(scenario,everything())

saveRDS(nmi_all,file = paste("simulation_fixWeight/",scenario_name[n],"_nmi.rds",sep = ""))
saveRDS(cluster_all,file = paste("simulation_fixWeight/",scenario_name[n],"_cluser.rds",sep = ""))
