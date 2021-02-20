# Scenario:
# * eff_size: c(5,5),c(3,7)
# * sub_ratio: rep(1/3,3),c(0.3,0.4,0.3)

# Under every scenario, we test U = c(20,30,40,50) separately,
# the focus is the sum of weight for first data type and the final preformance.

## get the repeat ID from cluster
setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
task_ID = as.integer(Sys.getenv("SGE_TASK_ID"))
N= 1:16
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
u_collect = seq(20,50,by = 10)

par = list()
m = 1
for (j in 1:length(sub_ratio)){
  for (i in 1:length(eff_size)){
    for (k in 1:length(u_collect)){
      par[[m]] = list("u" = u_collect[[k]],
                      "sub_ratio" = sub_ratio[[j]],
                      "eff_size" = eff_size[[i]]
                      )
      m = m+1
    }
  }
}

scenario_name = c(paste("sub111_",
                        c(paste("eff55_", paste("u",as.character(seq(20,50,by = 10)),sep = ""),sep = ""),
                          paste("eff37_", paste("u",as.character(seq(20,50,by = 10)),sep = ""),sep = "")),
                        sep = ""),
                  paste("sub343_",
                        c(paste("eff55_", paste("u",as.character(seq(20,50,by = 10)),sep = ""),sep = ""),
                          paste("eff37_", paste("u",as.character(seq(20,50,by = 10)),sep = ""),sep = "")),
                        sep = ""))
repeat_time = 20

## function to get the multiple_kernel, affinity matrix from data
get_info = function(data1,weight1,cores.ratio = 0.5,normalize = T){
  if(normalize){
    data1 = standardNormalization(data1)
  }
  mk = multiple.kernel(data1,cores.ratio = cores.ratio)
  af = affinityMatrix(dist2_w(as.matrix(data1),as.matrix(data1), weight = weight1))
  gap = estimateNumberOfClustersGivenGraph(af,2:5)$`Eigen-gap best`
  return(list(gap = gap, mk = mk, af = af))
}

## trycatch
table = NULL
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

  ## calculate multiple_kernels
  D_kernels = list()
  for(i in 1:110){
    if(i<=55){
      D_kernels[[i]] = data1_info$mk[[i]]
    }else if (i<=110){
      D_kernels[[i]] = data2_info$mk[[i-55]]
    }
  }
  ## for every # clusters
  for(k in 2:4){

    print(paste("start scenario:", scenario_name[[n]], as.character(k)))
    print(paste("repeat time:",as.character(j)))

    ##Get the SNF result
    gap_w = estimateNumberOfClustersGivenGraph(W,2:5)$`Eigen-gap best`
    cluster_SNF = spectralClustering(W, k)
    nmi_SNF = compare(cluster_SNF, truelabel, method = "nmi")
    rand_SNF = compare(cluster_SNF, truelabel, method = "rand")

    ## Get the SIMLR result
    sim_res = NULL
    t = 1
    # using trycatch to avoid eruption
    if(is.null(sim_res) & t <= 3){
      tryCatch({sim_res = SIMLR_multi_u(D_Kernels = D_kernels, c = k, u = par[[n]]$u, cores.ratio = 0.5)},
               error = function(e){cat("!!!!!!!!!!!!!!!!!!!!!!!!",conditionMessage(e))})
    }
    if(!is.null(sim_res)){
      data1_weight = sum(sim_res$alphaK[1:55])
      gap_sim = estimateNumberOfClustersGivenGraph(sim_res$S, 2:5)$`Eigen-gap best`
      nmi_SIMLR = compare(sim_res$y$cluster, truelabel, method = "nmi")
      rand_SIMLR  = compare(sim_res$y$cluster, truelabel, method = "rand")

      ## write the table
      evaluation = tibble(
        cluster = k,
        gap1 = data1_info$gap,
        gap2 = data2_info$gap,
        gap_w,
        gap_sim,
        nmi_SNF,
        nmi_SIMLR,
        rand_SNF,
        rand_SIMLR,
        data1_weight = data1_weight)
      table = rbind(table,evaluation)
    }
  }
}

table = table %>%
  mutate(scenario = scenario_name[n]) %>%
  select(scenario,everything())

saveRDS(table,file = paste("simulation_rho/",scenario_name[n],".rds",sep = ""))














