setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
task_id = as.integer(Sys.getenv("SGE_TASK_ID"))
N = 1:2
n = N[task_id]
library(clValid)
library(parallel)
library("grDevices")
library("Matrix")
library("dplyr")
library("SNFtool")
library("SIMLR")
library("abSNF")

source("./functions/SIMLR_multi_cluster.R")
source("./functions/simulation_function.R")
source("./functions/simulation_verify_cluster.R")
#source("./functions/repeated_simulation.R")

## parameter setting
distribution = list(c("norm","norm"),c("norm", "logit")) # 2 
subtype_ratio = list(c(1/3,1/3,1/3),c(0.1,0.3,0.6),c(0.1,0.1,0.8)) # 3 
effective_size = list(c(1,9),c(3,7),c(5,5)) # 3
sigma = 5

par = list()
m = 1
for(k in 1:length(distribution)){
  for(i in 1:length(subtype_ratio)){
    for(j in 1:length(effective_size)){
      par[[m]] =list(sub_ratio = subtype_ratio[[i]],
                     eff_size = effective_size[[j]],
                     distribution = distribution[[k]]
      )
      m = m+1
    }
  }
}

scenario_name = paste(
  c(paste("./simulation_data_var5_20/nn",
          c(paste("sub111_",paste("eff", c("19","37","55"), sep = ""),sep = ""),
            paste("sub136_",paste("eff", c("19","37","55"), sep = ""),sep = ""),
            paste("sub118_",paste("eff", c("19","37","55"), sep = ""),sep = "")),
          sep = "_"),
    paste("./simulation_data_var5_20/nl",
          c(paste("sub111_",paste("eff", c("19","37","55"), sep = ""),sep = ""),
            paste("sub136_",paste("eff", c("19","37","55"), sep = ""),sep = ""),
            paste("sub118_",paste("eff", c("19","37","55"), sep = ""),sep = "")),
          sep = "_")),
  paste("_var",as.character(sigma), sep = ""), 
  sep = "")



sim =  simulation_2(sub_ratio = par[[n]]$sub_ratio, 
                    eff_size = par[[n]]$eff_size, 
                    dist = par[[n]]$distribution, 
                    sigma = sigma)

evaluation = simulation_verify(K = 3,sim = sim,truelabel = sim$true_label,name = scenario_name[[n]])

save(evaluation, file = paste("./evaluation_",as.character(n),".Rdata",sep = ""))

## the simlr result is in the sim_res
# The SNF matrix and D_kernels list also saved

