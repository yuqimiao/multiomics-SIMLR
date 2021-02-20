## normalization cases
# In this simulation, we want to find the effec of normalization to the results

setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
task_id = as.integer(Sys.getenv("SGE_TASK_ID"))
N = 1:4
n = N[task_id]
library("grDevices")
library("Matrix")
library("dplyr")
library("stringr")
library("SNFtool")
library("SIMLR")
library("abSNF")

source("./functions/SIMLR_multi_cluster.R")
source("./functions/simulation_function.R")
source("./functions/simulation_verify_cluster.R")

## Set the basic scnario
distribution = c("normal","normal")
sub_ratio = rep(1/3,3)
eff_size = list(c(1,9),c(3,4))
sigma = list(c(4,18^2),c(12^2,12^2))
normalization = c(TRUE, FALSE)
repeat_time = 20

par = list()
m = 1
for(i in 1:length(eff_size)){
  for(j in 1:length(normalization)){
    par[[m]] = list(distribution = distribution,
                    sub_ratio = sub_ratio,
                    eff_size = eff_size[[i]],
                    sigma = sigma[[i]],
                    normalization = normalization[[j]])
    m = m+1
  }

}

scenario_name = c(paste("nn_111_19",c("norm","noNorm"), sep = "_"),paste("nn_111_34",c("norm","noNorm"), sep = "_"))

## No normalization
table = NULL
for(i in 1:repeat_time){
  sim = simulation_2(size = 150, sub_ratio = sub_ratio, eff_size = par[[n]]$eff_size,sigma = par[[n]]$sigma)
  for(k in 2:5){
    eva_tmp = NULL
    t = 0
    while(is.null(eva_tmp) & (t <= 3)){
      tryCatch({eva_tmp = simulation_verify_2data(K = k,sim = sim,truelabel = sim$truelabel,normalize = par[[n]]$normalization)},
                error = function(e){cat("!!!!!!!!!!!!!!!!!!!",conditionMessage(e))})
      t = t+1
    }
    table = rbind(table,eva_tmp)
  }
}
table = table %>%
  mutate(scenario = str_extract(scenario_name[[n]],pattern = "nn_[0-9a-zA-Z].*")) %>%
  select(scenario, everything())



save(table, file =paste("simulation_normalization/", scenario_name[[n]],".Rdata",sep = "") )

















