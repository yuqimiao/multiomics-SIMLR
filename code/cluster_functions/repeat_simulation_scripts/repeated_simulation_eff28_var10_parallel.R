## 1 time simulation
# input:
#   simulation scenario:
#     sub_ratio
#     eff_size
#     distrbution
#     scenario_name
#   repeat_time
#     
# Output:
#   stored_data: 
#     D_kernels:110 kernels
#     W: SNF matrix
#     sim_res: sim1,sim2,sim_multi
#   table_ave: averaged measures
#   table_std: standard derivation

setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")  
task_ID = as.integer(Sys.getenv("SGE_TASK_ID"))
N= 1:18
n = N[task_ID]
library("grDevices")
library("Matrix")
library("dplyr")
library("SNFtool")
library("SIMLR")
library("abSNF")

source("./functions/SIMLR_multi_cluster.R")
source("./functions/simulation_function.R")
source("./functions/simulation_verify_cluster.R")

simulation_repeat = function(sub_ratio, eff_size, distribution, scenario_name, repeat_time,sigma = 2){
  normalize = FALSE
  if(length(unique(distribution))!=1){
    normalize = TRUE
  }
  table = NULL
  for(i in 1: repeat_time){
    sim = simulation_2(sub_ratio = sub_ratio, eff_size = eff_size, dist = distribution,sigma = sigma)
    truelabel = sim$true_label
    save = FALSE
    if(i == repeat_time){
      save = TRUE
    }
    for (K in 2:5){
      tmp1 = NULL
      print(paste("start scenario:", scenario_name, as.character(K)))
      print(paste("repeat time:",as.character(i)))
      t = 0
      while(is.null(tmp1) & t<=3){
        tryCatch({tmp1 = simulation_verify(K,sim,truelabel,normalize = normalize,name = scenario_name,save = save)},
                 error = function(e){cat("!!!!!!!!!!!!!!!!!!!!!!!!",conditionMessage(e))})
        t = t+1
      }
      table = rbind(table, tmp1)
    }
  }
  
  table_ave = table %>% group_by(cluster)%>% summarise_all(mean)
  table_std = table %>% group_by(cluster)%>% summarise_all(sd)
  save(table, file = paste(scenario_name, "_all.Rdata",sep = ""))
  save(table_ave, file = paste(scenario_name, "_mean.Rdata",sep = ""))
  save(table_std, file = paste(scenario_name, "_std.Rdata",sep = ""))
  return(0)
}

## parameter setting
distribution = list(c("norm","norm"),c("norm", "logit")) # 2 
subtype_ratio = list(c(1/3,1/3,1/3),c(0.1,0.3,0.6),c(0.1,0.1,0.8)) # 3 
effective_size = list(c(2,8)) # 1
sigma = 10

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
  c(paste("./simulation_data_var10_50/nn",
          c(paste("sub111_",paste("eff", "28", sep = ""),sep = ""),
            paste("sub136_",paste("eff", "28", sep = ""),sep = ""),
            paste("sub118_",paste("eff", "28", sep = ""),sep = "")),
          sep = "_"),
    paste("./simulation_data_var10_50/nl",
          c(paste("sub111_",paste("eff", "28", sep = ""),sep = ""),
            paste("sub136_",paste("eff", "28", sep = ""),sep = ""),
            paste("sub118_",paste("eff", "28", sep = ""),sep = "")),
          sep = "_")),
  paste("_var",as.character(sigma), sep = ""), 
  sep = "")

repeat_time = 50
## repeat begin

print(paste("start scenario:", scenario_name[n]))
simulation_repeat(sub_ratio = par[[n]]$sub_ratio,
                  eff_size = par[[n]]$eff_size,
                  distribution = par[[n]]$distribution,
                  sigma = sigma,
                  scenario_name = scenario_name[n],
                  repeat_time = repeat_time)


    








