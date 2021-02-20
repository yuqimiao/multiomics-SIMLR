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
library("stringr")
library("SNFtool")
library("SIMLR")
library("abSNF")

source("./functions/SIMLR_multi_cluster.R")
source("./functions/simulation_function.R")
source("./functions/simulation_verify_cluster.R")

simulation_repeat = function(sub_ratio, eff_size, distribution, data_divide,scenario_name, repeat_time,sigma = 10,cores.ratio  = 0.5){
  normalize = FALSE
  if(length(unique(distribution))!=1){
    normalize = TRUE
  }
  table = NULL
  for(i in 1: repeat_time){
    sim = simulation_3(size = 240,sub_ratio = sub_ratio, eff_size = eff_size, data_divide = data_divide, dist = distribution,sigma = sigma)
    truelabel = sim$truelabel
    for (K in 2:7){
      tmp1 = NULL
      print(paste("start scenario:", scenario_name, as.character(K)))
      print(paste("repeat time:",as.character(i)))
      t = 0
      while(is.null(tmp1) & t<=3){
        tryCatch({tmp1 = simulation_verify_3data(K,sim,truelabel,normalize = normalize,cores.ratio  = cores.ratio)},
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
distribution = list(c("normal","normal","normal"),c("normal","normal","logit")) # 2 
subtype_ratio = list(rep(1/6,6),c(rep(0.1,2),rep(0.2,4)),c(rep(0.1,4),rep(0.3,2))) # 2
effective_size = list(c(1,3,6),c(2,4,4),c(5,5,5)) # 3
data_divide = c("123/456","145/236","124/356")
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
  c(paste("./simulation_data_var10_s6/nnn",
          c(paste("sub111111_",paste("eff", c("136","244","555"), sep = ""),sep = ""),
            paste("sub112222_",paste("eff", c("136","244","555"), sep = ""),sep = ""),
            paste("sub111133_",paste("eff", c("136","244","555"), sep = ""),sep = "")),
          sep = "_"),
    paste("./simulation_data_var10_s6/nnl",
          c(paste("sub111111_",paste("eff", c("136","244","555"), sep = ""),sep = ""),
            paste("sub112222_",paste("eff", c("136","244","555"), sep = ""),sep = ""),
            paste("sub111133_",paste("eff", c("136","244","555"), sep = ""),sep = "")),
          sep = "_")),
  paste("_var",as.character(sigma), sep = ""), 
  sep = "")

repeat_time = 50
## repeat begin

print(paste("start scenario:", scenario_name[n]))
simulation_repeat(sub_ratio = par[[n]]$sub_ratio,
                  eff_size = par[[n]]$eff_size,
                  distribution = par[[n]]$distribution,
                  data_divide = data_divide,
                  sigma = sigma,
                  scenario_name = scenario_name[n],
                  repeat_time = repeat_time)


    








