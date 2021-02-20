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


library("grDevices")
library("Matrix")
library("dplyr")
library("stringr")
library("SNFtool")
library("SIMLR")
library("abSNF")

source("./code/functions/SIMLR_multi.R")
source("./code/functions/simulation_function.R")
source("./code/functions/simulation_verify.R")


table_balance = NULL
for(i in 1:5){
  tmp = simulation_3(size = 150,sub_ratio = rep(1/3,3), eff_size = c(1,1,1), data_divide =
                       c("1/23","12/3","2/13"), sigma = 4)
  res = simulation_verify_3data(K =3, tmp,tmp$truelabel)
  table_balance = rbind(table_balance, res)
}

save(table_balance, file = "table_balance.Rdata")









