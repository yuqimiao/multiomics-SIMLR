setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")

# Goal: large scale parameter tunning for gl_simlr_clust:

# scenario setting:
## 2 changes: division_methods(data_divide): div4: 123, div2: 456
##            vars_of_feature(sigma): 1:14, 3:25, 5:36
library(parallel)
library(clValid)
library(dplyr)
library(SNFtool)
library(igraph)
library(Matrix)
library(quadprog)

source("./functions/gl_simlr2.0_cluster.R")
source("./functions/simulation_function.R")
set.seed(286) # set seed because we

scen7 = simulation_3(eff_size = c(1,1,1), data_divide = rep("1/2/3/4",3),sigma = c(1,3,5))
scen8 = simulation_3(eff_size = c(1,1,1), sigma = c(1,3,5))

mk_list = lapply(list(scen7,scen8), function(x){unlist(lapply(x[1:3], multi_kernels_gl))})

names(mk_list) = c('scen7','scen8')

# par setting
par = list()
beta_list = c(0.01,0.03,0.05,0.1,0.3)
gamma_list = c(0.01,0.03,0.05,0.1)
rho_list = c(0.01,0.03,0.05,0.1,0.3,0.5,1)
alpha0_list = alpha_list = c(1,3,5)


for(i in 1:length(beta_list)){
  for(j in 1:length(gamma_list)){
    for(k in 1:length(rho_list)){
      for(l in 1:length(alpha0_list)){
        for(m in 1:length(alpha_list)){
          par_name = paste("beta", beta_list[i],
                           "_gamma", gamma_list[j],
                           "_rho", rho_list[k],
                           "_alpha0", alpha0_list[l],
                           "_alpha", alpha_list[m],
                           sep = "")
          par[[par_name]] = list(beta = beta_list[i],
                                 gamma = gamma_list[j],
                                 rho = rho_list[k],
                                 alpha0 = alpha0_list[l],
                                 alpha = alpha_list[m])
        }
      }
    }
  }
}

save(mk_list, file = "simu_0213/mk_list.Rdata")
save(par, file = "simu_0213/par.Rdata")
