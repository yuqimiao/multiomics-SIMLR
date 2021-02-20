## package loading
library(tidyverse)
library(parallel)
library(clValid)
library(abSNF)
library(SNFtool)
library(igraph)
library(Matrix)
library(CancerSubtypes)
library(rlist)

source("code/functions/simulation_function.R")
source("code/R/compute.multiple.kernel.R")
source("code/functions/simulation_verify.R")
source("./code/functions/gl_simlr2.0.R")
source("./code/functions/SIMLR_multi.R")
source("./code/functions/SIMLR_no_weights.R")
library(RNOmni)

## weight derivation

all_ind_weight = list()
for(i in 1:100){
  data = simulation_3(eff_size = c(5,5,0),
                      sigma =c(10,10,100),
                      data_divide = c("12/34","13/24","14/23"))
  data_list = data[1:3]
  ## construct kernels
  dist_list = NA
  dist_kernel_list = lapply(data_list, FUN = function(x) D_kernels(x_fun = x) )
  kernel_list = lapply(dist_kernel_list, FUN = function(x){x[[1]]})
  ## calculate square of rank normalization,
  rank2_list = lapply(kernel_list, FUN = function(x) {
    t(apply(x,MARGIN = 1, function(y){2^(RNOmni::rankNorm(y))}))
  })
  ## normalization to get weights
  ind_weight = lapply(rank2_list, function(x) x/Reduce("+",rank2_list))
  all_ind_weight[[i]] = ind_weight
}
save(all_ind_weight, file = "./data/ind_weight_simulation/ind_weight_550_10_div2.Rdata")
