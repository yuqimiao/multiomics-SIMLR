setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
task_ID = as.integer(Sys.getenv("SGE_TASK_ID"))
N= 1:3087
n = N[task_ID]

# Goal: large scale parameter tunning for gl_simlr_clust
library(parallel)
library(clValid)
library(dplyr)
library(SNFtool)
library(igraph)
library(Matrix)
library(quadprog)

load(file = "simu_0206/mk_list.Rdata")
load(file = "simu_0206/par.Rdata")
source("./functions/gl_simlr2.0_cluster.R")
# source("code/functions/simulation_function.R")
# source("code/functions/gl_simlr2.0.R")
# source("code/functions/gl_simlr.R")


par1 = par[[n]]
res = lapply(mk_list, FUN = function(x){
  tryCatch({Crazy_5term_clust(kernel_list = x,
                              c=4,
                              alpha0 = par1$alpha0,
                              alpha = par1$alpha,
                              rho = par1$rho,
                              beta = par1$beta,
                              gamma = par1$gamma,
                              normalize = 0)},
           error = function(e) {cat(conditionMessage(e))})
})

nmi_ls = lapply(res, FUN = function(x){
  cluster = spectralClustering(as.matrix(x$S), K = 4)
  nmi =igraph::compare(rep(1:4, each = 50), cluster, method = "nmi")
  nmi
})

unmatched_mk = names(mk_list)[-match(names(nmi_ls),names(mk_list))]
for(i in unmatched_mk){
  nmi_ls[[i]] = NA
}

save(res, file = paste("simu_0206/res_",names(par)[n], ".Rdata", sep = ""))
save(nmi_ls, file = paste("simu_0206/nmi_",names(par)[n], ".Rdata", sep = ""))

