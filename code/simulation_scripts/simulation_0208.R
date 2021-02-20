setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
task_ID = as.integer(Sys.getenv("SGE_TASK_ID"))
N= 1:300
n = N[task_ID]

# Goal: large scale parameter tunning for gl_simlr_clust
library(parallel)
library(clValid)
library(dplyr)
library(purrr)
library(SNFtool)
library(igraph)
library(Matrix)
library(quadprog)

source("functions/CIMLR.R")
source("./R/compute.multiple.kernel.cimlr.R")
source("./functions/gl_simlr2.0_cluster.R")
source("./functions/simulation_function.R")
# simulation
scen1 = simulation_3(eff_size = c(1,1,1), data_divide = rep("1/2/3/4",3),sigma = rep(1,3))
scen2 = simulation_3(eff_size = c(1,1,1), data_divide = rep("1/2/3/4",3),sigma = rep(3,3))
scen3 = simulation_3(eff_size = c(1,1,1), data_divide = rep("1/2/3/4",3),sigma = rep(5,3))
scen4 = simulation_3(eff_size = c(1,1,1), sigma = rep(1,3))
scen5 = simulation_3(eff_size = c(1,1,1), sigma = rep(3,3))
scen6 = simulation_3(eff_size = c(1,1,1), sigma = rep(5,3))

data_list = list(scen1,scen2,scen3,scen4,scen5,scen6)
names(data_list) = c('scen1','scen2','scen3','scen4','scen5','scen6')

# for gl_simlr

## multiple_kernel
mk_list = lapply(data_list, function(x){unlist(lapply(x[1:3], multi_kernels_gl))})

## best_par tunned from simu_seed286
load("simu_0207/best_tuned_par.Rdata")

## implement
res_gl = lapply(names(mk_list), FUN = function(x){
  par1 = best_par[[x]]
  nmi_f = 0
  nmi_ls_tmp = NULL
  m = 1
  while(nmi_f!=1 & m <= 30) {
    tryCatch({ res = Crazy_5term_clust(kernel_list = mk_list[[x]],
                                       c=4,
                                       alpha0 = par1[m,]$alpha0,
                                       alpha = par1[m,]$alpha,
                                       rho = par1[m,]$rho,
                                       beta = par1[m,]$beta,
                                       gamma = par1[m,]$gamma,
                                       normalize = 0)
    print(length(res))
    cluster = spectralClustering(as.matrix(res$S), K = 4)
    nmi_f = igraph::compare(rep(1:4, each = 50), cluster, method = "nmi")
    nmi_ls_tmp = c(nmi_ls_tmp, nmi_f)
    if(nmi_f == max(nmi_ls_tmp)){
      res_max = res
      res_max[["par"]] = par1[m,]
    }
    if(dim(par1)[1] == 1){ break }
    m = m+1},
    error = function(e) {cat(conditionMessage(e))})
  }
  return(res_max)
})

save(res_gl, file = paste("simu_0208/res_gl",n, ".Rdata", sep = ""))

## evaluate
nmi_gl_ls = map_dbl(res_gl, function(x){
  cluster = spectralClustering(as.matrix(x$S), K = 4)
  nmi =igraph::compare(rep(1:4, each = 50), cluster, method = "nmi")
  nmi
})
unmatched_mk = names(mk_list)[-match(names(nmi_gl_ls),names(mk_list))]
for(i in unmatched_mk){
  nmi_gl_ls[[i]] = NA
}


# for CIMLR
res_cim = lapply(data_list, FUN = function(x){
  tryCatch({
    CIMLR(lapply(x[1:3], FUN = function(x) t(x)),4,cores.ratio = 0.3)
  },
  error = function(e) {cat(conditionMessage(e))})
})
save(res_cim, file = paste("simu_0208/res_cim",n, ".Rdata", sep = ""))

## evaluate
nmi_cim_ls = map_dbl(res_cim, function(x){
  nmi =igraph::compare(rep(1:4, each = 50), x$y$cluster, method = "nmi")
  nmi
})


# final evaluation stats
gl_par_table = map_dfr(res_gl, .f = function(x) x$par)%>% select(-ID, -scenario, -nmi)

eva_table = as_tibble(cbind(nmi_gl_ls,nmi_cim_ls)) %>%
  mutate(scenario = c('scen1','scen2','scen3','scen4','scen5','scen6'),
         w1_gl = map_dbl(res_gl, .f = function(x) {sum(x$w[1:22])}),
         w2_gl = map_dbl(res_gl, .f = function(x) {sum(x$w[23:44])}),
         w1_cim = map_dbl(res_cim, .f = function(x) {sum(x$alphaK[1:55])}),
         w2_cim = map_dbl(res_cim, .f = function(x) {sum(x$alphaK[56:110])}))

eva_table = cbind(eva_table, gl_par_table)
save(eva_table, file = paste("simu_0208/eva_table", n,".Rdata", sep = ""))

