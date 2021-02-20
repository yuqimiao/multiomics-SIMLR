## This is to implement the SNF, SIMLR and GL_SIMLR to the simulated data

setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
task_ID = as.integer(Sys.getenv("SGE_TASK_ID"))
N= 1:2
n = N[task_ID]
set.seed(1234)

source("./functions/SIMLR_multi_cluster.R")
source("./functions/gl_simlr_cluster.R")
source("./functions/dist_kernels.R")

library(abSNF)
library(SNFtool)
library(MiSPU)
data("throat.tree")
data("dd")
repeat_time = 50
sample_size = 300
# set_par
sub_list = list(rep(1/3,3), c(0.2,0.3,0.5))
# pool_list = c("micro_case_pool.Rdata","micro_case_pool2.Rdata")

# m = 1
# par_list = NULL
# for(i in 1:length(sub_list)){
#   for(j in 1:length(pool_list)){
#     par_list[[m]] = list(subr = sub_list[[i]],
#                          effpool = pool_list[[j]])
#     m = m+1
#   }
# }

# sub_name = c("111_1","111_2", "235_1", "235_2")
sub_name = c("111", "235")

compare_tib = NULL

# sub_ratio = par_list[[n]]$subr
# load(par_list[[n]]$effpool)
sub_ratio = sub_list[[n]]
load("microb_multisim.Rdata")

throat.tree$tip.label = paste0('OTU', throat.tree$tip.label)
OTUnames = throat.tree$tip.label
UniFrac_BC = function(OTUtab = OTUtab) {return(as.matrix(stats::dist(OTUtab, method = 'manhattan'))/2)}
w_list = list()
cluster_list = list()
for(j in 1:repeat_time){
  ## Sampling from the pool of patients
  sample_ind = NULL
  for(i in 1:3){
    n_sub = sample_size*sub_ratio[i]
    sample_ind[[i]] = sample(which(dataset_abun$Subtype==i),n_sub)
  }
  samples = dataset_abun[unlist(sample_ind),OTUnames]
  true_subtype = dataset_abun[unlist(sample_ind),ncol(dataset_abun)]

  ## calculate the distance and distance kernels
  unifracs = GUniFrac::GUniFrac(otu.tab = samples, throat.tree, alpha = 1)$unifracs
  D_BC = UniFrac_BC(samples) # Bray-Curtis distance
  D_U = unifracs[,,"d_UW"] # Unweighted UniFrac distance
  D_W = unifracs[,,"d_1"] # Weighted UniFrac distance

  D_list = list(D_BC = D_BC, D_U = D_U, D_W = D_W)
  DK_list = lapply(D_list, FUN = dist_kernel)

  ## estimate cluster using eigen gap from mean dk
  DK0 = matrix(0, sample_size,sample_size)
  for (i in 1:length(DK_list)){
    DK0 = DK0+as.matrix(DK_list[[i]])
  }
  DK0 = DK0/length(DK_list)

  c = estimateNumberOfClustersGivenGraph(DK0,NUMC = 2:20)$`Eigen-gap best`

  ## compare the performance

  ### SNF
  aff_list = NULL
  for(i in 1:length(D_list)){
    aff_list[[i]] = affinityMatrix(D_list[[i]],K = 30)
  }
  W = SNF(aff_list,K = 30)
  cluster_snf = spectralClustering(W,K = c)

  ## implementing
  simlr_res = SIMLR_multi(DK_list, c = c)

  ## implementing gl_simlr with both default hyper and estimated hyper
  gl_res = GL_SIMLR_1(dk_list = DK_list, c = c, B = 30,gamma = 5)
  gl_res2 = GL_SIMLR_4(D_list = D_list, c = c, B = 30,normal_type = 1)

  ## save weight and cluster result
  w_mtd = list(gl1 = gl_res$w_list, gl2 = gl_res2$w_list, simlr = simlr_res$alphaK)
  w_list[[j]] = w_mtd
  clustr_mtd = list(gl1 = gl_res$cluster, gl2 = gl_res2$cluster, simlr = simlr_res$y$cluster, snf = cluster_snf)
  cluster_list[[j]] = clustr_mtd
  ## compare methods tibble
  comp = tibble(
    est_cluster = c,
    iter = j,
    sub_ratio = sub_name[n],
    rand_snf = igraph::compare(cluster_snf, true_subtype, method = "rand"),
    rand_simlr = igraph::compare(simlr_res$y$cluster, true_subtype, method = "rand"),
    rand_gldef = igraph::compare(gl_res$cluster, true_subtype, method = "rand"),
    rand_glest = igraph::compare(gl_res2$cluster, true_subtype, method = "rand"),
    nmi_snf = igraph::compare(cluster_snf, true_subtype, method = "nmi"),
    nmi_simlr = igraph::compare(simlr_res$y$cluster, true_subtype, method = "nmi"),
    nmi_gldef = igraph::compare(gl_res$cluster, true_subtype, method = "nmi"),
    nmi_glest = igraph::compare(gl_res2$cluster, true_subtype, method = "nmi")
  )
  compare_tib = rbind(compare_tib,comp)
}

save(compare_tib, file = paste("simulation_microgl/",sub_name[n],".Rdata",sep = ""))
save(w_list, file = paste("simulation_microgl/",sub_name[n],"w_list.Rdata",sep = ""))
save(cluster_list, file = paste("simulation_microgl/",sub_name[n],"cluster_list.Rdata",sep = ""))

