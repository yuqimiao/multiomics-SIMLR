## This is to implement the SNF, SIMLR and GL_SIMLR to the simulated data
source("./code/microb_apply/data_generation.R")
source("./code/functions/SIMLR_multi.R")
source("./code/functions/gl_simlr.R")
source("./code/functions/dist_kernels.R")
library(abSNF)
library(SNFtool)

repeat_time = 30
sample_size = 300
sub_list = list(rep(1/3,3), c(0.2,0.3,0.5))
sub_name = c("111", "235")
compare_tib = NULL
for(t in 1:2){
  sub_ratio = sub_list[[t]] ## balanced case
  for(j in 1:repeat_time){
    ## Sampling from the pool of patients
    sample_ind = NULL
    for(i in 1:3){
      n_sub = sample_size*sub_ratio[i]
      sample_ind[[i]] = sample(which(dataset_abun$Subtype==i),n_sub)
    }
    samples = dataset_abun[unlist(sample_ind),-ncol(dataset_abun)]
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
    gl_res = GL_SIMLR_1(dk_list = DK_list, c = c, B = 30)
    gl_res2 = GL_SIMLR_4(dk_list = DK_list, c = c, B = 30)


    comp = tibble(
      est_cluster = c,
      iter = j,
      sub_ratio = sub_name[t],
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
}

## find clusters that can separate the case subtype into
c_vec = NULL
repeat_time = 20
sample_size = 150
source("./code/microb_apply/data_generation.R")
for(j in 1:repeat_time){
  ## Sampling from the pool of patients
  sample_ind = NULL
  for(i in 1:3){
    n_sub = sample_size*sub_ratio[i]
    sample_ind[[i]] = sample(which(dataset_abun$Subtype==i),n_sub)
  }
  samples = dataset_abun[unlist(sample_ind),-ncol(dataset_abun)]
  true_subtype = dataset_abun[unlist(sample_ind),ncol(dataset_abun)]

  ## calculate the distance and distance kernels
  unifracs = GUniFrac::GUniFrac(otu.tab = samples, throat.tree, alpha = 1)$unifracs
  D_BC = UniFrac_BC(samples) # Bray-Curtis distance
  D_U = unifracs[,,"d_UW"] # Unweighted UniFrac distance
  D_W = unifracs[,,"d_1"]

  D_list = list(D_BC = D_BC, D_U = D_U, D_W = D_W)
  DK_list = lapply(D_list, FUN = dist_kernel)

  ## estimate cluster using eigen gap from mean dk
  DK0 = matrix(0, sample_size,sample_size)
  for (i in 1:length(DK_list)){
    DK0 = DK0+as.matrix(DK_list[[i]])
  }
  DK0 = DK0/length(DK_list)

  c = estimateNumberOfClustersGivenGraph(DK0,NUMC = 2:20)$`Eigen-gap best`
  c_vec = c(c_vec, c)
}




