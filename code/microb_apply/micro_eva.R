library('GUniFrac')
library(igraph)
UniFrac_BC = function(OTUtab = OTUtab) {return(as.matrix(stats::dist(OTUtab, method = 'manhattan'))/2)}

micro_cluster_eva = function(sample_size = 300,
                             sub_ratio = rep(1/3,3),
                             sub_name = "111",
                             dataset_abun,
                             tree = throat.tree,
                             weight_boosting = F,
                             pv = NA,
                             threshold = 0.5){
  n_leaf = NA
  sample_ind = NULL
  for(i in 1:length(sub_ratio)){
    n_sub = sample_size*sub_ratio[i]
    sample_ind[[i]] = sample(which(dataset_abun$Subtype==i),n_sub)
  }
  samples = dataset_abun[unlist(sample_ind),-ncol(dataset_abun)]
  # samples = apply(samples, MARGIN = 2, FUN = function(x)(x-mean(x))/sd(x))
  # samples = samples*-log10(pv)/rowSums(samples*-log10(pv))
  true_subtype = dataset_abun[unlist(sample_ind),ncol(dataset_abun)]
  if(class(true_subtype) != "numeric"){true_subtype = true_subtype[[1]]}

  ## calculate the distance and distance kernels
  if(weight_boosting){
    tree$edge.length = tree$edge.length * -log10(ifelse(pv>threshold,1,pv))
    n_leaf = sum(pv>threshold)
  }
  unifracs = GUniFrac::GUniFrac(otu.tab = samples, tree)$unifracs
  D_BC = UniFrac_BC(samples) # Bray-Curtis distance
  D_U = unifracs[,,"d_UW"] # Unweighted UniFrac distance
  D_W = unifracs[,,"d_1"] # Weighted UniFrac distance

  D_list = list(D_BC = D_BC, D_U = D_U, D_W = D_W)
  DK_list = lapply(D_list, FUN = function(x) {dist_kernel(x, is.dist = T)})

  mds_list = lapply(D_list, function(x){
    mds = cmdscale(dist,2)
    return(mds)
  })
  g_sim1 = as_tibble(list.rbind(mds_list)) %>%
    mutate(dist_type = rep(names(D_list), each = 200),
           true = rep(true_subtype, 3)) %>%
    ggplot(aes(x = V1, y = V2, color = factor(true)))+
    geom_point()+
    facet_grid(data_type~., scale = "free")

  g_sim1
  ggsave(filename = "micro_samples_mds.jpeg")

  ## estimate cluster using eigen gap from mean dk
  DK0 = matrix(0, sample_size,sample_size)
  for (i in 1:length(DK_list)){
    DK0 = DK0+as.matrix(DK_list[[i]])
  }
  DK0 = DK0/length(DK_list)
  S0 = max(max(DK0))-DK0

  c = estimateNumberOfClustersGivenGraph(S0,NUMC = 2:20)$`Eigen-gap best`

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
  gl_res = gl_simlr(dist_list = D_list, c = c, B = 30)
  # gl_res2 = GL_SIMLR_4(D_list = D_list, c = c, B = 30,normal_type = 1,rho = 0.05)
  # gl_tune_res = GL_SIMLR_tunning(D_list = D_list, c = c,gamma_grid = NA,rho_grid = seq(0.001,0.1, length.out = 10))

  comp = tibble(
    est_cluster = c,
    iter = 1,
    sub_ratio = sub_name,
    n_leaf = n_leaf,
    rand_snf = igraph::compare(cluster_snf, true_subtype, method = "rand"),
    rand_simlr = igraph::compare(simlr_res$y$cluster, true_subtype, method = "rand"),
    rand_gldef = igraph::compare(gl_res$cluster, true_subtype, method = "rand"),
    # rand_glest = igraph::compare(gl_res2$cluster, true_subtype, method = "rand"),
    nmi_snf = igraph::compare(cluster_snf, true_subtype, method = "nmi"),
    nmi_simlr = igraph::compare(simlr_res$y$cluster, true_subtype, method = "nmi"),
    nmi_gldef = igraph::compare(gl_res$cluster, true_subtype, method = "nmi"),
    # nmi_glest = igraph::compare(gl_res2$cluster, true_subtype, method = "nmi")
  )


  # data = tibble(true = true_subtype, snf = cluster_snf, simlr = simlr_res$y$cluster, gl1 = gl_res$cluster, gl2 = gl_res2$cluster)
  data = tibble(true = true_subtype, snf = cluster_snf, simlr = simlr_res$y$cluster, gl1 = gl_res$cluster)
  # w_mtd = list(gl1 = gl_res$w_list, gl2 = gl_res2$w_list, simlr = simlr_res$alphaK)
  w_mtd = list(gl1 = gl_res$w, simlr = simlr_res$alphaK)

  simlarity_matrix = list(snf_sim = W, gl_sim = gl_res$k,simlr_sim = simlr_res$S)

  return(list(comp = comp, data = data,w_mtd = w_mtd,simlarity_matrix = simlarity_matrix))
}


