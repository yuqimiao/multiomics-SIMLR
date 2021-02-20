# setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
# .libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
# task_ID = as.integer(Sys.getenv("SGE_TASK_ID"))
# N= 1:100
# n = N[task_ID]
# source("./functions/simulation_function.R")
# source("./functions/gl_simlr.R")


library(SNFtool)
library(abSNF)
library(dplyr)
library(stringr)
library(igraph)
source("code/functions/simulation_function.R")
source("code/functions/gl_simlr.R")

par = list()
sub_ratio_ls = list(rep(1/3,3), c(0.2,0.3,0.5))
eff_size_ls = list(c(1,1),c(1,3),c(1,5))
sigma_ls = c(1,10,100)
m = 1
for(i in 1:length(sub_ratio_ls)){
  for(j in 1:length(eff_size_ls)){
    for(k in 1:length(sigma_ls)){
      par[[m]] = list(sub_ratio = sub_ratio_ls[[i]],
                      eff_size = eff_size_ls[[j]],
                      sigma = sigma_ls[k])
      m = m+1
    }
  }
}
scenario_name = c(paste("sub111", c(paste("eff11", paste("sig",sigma_ls, sep = ""),sep = "_"),
                                  paste("eff13", paste("sig",sigma_ls, sep = ""),sep = "_"),
                                  paste("eff15", paste("sig",sigma_ls, sep = ""),sep = "_")),
                      sep = "_"),
                  paste("sub235", c(paste("eff11", paste("sig",sigma_ls, sep = ""),sep = "_"),
                                    paste("eff13", paste("sig",sigma_ls, sep = ""),sep = "_"),
                                    paste("eff15", paste("sig",sigma_ls, sep = ""),sep = "_")),
                        sep = "_")
                  )

## data generation
eva_nmi = NULL
eva_rand = NULL
gl_weight = NULL
sim_list_all = list()
K = 30
for(m in 1:length(par)){
  data = simulation_2(sub_ratio = par[[m]]$sub_ratio,
                      sigma = rep(par[[m]]$sigma,2),
                      eff_size = par[[m]]$eff_size)
  truth = data$truelabel
  #build similarity matrix

  ## Single data
  data1 = data$data1
  dist1 = dist2(data1,data1)^(1/2)
  W1 = affinityMatrix(dist1, K)

  data2 = data$data2
  dist2 = dist2(data2,data2)^(1/2)
  W2 = affinityMatrix(dist2, K)

  ## concatenate
  data_c = cbind(data$data1, data$data2)
  dist_c = dist2(data_c,data_c)^(1/2)
  W_c = affinityMatrix(dist_c, K)

  ## SNF
  dist21 = dist2(data$data1,data$data1)^(1/2)
  dist22 = dist2(data$data2,data$data2)^(1/2)
  W21 = affinityMatrix(dist21,K)
  W22 = affinityMatrix(dist22,K)
  W_s = SNF(list(W21,W22), K)

  ##gl-simlr
  W_g = GL_SIMLR_4(data_list = list(t(data$data1),t(data$data2)), c = 3)
  # clustering result

  ## single data
  cluster_kmeans1 = kmeans(data1,3)$cluster
  cluster_hclust1 = cutree(hclust(as.dist(dist1),method = "complete"),3)
  cluster_spec1 = spectralClustering(W1, 3)

  cluster_kmeans2 = kmeans(data2,3)$cluster
  cluster_hclust2 = cutree(hclust(as.dist(dist2),method = "complete"),3)
  cluster_spec2 = spectralClustering(W2, 3)

  ## concatenate data
  cluster_kmeans_c = kmeans(data_c,3)$cluster
  cluster_hclust_c = cutree(hclust(as.dist(dist_c),method = "complete"),3)
  cluster_spec_c = spectralClustering(W_c, 3)

  ## SNF
  cluster_spec_s = spectralClustering(W_s,3)
  ## gl_simlr
  cluster_spec_g = W_g$cluster
  #evaluation
  eva_table_nmi = tibble(scenario = scenario_name[m],
                         nmi_kmean1 = compare(truth, cluster_kmeans1, method = "nmi"),
                         nmi_kmean2 = compare(truth, cluster_kmeans2, method = "nmi"),
                         nmi_kmean_c = compare(truth, cluster_kmeans_c, method = "nmi"),
                         nmi_hclust1 = compare(truth, cluster_hclust1, method = "nmi"),
                         nmi_hclust2 = compare(truth, cluster_hclust2, method = "nmi"),
                         nmi_hclust_c = compare(truth, cluster_hclust_c, method = "nmi"),
                         nmi_spec1 = compare(truth, cluster_spec1, method = "nmi"),
                         nmi_spec2 = compare(truth, cluster_spec2, method = "nmi"),
                         nmi_spec_c = compare(truth, cluster_spec_c, method = "nmi"),
                         nmi_spec_s = compare(truth, cluster_spec_s, method = "nmi"),
                         nmi_spec_g = compare(truth, cluster_spec_g, method = "nmi")
                        )
  eva_table_rand = tibble(scenario = scenario_name[m],
                          rand_kmean1 = compare(truth, cluster_kmeans1, method = "rand"),
                          rand_kmean2 = compare(truth, cluster_kmeans2, method = "rand"),
                          rand_kmean_c = compare(truth, cluster_kmeans_c, method = "rand"),
                          rand_hclust1 = compare(truth, cluster_hclust1, method = "rand"),
                          rand_hclust2 = compare(truth, cluster_hclust2, method = "rand"),
                          rand_hclust_c = compare(truth, cluster_hclust_c, method = "rand"),
                          rand_spec1 = compare(truth, cluster_spec1, method = "rand"),
                          rand_spec2 = compare(truth, cluster_spec2, method = "rand"),
                          rand_spec_c = compare(truth, cluster_spec_c, method = "rand"),
                          rand_spec_s = compare(truth, cluster_spec_s, method = "rand"),
                          rand_spec_g = compare(truth, cluster_spec_g, method = "rand")
  )
  eva_nmi = rbind(eva_nmi, eva_table_nmi)
  eva_rand = rbind(eva_rand, eva_table_rand)
  ## save similarity matrix
  # sim_list = list(W1 = W1, W2 = W2, W_c = W_c, W_s = W_s)
  # sim_list_all[[m]] = sim_list
  gl_weight = rbind(gl_weight,W_g$w_list)
}
# res = list(eva_nmi = eva_nmi,eva_rand = eva_rand, sim_list_all = sim_list_all)
res = list(eva_nmi = eva_nmi,eva_rand = eva_rand, gl_weight = gl_weight)

# save(res,
#      file = paste("simu_snf/sim_res",n,".Rdata", sep = ""))



as_tibble(res$gl_weight) %>%
  mutate(scenario = str_replace_all(scenario_name,"[a-z]","")) %>%
  pivot_longer(
    V1:V2,
    names_to = "data",
    values_to = "weight"
  ) %>%
  ggplot(aes(x = scenario, y = weight, color = data, group = data))+
  geom_line()



