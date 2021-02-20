setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
task_ID = as.integer(Sys.getenv("SGE_TASK_ID"))
N= 1:100
n = N[task_ID]

library(SNFtool)
library(abSNF)
library(dplyr)
library(stringr)
library(igraph)
library(abSNF)
source("./functions/simulation_function.R")
source("./functions/gl_simlr_cluster.R")
# source("code/functions/simulation_function.R")
# source("code/functions/gl_simlr.R")

par = list()
scenario_name = list()
data_divide_list = list(c("12/34","13/24","14/23"), rep("1/2/3/4",3))
sub_ratio_ls = list(rep(1/4,4), c(0.1,0.2,0.3,0.4))
eff_size_ls = list(c(1,1,0),c(1,3,3),c(1,3,5))
sigma_ls = c(1,10,100)

name_list = list(
  ratio = c("sub1111", "sub1234"),
  eff = c("eff110","eff133","eff135"),
  sigma = paste("sig", sigma_ls, sep = ""),
  divide = c("div2", "div4")
)

m = 1
for(i in 1:length(sub_ratio_ls)){
  for(j in 1:length(eff_size_ls)){
    for(k in 1:length(sigma_ls)){
      for(l in 1:length(data_divide_list)){
        par[[m]] = list(sub_ratio = sub_ratio_ls[[i]],
                        eff_size = eff_size_ls[[j]],
                        sigma = sigma_ls[k],
                        data_divide = data_divide_list[[l]])
        scenario_name[[m]] = paste(name_list$ratio[i],
                                   name_list$eff[j],
                                   name_list$sigma[k],
                                   name_list$divide[l], sep = "_")
        m = m+1
      }
    }
  }
}

scenario_name = unlist(scenario_name)

## data generation
eva_nmi = NULL
eva_rand = NULL
sim_list_all = list()
gl_weight = NULL
global_local = NULL
GL_trace = NULL
K = 30
tunning_res = NULL
tunning_best = NULL
for(m in 1:length(par)){
  data = simulation_3(sub_ratio = par[[m]]$sub_ratio,
                      sigma = rep(par[[m]]$sigma,3),
                      eff_size = par[[m]]$eff_size,
                      data_divide = par[[m]]$data_divide)
  truth = data$truelabel
  #build similarity matrix

  ## Single data
  data1 = data$data1
  dist1 = dist2(data1,data1)^(1/2)
  W1 = affinityMatrix(dist1, K)

  data2 = data$data2
  dist2 = dist2(data2,data2)^(1/2)
  W2 = affinityMatrix(dist2, K)

  data3 = data$data3
  dist3 = dist2(data3,data3)^(1/2)
  W3 = affinityMatrix(dist3, K)

  ## concatenate
  data_c = cbind(data$data1, data$data2,data$data3)
  dist_c = dist2(data_c,data_c)^(1/2)
  W_c = affinityMatrix(dist_c, K)

  ## SNF
  dist21 = dist2(data$data1,data$data1)^(1/2)
  dist22 = dist2(data$data2,data$data2)^(1/2)
  dist23 = dist2(data$data3,data$data3)^(1/2)
  W21 = affinityMatrix(dist21,K)
  W22 = affinityMatrix(dist22,K)
  W23 = affinityMatrix(dist23,K)
  W_s = SNF(list(W21,W22,W23), K)

  ## abSNF
  dist21 = dist2_w(data$data1,data$data1,weight = data$weight1)^(1/2)
  dist22 = dist2_w(data$data2,data$data2,weight = data$weight2)^(1/2)
  dist23 = dist2_w(data$data3,data$data3,weight = data$weight3)^(1/2)
  W21 = affinityMatrix(dist21,K)
  W22 = affinityMatrix(dist22,K)
  W23 = affinityMatrix(dist23,K)
  W_abs = SNF(list(W21,W22,W23), K)

  ##gl-simlr
  W_g = GL_SIMLR_4(data_list = list(t(data$data1),t(data$data2),t(data$data3)), c = 4,normal_type = 0)
  # W_g_nofnorm = GL_SIMLR_4(data_list = list(t(data$data1),t(data$data2),t(data$data3)), c = 4,F_norm = F)
  W_g_nbeta = GL_SIMLR_4(data_list = list(t(data$data1),t(data$data2),t(data$data3)), c = 4, beta = 0)
  W_g_tunn = GL_SIMLR_tunning(data_list = list(t(data$data1),t(data$data2),t(data$data3)),
                              c = 4,
                              # gamma_grid = seq(20,30,by = 5),
                              # beta_grid = seq(20,30,by = 5),
                              # rho_grid = seq(0.01,0.05,by = 0.02),
                              rho_grid = 0.1,
                              gamma_grid = 10,
                              beta_grid = 10,
                              true_label = data$truelabel)
  # clustering result

  ## single data
  cluster_kmeans1 = kmeans(data1,4)$cluster
  cluster_hclust1 = cutree(hclust(as.dist(dist1),method = "complete"),4)
  cluster_spec1 = spectralClustering(W1, 4)

  cluster_kmeans2 = kmeans(data2,4)$cluster
  cluster_hclust2 = cutree(hclust(as.dist(dist2),method = "complete"),4)
  cluster_spec2 = spectralClustering(W2, 4)

  cluster_kmeans3 = kmeans(data3,4)$cluster
  cluster_hclust3 = cutree(hclust(as.dist(dist3),method = "complete"),4)
  cluster_spec3 = spectralClustering(W3, 4)

  ## concatenate data
  cluster_kmeans_c = kmeans(data_c,4)$cluster
  cluster_hclust_c = cutree(hclust(as.dist(dist_c),method = "complete"),4)
  cluster_spec_c = spectralClustering(W_c, 4)

  ## SNF
  cluster_spec_s = spectralClustering(W_s,4)
  ## SNF
  cluster_spec_abs = spectralClustering(W_abs,4)
  ## gl_simlr
  cluster_spec_g = W_g$cluster
  # cluster_spec_g_nofnorm = W_g_nofnorm$cluster
  cluster_spec_g_nobeta = W_g_nbeta$cluster
  cluster_spec_g_tunn = W_g_tunn$res_opt$cluster
  #evaluation
  eva_table_nmi = tibble(scenario = scenario_name[m],
                         nmi_kmean1 = compare(truth, cluster_kmeans1, method = "nmi"),
                         nmi_kmean2 = compare(truth, cluster_kmeans2, method = "nmi"),
                         nmi_kmean3 = compare(truth, cluster_kmeans3, method = "nmi"),
                         nmi_kmean_c = compare(truth, cluster_kmeans_c, method = "nmi"),
                         nmi_hclust1 = compare(truth, cluster_hclust1, method = "nmi"),
                         nmi_hclust2 = compare(truth, cluster_hclust2, method = "nmi"),
                         nmi_hclust3 = compare(truth, cluster_hclust3, method = "nmi"),
                         nmi_hclust_c = compare(truth, cluster_hclust_c, method = "nmi"),
                         nmi_spec1 = compare(truth, cluster_spec1, method = "nmi"),
                         nmi_spec2 = compare(truth, cluster_spec2, method = "nmi"),
                         nmi_spec3 = compare(truth, cluster_spec3, method = "nmi"),
                         nmi_spec_c = compare(truth, cluster_spec_c, method = "nmi"),
                         nmi_spec_s = compare(truth, cluster_spec_s, method = "nmi"),
                         nmi_spec_abs = compare(truth, cluster_spec_abs, method = "nmi"),
                         nmi_spec_g = compare(truth, cluster_spec_g, method = "nmi"),
                         # nmi_spec_g_nofnorm = compare(truth, cluster_spec_g_nofnorm, method = "nmi"),
                         nmi_spec_g_nobeta = compare(truth, cluster_spec_g_nobeta, method = "nmi"),
                         nmi_spec_gt = compare(truth, cluster_spec_g_tunn, method = "nmi")
  )
  eva_table_rand = tibble(scenario = scenario_name[m],
                          rand_kmean1 = compare(truth, cluster_kmeans1, method = "rand"),
                          rand_kmean2 = compare(truth, cluster_kmeans2, method = "rand"),
                          rand_kmean3 = compare(truth, cluster_kmeans3, method = "rand"),
                          rand_kmean_c = compare(truth, cluster_kmeans_c, method = "rand"),
                          rand_hclust1 = compare(truth, cluster_hclust1, method = "rand"),
                          rand_hclust2 = compare(truth, cluster_hclust2, method = "rand"),
                          rand_hclust3 = compare(truth, cluster_hclust3, method = "rand"),
                          rand_hclust_c = compare(truth, cluster_hclust_c, method = "rand"),
                          rand_spec1 = compare(truth, cluster_spec1, method = "rand"),
                          rand_spec2 = compare(truth, cluster_spec2, method = "rand"),
                          rand_spec3 = compare(truth, cluster_spec3, method = "rand"),
                          rand_spec_c = compare(truth, cluster_spec_c, method = "rand"),
                          rand_spec_s = compare(truth, cluster_spec_s, method = "rand"),
                          rand_spec_abs = compare(truth, cluster_spec_abs, method = "rand"),
                          rand_spec_g = compare(truth, cluster_spec_g, method = "rand"),
                          # rand_spec_g_nofnorm = compare(truth, cluster_spec_g_nofnorm, method = "rand"),
                          rand_spec_g_nobeta = compare(truth, cluster_spec_g_nobeta, method = "rand"),
                          rand_spec_gt = compare(truth, cluster_spec_g_tunn, method = "rand")
  )
  eva_nmi = rbind(eva_nmi, eva_table_nmi)
  eva_rand = rbind(eva_rand, eva_table_rand)
  ## save similarity matrix
  # sim_list = list(W1 = W1, W2 = W2, W_c = W_c, W_s = W_s)
  # sim_list_all[[m]] = sim_list
  gl_weight = rbind(gl_weight,cbind(scenario = rep(scenario_name[m], 3),
                                    type = 1:3,
                                    GL = W_g$w_list,
                                    # GL_nfnorm = W_g_nofnorm$w_list,
                                    GL_nbeta = W_g_nbeta$w_list,
                                    GL_tunn = W_g_tunn$res_opt$w_list))
  global_local = rbind(global_local, cbind(scenario = rep(scenario_name[m], 3),
                                           type = 1:3,
                                           W_g$global_local_term,
                                           # W_g_nofnorm$global_local_term,
                                           W_g_nbeta$global_local_term,
                                           W_g_tunn$res_opt$global_local_term)
  )
  GL_trace = rbind(GL_trace, cbind(scenario = rep(scenario_name[m], 3),
                                   type = 1:3,
                                   GL = W_g$GL_trace,
                                   # GL_nfnorm = W_g_nofnorm$GL_trace,
                                   GL_nbeta = W_g_nbeta$GL_trace,
                                   GL_tunn = W_g_tunn$res_opt$GL_trace))
  tunning_res = rbind(tunning_res, as_tibble(W_g_tunn$tunning_crit) %>% mutate(scenario = scenario_name[i]))
  tunning_best = rbind(tunning_best, as_tibble(list(rho_opt = W_g_tunn$rho_opt,
                                                    gamma_opt = W_g_tunn$gamma_opt,
                                                    beta_opt = W_g_tunn$beta_opt,
                                                    gl_opt = W_g_tunn$gl_opt,
                                                    GL_trace_opt = W_g_tunn$GL_trace_opt,
                                                    scenario = scenario_name[m]))
  )
}
names(gl_weight) = c("scenario", "type","GL","GL_nbeta","GL_tunn")
names(GL_trace) = c("scenario", "type","GL","GL_nbeta","GL_tunn")
names(global_local) = c("scenario", "type","global_GL","global_GL_nbeta","global_GL_tunn","local_GL","local_GL_nbeta","local_GL_tunn")
# res = list(eva_nmi = eva_nmi,eva_rand = eva_rand, sim_list_all = sim_list_all)
res = list(eva_nmi = eva_nmi,
           eva_rand = eva_rand,
           GL_trace = GL_trace,
           gl_weight = gl_weight,
           global_local = global_local,
           tunning_res = tunning_res,
           tunning_best = tunning_best)
#
save(res,
     file = paste("simu_1224/sim_res",n,".Rdata", sep = ""))







