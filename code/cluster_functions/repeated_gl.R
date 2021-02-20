setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
.libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
task_ID = as.integer(Sys.getenv("SGE_TASK_ID"))
N= 1:18
n = N[task_ID]
library("grDevices")
library("Matrix")
library("dplyr")
library("stringr")
library("SNFtool")
library("SIMLR")
library("abSNF")

source("./functions/simulation_function.R")
source("./functions/gl_simlr_cluster.R")

## simulation scenario
distribution = list(c("normal","normal","normal"),c("normal","normal","logit")) # 2
subtype_ratio = list(rep(0.25,4),c(0.1,0.2,0.3,0.4)) # 2
effective_size = list(c(1,3,6)/10,c(2,3,5),c(5,5,5)) # 3
data_divide = c("12/34","13/24","14/23")


par = list()
m = 1
for(k in 1:length(distribution)){
  for(i in 1:length(subtype_ratio)){
    for(j in 1:length(effective_size)){
      par[[m]] =list(sub_ratio = subtype_ratio[[i]],
                     eff_size = effective_size[[j]],
                     distribution = distribution[[k]]
      )
      m = m+1
    }
  }
}

scenario_name =
  c(paste("nnn",
          c(paste("sub1111_",paste("eff", c("136","244","555"), sep = ""),sep = ""),
            paste("sub1234_",paste("eff", c("136","244","555"), sep = ""),sep = "")),
          sep = "_"),
    paste("nnl",
          c(paste("sub1111_",paste("eff", c("136","244","555"), sep = ""),sep = ""),
            paste("sub1234_",paste("eff", c("136","244","555"), sep = ""),sep = "")),
          sep = "_"))



repeat_time = 50
eva_table = NULL
cluster_res = NULL
for(i in 1:repeat_time){
  ## gl_simlr
  sim = simulation_3(eff_size = par[[n]]$eff_size, sub_ratio = par[[n]]$sub_ratio, dist = par[[n]]$distribution,sigma = rep(100,3) )
  data_list = list(data1 = standardNormalization(t(sim$data1)),
                   data2 = standardNormalization(t(sim$data2)),
                   data3 = standardNormalization(t(sim$data3)))
  gl_result = GL_SIMLR_tunning(data_list = data_list)$best_tunn
  nmi_gl = gl_result$nmi
  rand_gl = gl_result$rand
  # SNF
  aff_mat = NULL
  for(i in 1:length(data_list)){
    aff_mat[[i]] = affinityMatrix(dist2_w(t(data_list[[i]]),t(data_list[[i]]),weight = sim[[i+3]]))
  }
  snf_res = SNF(aff_mat)
  cluster_snf = spectralClustering(snf_res,4)
  cluster_list = list("snf" = cluser_snf, "gl" = gl_result$cluster)
  cluster_res[[i]] = cluster_list
  nmi_snf = igraph::compare(cluster_snf,sim$truelabel,method = "nmi")
  rand_snf = igraph::compare(cluster_snf,sim$truelabel,method = "rand")
  eva = c(scenario = scenario_name[[n]],
          nmi_gl = nmi_gl,
          rand_gl = rand_gl,
          nmi_snf = nmi_snf,
          rand_snf = rand_snf,
          w1 = gl_result$w1,
          w2 = gl_result$w2,
          beta = gl_result$beta,
          gamma = gl_result$gamma,
          rho = gl_result$rho)
  eva_table = rbind(eva_table, eva)
}

save(eva_table, file = paste("./simulation_gl/",scenario_name[[n]],".Rdata", sep = ""))
save(cluster_res, file = paste("./simulation_gl/",scenario_name[[n]],"_cluster",".Rdata", sep = ""))

