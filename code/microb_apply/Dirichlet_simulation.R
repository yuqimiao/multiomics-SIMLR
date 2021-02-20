library(DAFOT)
# BiocManager::install("metagenomeFeatures")
library(metagenomeFeatures)
gg85 <- get_gg13.8_85MgDb()
gamma_16S <- mgDb_select(gg85, type = "all", keys = "Gammaproteobacteria", keytype = "Class")
Tree=gamma_16S$tree
Tree$tip.label<-1:length(Tree$tip.label)
alphaP=c(rep(1,length(Tree$tip.label)))## only simulate leaves, no internal nodes
alphaQ=c(rep(1,length(Tree$tip.label)))## only simulate leaves, no internal nodes, remove rep(0,Tree$Nnode)
alphaQ[1:10]=alphaQ[1:10]+2
alphaQ[21:30]=alphaQ[21:30]-0.5

DataPQ<-DataGenerating(200,200,alphaP,alphaQ,10000)
# DAFOT(DataPQ$P,DataPQ$Q,Tree,100,0.05)

## Format transformation
sub1 = t(DataPQ$P)
sub2 = t(DataPQ$Q)
pool = rbind(sub1,sub2)
colnames(pool) = Tree$tip.label
pool = as_tibble(pool) %>% mutate(Subtype = rep(c(1,2), each = 200))


mds_new = sample_mds(pool,sub_ratio = rep(1/2,2), sample_size = 200, tree = Tree)
plot_generate(mds_new,names = "new_micro")

eva_new = micro_cluster_eva(sample_size = 200, sub_ratio = rep(1/2,2),dataset_abun = pool, tree = Tree)

t(eva_new$comp)
eva_new$w_mtd
eva_new$gl_tune_res$best_tunn
eva_new$gl_tune_res$tunning_all
## cluster ==> phylo-related/unrelated/pres-abs/abun

simlr_list = lapply(eva_new$simlarity_matrix, FUN = function(x)(apply(as.matrix(x),MARGIN = 1,FUN = function(x){(x-min(x))/(max(x)-min(x))})))
simlr_list_d0 = lapply(eva_new$simlarity_matrix, function(x){
  x = as.matrix(x)
  diag(x) = 0
  return(x)})

ydata_snf = tsne(simlr_list_d0[[1]])
ydata_gl = tsne(simlr_list_d0[[2]])
ydata_simlr = tsne(simlr_list_d0[[3]])
# ydata_gl = tsne(as.matrix(micro_mn$simlarity_matrix$gl2_sim))
# ydata_snf = tsne(as.matrix(micro_mn$simlarity_matrix$snf_sim))
# ydata_simlr = tsne(as.matrix(micro_mn$simlarity_matrix$simlr_sim))

ydata = as_tibble(rbind(ydata_gl,ydata_snf,ydata_simlr)) %>%
  mutate(method = rep(c("gl","snf","simlr"), each = 200),
         true_label = rep(rep(c(1,2), each = 100), 3))
g1 = ydata %>%
  filter(method == "gl") %>%
  ggplot(aes(V1,V2, color = factor(true_label))) +
  geom_point(size = 1, alpha = 0.7)+
  theme(legend.position = "bottom")

g2 = ydata %>%
  filter(method == "snf") %>%
  ggplot(aes(V1,V2, color = factor(true_label))) +
  geom_point(size = 1, alpha = 0.7)+
  theme(legend.position = "bottom")

g3 = ydata %>%
  filter(method == "simlr") %>%
  ggplot(aes(V1,V2, color = factor(true_label))) +
  geom_point(size = 1, alpha = 0.7)+
  theme(legend.position = "bottom")

library(patchwork)
g1+g2+g3
