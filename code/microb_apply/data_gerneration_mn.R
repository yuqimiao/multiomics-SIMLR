# load data
library(MiSPU)
library(tidyverse)
library(rlist)
data(throat.tree)
data(dd)
throat.tree$tip.label = paste0('OTU', throat.tree$tip.label)
names(dd$pi) = paste0('OTU', names(dd$pi))
names(dd$theta) = paste0('OTU', names(dd$theta))

# do clustering
ncluster = 20
clustering = pam(as.dist(cophenetic(throat.tree)), ncluster, diss = TRUE)$clustering
p_est = dd$pi # mle of prob
p_est = p_est[names(clustering)]
p_clus = tapply(p_est, clustering, sum) # abundance level of the 20 clusters
# OTU info and distribution for each cluster
OTU_info = tibble(number = table(clustering),abundance = p_clus) %>%
  mutate(cluster = seq_along(number)) %>%
  dplyr::select(cluster, everything())

## original distribution
OTU_tibble = tibble(p_est, clustering) %>%
  mutate(OTU = names(p_est))

OTU_tibble %>%
  ggplot(aes(x = factor(clustering), y = p_est))+
  geom_boxplot()

## narrower distribution
OTU_tibble %>%
  filter(p_est <=0.002) %>%
  ggplot(aes(x = factor(clustering), y = p_est))+
  geom_boxplot()

## some parameters for distribution of OTU counts
theta = dd$theta
gplus = (1-theta)/theta
g_est = p_est*gplus

# a function to generate OTU table
sim_data = function(nSam = 100, mu = 1000, size = 25)
{
  comm = matrix(0, nSam, length(g_est))
  comm.p = comm
  rownames(comm) = 1:nrow(comm)
  colnames(comm) = names(g_est)
  nSeq = rnbinom(nSam, mu = mu, size = size)
  for (i in 1:nSam) {
    comm.p[i, ] = MiSPU::rdirichlet(1, g_est)[1, ]
    comm[i, ] = rmultinom(1, nSeq[i], prob = comm.p[i, ])[,1]
  }
  return(comm)
}
# simulate data
set.seed(1234)
OTUtab = sim_data(nSam = 10^4)
## check non-0 features in each cluster
zero_prop = lapply(1:20,FUN = function(i){
  sapply(names(clustering)[clustering == i], function(x){
    sum(OTUtab[,x] == 0)/10000
  })
})
OTU_info = OTU_info %>%
  mutate(zero_prop = unlist(lapply(zero_prop, mean))) %>%
  arrange(desc(abundance))
view(OTU_info)
## for unrelated OTU with pre/abs
un_abu_feature = list()
sort_OTU_cluster = clustering[names(sort(p_est,decreasing = T))]
flag = 0
i = 1
for(m in 1:10){
  if(!flag){
    un_abu_feature[[m]] = names(sort_OTU_cluster)[i]
    i = i+1
  }
  if(i>1){
    mask1 = (sort_OTU_cluster[i] %in% sort_OTU_cluster[1:(i-1)])
    mask2 = sort_OTU_cluster[i] %in% c(5,7)
    if(mask1|mask2){flag = 1}
  }
  while(flag){
    i = i+1
    mask1 = (sort_OTU_cluster[i] %in% sort_OTU_cluster[1:(i-1)])
    mask2 = sort_OTU_cluster[i] %in% c(5,7)
    if(!mask1&!mask2){flag = 0}
  }
}
un_abu_feature = unlist(un_abu_feature)
#################################################################################################################
### Set info OTUs
info_OTUs_ls = vector('list', length = 3) # list of info OTUs
names(info_OTUs_ls) = c('Phy-related_abund', 'Phy-related_pres', 'Phy-unrelated_abun')

# 1st set of info OTUs: the 7th largest cluster. We will use abundance of these OTUs
info_OTUs_ls[[1]] = names(which(clustering == names(p_clus)[7]))
length(info_OTUs_ls[[1]]) # 57 OTUs
sum(p_est[info_OTUs_ls[[1]]]) # abundance level 0.1038962

# 2nd set of info OTUs: the 5th largest cluster. We will use presence/absence of these OTUs
info_OTUs_ls[[2]] = names(which(clustering == names(p_clus)[5]))
length(info_OTUs_ls[[2]]) # 29 OTUs
sum(p_est[info_OTUs_ls[[2]]]) # abundance level 0.046

# 3rd set of info OTUs: some phylo-unrelated OTUs with large abundances not in the previous 2 sets
# We will use abundance of these OTUs
info_OTUs_ls[[3]] = un_abu_feature
length(info_OTUs_ls[[3]]) # 10 OTUs
sum(p_est[info_OTUs_ls[[3]]]) # abundance level 0.2544
length(unique(clustering[info_OTUs_ls[[3]]])) # the 10 info OTUs are scattered in 10 clusters, so phylo-unrelated and also not in 5/7

### Generate data
set.seed(123)
OTUnames = throat.tree$tip.label
dataset = matrix(NA, nrow = 10^4, ncol = length(OTUnames)+1)
colnames(dataset) = c(OTUnames, 'Subtype')

################################################################################

## set the betas for each category
abun_generate = function(beta0 = -3,beta_group = c(80,4,65)){
  print(c(beta0,beta_group))
  beta_list = list()
  for (k in 1:3) {
    beta = rep(0, length(OTUnames))
    names(beta) = OTUnames
    beta[info_OTUs_ls[[k]]] = rep(beta_group[k], length(info_OTUs_ls[[k]]),replace = T)
    # for(j in (1:3)[-k]){
    #   beta[info_OTUs_ls[[j]]] = rep(-5, length(info_OTUs_ls[[j]]))
    # }
    # beta0[info_OTUs_ls[[k]]] = sample(c(beta_group[k],-beta_group[k]), length(info_OTUs_ls[[k]]), replace=T)
    beta_list[[k]] = beta
  }

  ## Calculation of category probability for every individual
  eta_list = NULL
  for (k in 1:3) {
    if (k == 2) {
      OTUtab[OTUtab>0] = 1 # change abundance to presence information
      input_OTUtab = OTUtab
    } else {
      input_OTUtab = OTUtab/rowSums(OTUtab) # convert OTU counts to abundance
    }
    # input_OTUtab = OTUtab/rowSums(OTUtab)
    # input_OTUtab_sc = scale(input_OTUtab)
    input_OTUtab_sc = input_OTUtab - mean(input_OTUtab) # convert OTU counts to abundance
    eta_list = cbind(eta_list, input_OTUtab_sc %*% beta_list[[k]] + beta0)
  }
  par(mfrow = c(1,3))
  hist(input_OTUtab_sc[,info_OTUs_ls[[1]]],main = "subtype1", xlab = "centered OTU%")
  hist(input_OTUtab_sc[,info_OTUs_ls[[2]]],main = "subtype2", xlab = "centered OTU%")
  hist(input_OTUtab_sc[,info_OTUs_ls[[3]]],main = "subtype3", xlab = "centered OTU%")

  par(mfrow = c(1,3))
  hist(eta_list[,1],main = "subtype1", xlab = "eta")
  hist(eta_list[,2],main = "subtype2", xlab = "eta")
  hist(eta_list[,3],main = "subtype3", xlab = "eta")


  expeta_tibble = as_tibble(eta_list) %>%
    mutate_all(exp) %>%
    mutate(V4 = 1,
           denom = 1+V1+V2+V3
    )
  P_tibble = apply(expeta_tibble[1:4], MARGIN = 2, FUN = function(x) x/expeta_tibble$denom)
  colMeans(P_tibble)
  Y = apply(P_tibble, MARGIN = 1, FUN = function(x) which.max(x))
  print(table(Y))
  dataset = data.frame(cbind(OTUtab,Y))
  colnames(dataset) = c(OTUnames, 'Subtype')
  dataset_abun = dataset # the dataset with OTU abundances
  dataset_abun[, OTUnames] = dataset_abun[, OTUnames]/rowSums(dataset_abun[, OTUnames])

  return(list(dataset=dataset, dataset_abun = dataset_abun))
}

## the parameter of betas should be 2 levels
abun_generate_2 = function(beta0 = -7,beta_group = c(80,4,65)){
  print(c(beta0,beta_group))
  beta_list = list()
  # 2 layers of beta, beta_kt
  #k is for the beta of different types
  for (k in 1:3) {
    # t is for betas of different categories
    beta_categ = list()
    for(t in 1:3){
      beta = rep(0, length(info_OTUs_ls[[k]]))
      names(beta) = info_OTUs_ls[[k]]
      beta[info_OTUs_ls[[k]]] = beta_group[k]*rnorm(length(info_OTUs_ls[[k]]),mean = 1, sd = 10)
      beta_categ[[t]] = beta
    }

    # for(j in (1:3)[-k]){
    #   beta[info_OTUs_ls[[j]]] = rep(-5, length(info_OTUs_ls[[j]]))
    # }
    # beta0[info_OTUs_ls[[k]]] = sample(c(beta_group[k],-beta_group[k]), length(info_OTUs_ls[[k]]), replace=T)
    beta_list[[k]] = beta_categ
  }
  beta_all = list.cbind(lapply(1:3,function(x){unlist(lapply(1:3, function(y) {beta_list[[y]][[x]]}))}))
  ## Calculation of category probability for every individual
  eta_list = list()
  input_OTUtab = OTUtab/rowSums(OTUtab)
  OTU_pa = OTUtab
  OTU_pa[OTU_pa>0] = 1
  input_OTUtab[,info_OTUs_ls[[2]]] = OTU_pa[,info_OTUs_ls[[2]]] # change abundance to presence information
  input_OTUtab_sc = input_OTUtab - rowMeans(input_OTUtab) # convert OTU counts to abundance

  eta_list = input_OTUtab_sc[,unlist(info_OTUs_ls)] %*% as.matrix(beta_all)+beta0

  expeta_tibble = as_tibble(eta_list) %>%
    mutate_all(exp) %>%
    mutate(V4 = 1,
           denom = 1+V1+V2+V3
    )
  P_tibble = apply(expeta_tibble[1:4], MARGIN = 2, FUN = function(x) x/expeta_tibble$denom)
  colMeans(P_tibble)
  Y = apply(P_tibble, MARGIN = 1, FUN = function(x) which.max(x))
  print(table(Y))
  dataset = data.frame(cbind(OTUtab,Y))
  colnames(dataset) = c(OTUnames, 'Subtype')
  dataset_abun = dataset # the dataset with OTU abundances
  dataset_abun[, OTUnames] = dataset_abun[, OTUnames]/rowSums(dataset_abun[, OTUnames])

  return(list(dataset=dataset, dataset_abun = dataset_abun))
}

# dataset_abun = as_tibble(dataset_abun)
# save(dataset_abun, file = "microb_multisim.Rdata")

pv_calc = function(dataset_abun){
  abun_tibble = as_tibble(dataset_abun)
  abun_tibble_pv = abun_tibble %>%
    mutate(Subtype = ifelse(Subtype == 4,0,1)) %>%
    pivot_longer(
      OTU1883:OTU2582,
      names_to = "OTU",
      values_to = "abundance"
    ) %>%
    nest(data = c(Subtype,abundance)) %>%
    mutate(model = map(data, ~lm(Subtype~abundance, data = .x)),
           model = map(model, broom::tidy)) %>%
    unnest(model)

  info_rank = abun_tibble_pv %>% filter(term == "abundance") %>%
    arrange(p.value) %>%
    mutate(order = seq_along(term)) %>%
    filter(OTU %in% unlist(info_OTUs_ls)) %>%
    summarise(info_rank_mean = mean(order),
              info_rank_var = var(order),
              info_rank_media = median(order))

  pv = abun_tibble_pv%>% filter(term == "abundance") %>% pull(p.value)
  return(list(pv = pv, info_rank = info_rank))

}

################################################################################


