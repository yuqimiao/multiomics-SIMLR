library('MiSPU')

# load data
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
p_clus = sort(tapply(p_est, clustering, sum), decreasing = T) # abundance level of the 20 clusters

# some parameters for distribution of OTU counts
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

n = 300 # sample size
OTUtab = sim_data(nSam = n) # generate OTU table

# generate some binray outcomes
info_OTUs = names(which(clustering == 14)) # consider the OTUs in the largest cluster as informative OTUs
beta = rep(1, length(info_OTUs)) # effect size
scaled_OTUtab = OTUtab/rowSums(OTUtab)
eta = scale(scaled_OTUtab[, info_OTUs] %*% beta)
prob = 1/(1 + exp(-eta))
Y = rbinom(n, 1, prob) # final binary outcome

# some distance based on OTUs and tree
dist_ls = MiSPU::GUniFrac(OTUtab[1:5,], throat.tree)


#################################################################################################################
### Set info OTUs
info_OTUs_ls = vector('list', length = 3) # list of info OTUs
names(info_OTUs_ls) = c('Phy-related_abund', 'Phy-related_pres', 'Phy-unrelated_abun')

# 1st set of info OTUs: the 2nd largest cluster. We will use abundance of these OTUs
info_OTUs_ls[[1]] = names(which(clustering == names(p_clus)[2]))
length(info_OTUs_ls[[1]]) # 57 OTUs
sum(p_est[info_OTUs_ls[[1]]]) # abundance level 0.1038962

# 2nd set of info OTUs: the 12th largest cluster. We will use presence/absence of these OTUs
info_OTUs_ls[[2]] = names(which(clustering == names(p_clus)[12]))
length(info_OTUs_ls[[2]]) # 13 OTUs
sum(p_est[info_OTUs_ls[[2]]]) # abundance level 0.03

# 3rd set of info OTUs: some phylo-unrelated OTUs with large abundances not in the previous 2 sets
# We will use abundance of these OTUs
info_OTUs_ls[[3]] = names(sort(p_est[!names(p_est) %in% unlist(info_OTUs_ls[1:2])], decreasing = T)[11:20])
length(info_OTUs_ls[[3]]) # 10 OTUs
sum(p_est[info_OTUs_ls[[3]]]) # abundance level 0.1244
length(unique(clustering[info_OTUs_ls[[3]]])) # the 10 info OTUs are scattered in 8 clusters, so phylo-unrelated

### Generate data
set.seed(123)
OTUnames = throat.tree$tip.label
dataset = matrix(NA, nrow = 10^4, ncol = length(OTUnames)+1)
colnames(dataset) = c(OTUnames, 'Subtype')
for (k in 1:3) {
  info_OTUs = info_OTUs_ls[[k]]
  beta = rep(1, length(info_OTUs)) # effect size for each OTU
  beta_group = 1 # effect size for the set of OTUs
  OTUtab = sim_data(nSam = 10^4)
  if (k == 2) {
    OTUtab[OTUtab>0] = 1 # change abundance to presence information
    input_OTUtab = OTUtab
  } else {
    input_OTUtab = OTUtab/rowSums(OTUtab) # convert OTU counts to abundance
  }
  eta = scale(input_OTUtab[, info_OTUs] %*% beta)
  prob = 1/(1 + exp(-eta * beta_group))
  Y = rbinom(10^4, 1, prob) # final binary outcome

  case_id = sample(which(Y==1), 3333+floor((k-1)/2), replace = F) # find 3333 cases (or 3334 cases for the last set of info OTUs)
  control_id = sample(which(Y==0), 3333+floor((k-1)/2), replace = F)
  dataset[((k-1)*3333+1):((k-1)*3333+length(case_id)),] = cbind(OTUtab[case_id,], rep(k, length(case_id)))
}
dataset = data.frame(dataset) # the dataset with OTU counts
dataset_abun = dataset # the dataset with OTU abundances
dataset_abun[, OTUnames] = dataset_abun[, OTUnames]/rowSums(dataset_abun[, OTUnames])

### Get microbial distances
library('GUniFrac')
UniFrac_BC = function(OTUtab = OTUtab) {return(as.matrix(stats::dist(OTUtab, method = 'manhattan'))/2)}

sample_id = 1:100
unifracs = GUniFrac::GUniFrac(otu.tab = dataset_abun[sample_id, OTUnames], throat.tree, alpha = 1)$unifracs

D_BC = UniFrac_BC(dataset_abun[sample_id, OTUnames]) # Bray-Curtis distance
D_U = unifracs[,,"d_UW"] # Unweighted UniFrac distance
D_W = unifracs[,,"d_1"] # Weighted UniFrac distance

### Some important variables
info_OTUs_ls
dataset_abun[c(1:5, 4001:4005, 8001:8005), c(1:5, 857)]

## saving

# save(dataset_abun, file = "micro_case_pool2.Rdata")

