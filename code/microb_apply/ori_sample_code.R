library('MiSPU')

# load data
data(throat.tree)
data(dd)

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

n = 100 # sample size
OTUtab = sim_data(nSam = n) # generate OTU table

# generate some binray outcomes
info_OTUs = names(which(clustering == 14)) # consider the OTUs in the largest cluster as informative OTUs
beta = rep(1, length(info_OTUs)) # effect size
scaled_OTUtab = OTUtab/rowSums(OTUtab)
eta = scale(scaled_OTUtab[, info_OTUs] %*% beta)
prob = 1/(1 + exp(-eta)) 
Y = rbinom(n, 1, prob) # final binary outcome

# some distance based on OTUs and tree
dist_ls = GUniFrac(OTUtab, throat.tree)
