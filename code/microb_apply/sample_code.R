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
  nSeq = rnbinom(nSam, mu = mu, size = size) ## using rbinom to obtain a 25
  for (i in 1:nSam) {
    comm.p[i, ] = MiSPU::rdirichlet(1, g_est)[1, ]
    comm[i, ] = rmultinom(1, nSeq[i], prob = comm.p[i, ])[,1]
  }
  return(comm)
}

n = 600 # sample size: pnly take 50% as cases and divide nto 3 clusters
c = 3
OTUtab = sim_data(nSam = n) # generate OTU table

## Creat different cases based on different signals set
info_OTUs = names(which(clustering == 14))
scaled_OTUtab = OTUtab/rowSums(OTUtab)

### the coefficient of info_OTUs
set.seed(1)
signal_cut = floor(rep(1/c,c)*length(info_OTUs))

if((length(info_OTUs)%% c)){
  for(i in 1:(length(info_OTUs)%% c)){
    signal_cut[i] = signal_cut[i]+1
  }
}

assign = rep(0,ncol(scaled_OTUtab))
assign[sample(which(clustering == 14))] = rep(1:c, signal_cut)
beta_list = list()
pool = cbind(matrix(0,n,c),scaled_OTUtab)
colnames(pool) = c(paste("y",1:c, sep = ""),colnames(scaled_OTUtab))

for(i in 1:c){
  beta = rep(0, ncol(scaled_OTUtab))
  beta[assign ==i] = 5
  eta = scale(scaled_OTUtab %*% beta)
  prob = 1/(1 + exp(-eta))
  pool[seq((i-1)*n/c+1, (i*n/c), by = 1),i] = rbinom(n/c,1,prob)
}


as_tibble(pool) %>%
  pivot_longer(
    y1:y3,
    names_to = "class",
    values_to = "status"
  ) %>%
  group_by(class) %>%
  summarise(case = sum(status))

############
## case extraction with the scaled OTUtable

# generate some binray outcomes
info_OTUs = names(which(clustering == 14)) # consider the OTUs in the largest cluster as informative OTUs

beta = rep(1, length(info_OTUs)) # effect size
scaled_OTUtab = OTUtab/rowSums(OTUtab)
eta = scale(scaled_OTUtab[, info_OTUs] %*% beta)
prob = 1/(1 + exp(-eta))
Y = rbinom(n, 1, prob) # final binary outcome
table(Y)

# some distance based on OTUs and tree
dist_ls = GUniFrac(OTUtab, throat.tree)
