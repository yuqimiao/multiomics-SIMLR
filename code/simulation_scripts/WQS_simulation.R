# setwd("/ifs/scratch/msph/biostat/sw2206/yuqi")
# .libPaths("/ifs/scratch/msph/biostat/sw2206/yuqi/R_libs")
# task_ID = as.integer(Sys.getenv("SGE_TASK_ID"))
# N= 1:100
# n = N[task_ID]

# association simulation
library(wqs)
library(bkmr)
# Data simulation
M = 10
ind = c(1,2)
inter_ind = 1:2
threshold = 0.1
beta1 = 10
beta_int = 1/2
# correlation structure (pre-specify)
Sigma <- diag(1, M, M)
Sigma[1, 3] <- Sigma[3, 1] <- 0.95
Sigma[2, 3] <- Sigma[3, 2] <- 0.3
Sigma[1, 2] <- Sigma[2, 1] <- 0.1

SimData2 = function (n = 100, M = 5, sigsq.true = 0.5, beta.true = 2, hfun = 3,
                     Zgen = "corr", Sigma, ind = 1:2,inter_ind = 1:2,
                     family = "gaussian", beta1 = 1, beta_int = 1/2) {
  stopifnot(n > 0, M > 0, sigsq.true >= 0, family %in% c("gaussian",
                                                         "binomial"))
  if (family == "binomial") {
    sigsq.true <- 1
  }

  HFun = function (z, ind, inter_ind, beta1, beta_int) {
    # the interaction term can only take 2 variables at 1 function
    4 * plogis((beta1 * (sum(z[ind])) + beta_int * (z[inter_ind[1]] * z[inter_ind[2]])),
    0, 0.3)
  }

  if (Zgen == "unif") {
    Z <- matrix(runif(n * M, -2, 2), n, M)
  }else if (Zgen == "norm") {
    Z <- matrix(rnorm(n * M), n, M)
  }else if (Zgen == "corr") {
    if (M < 3) {
      stop("M must be an integer > 2 for Zgen = 'corr'")
    }
    Z <- MASS::mvrnorm(n = n, mu = rep(0, M), Sigma = Sigma)
  }

  colnames(Z) <- paste0("z", 1:M)
  X <- cbind(3 * cos(Z[, 1]) + 2 * rnorm(n))
  eps <- rnorm(n, sd = sqrt(sigsq.true))
  h <- apply(Z, 1, function(x) {HFun(x, ind = ind, inter_ind = inter_ind, beta1 = beta1, beta_int = beta_int)})
  mu <- X * beta.true + h
  y <- drop(mu + eps)
  if (family == "binomial") {
    ystar <- y
    y <- ifelse(ystar > 0, 1, 0)
  }
  dat <- list(n = n, M = M, sigsq.true = sigsq.true, beta.true = beta.true,
              Z = Z, h = h, X = X, y = y, hfun = hfun, HFun = HFun,
              family = family)
  if (family == "binomial") {
    dat$ystar <- ystar
  }
  dat
}

data2 = SimData2(M = M, ind = ind, inter_ind = inter_ind,Sigma = Sigma,beta1 = beta1,beta_int = beta_int)
cor(data2$Z)

repeat_time = 10
fit_weights = NULL
fit_res = list()
num_correct = NULL
num_incorrect = NULL
for(i in 1:repeat_time){
  fitwqs = wqs.est(y.train = data2$y, x.train = data2$Z,z.train = data2$X)
  fit_weights = rbind(fit_weights, fitwqs$weights)
  fit_res[[i]] = summary(fitwqs$fit)
  selected = fitwqs$weights[fitwqs$weights>threshold]
  num_correct = c(num_correct, sum(sapply(1:length(ind), FUN = function(x) {sum(paste("w",ind[x], sep = "") == names(selected))})))
  num_incorrect = c(num_incorrect, length(selected) - num_correct[i])
}
