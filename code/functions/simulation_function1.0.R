# simulation1
# * 2 types of data: GE and MI
# * 150 samples
# * 1000 features with 5% effective variables
# * 3 subtypes:A,B,C, each 50 separately
# * criteria for effective var
# * GE -- split A and BC:
#   * GE_A ~ N(-1,2)
# * GE_B, GE_C ~ N(1,2)
# * MI -- split AB and C:
#   * MI_A,MI_B ~ N(1,2)
# * MI_C~ N(-1,2)
# noise comes from N(0,2)

# function: simulation_1
# perform 1 simulation
# input:
#   1. x -- sample size
#   2. eff_ratio -- % of effective variables within all features
#   3. eff_size -- the effective size of the effective variables
#   4. sub_ratio -- a vector of the ratio of each subtype within sample
# output:
#   GE and MI, rows are samples, columns are features


# 3 data type simulation
# adjust number of subtype based on given ratio

simulation_3 = function(size = 200,sub_ratio =rep(0.25,4),eff_size = rep(5,3),sigma = rep(100,3),data_divide =c("12/34","13/24","14/23"),dist = rep("normal",3),n1 = 50*0.6,n2 = 950*0.8,r1 = 1,r2 = 3,uninfo_r1 = 0,uninfo_r2 = 1){
  ## calculate index_cut for given subtype ratio
  index_cut = cumsum(sub_ratio)*size
  index = list()
  for(i in 1:length(index_cut)){
    if(i == 1){
      index[[i]] = 1:index_cut[1]
    }else{
      index[[i]] = (index_cut[[i-1]]+1):index_cut[[i]]
    }
  }

  data1 = get_data(size, index, data_divide[1], dist[1],eff_size[1],sigma[1])
  data2 = get_data(size, index, data_divide[2], dist[2],eff_size[2],sigma[2])
  data3 = get_data(size, index, data_divide[3], dist[3],eff_size[3],sigma[3])
  weight1 = get_weight(n1 = n1,n2 = n2,r1 = r1,r2 = r2,uninfo_r1 = uninfo_r1,uninfo_r2 = uninfo_r2)
  weight2 = get_weight(n1 = n1,n2 = n2,r1 = r1,r2 = r2,uninfo_r1 = uninfo_r1,uninfo_r2 = uninfo_r2)
  weight3 = get_weight(n1 = n1,n2 = n2,r1 = r1,r2 = r2,uninfo_r1 = uninfo_r1,uninfo_r2 = uninfo_r2)
  truelabel = NULL
  for (i in 1:length(index)){
    truelabel = c(truelabel, rep(i,length(index[[i]])))
  }
  return(list(data1 = data1, data2 = data2, data3 = data3, weight1 = weight1, weight2 = weight2, weight3 = weight3,truelabel = truelabel) )
}

simulation_2 = function(size = 150,sub_ratio =  rep(1/3,3),eff_size = rep(1,2),sigma = rep(2^2,2),data_divide =c("1/23","12/3"),dist = rep("normal",2),n1 = 50*0.6,n2 = 950*0.8,r1 = 1,r2 = 3,uninfo_r1 = 0,uninfo_r2 = 1){
  ## calculate index_cut for given subtype ratio
  index_cut = cumsum(sub_ratio)*size
  index = list()
  for(i in 1:length(index_cut)){
    if(i == 1){
      index[[i]] = 1:index_cut[1]
    }else{
      index[[i]] = (index_cut[[i-1]]+1):index_cut[[i]]
    }
  }

  data1 = get_data(size, index, data_divide[1], dist[1],eff_size[1],sigma[1])
  data2 = get_data(size, index, data_divide[2], dist[2],eff_size[2],sigma[2])

  weight1 = get_weight(n1 = n1,n2 = n2,r1 = r1,r2 = r2,uninfo_r1 = uninfo_r1,uninfo_r2 = uninfo_r2)
  weight2 = get_weight(n1 = n1,n2 = n2,r1 = r1,r2 = r2,uninfo_r1 = uninfo_r1,uninfo_r2 = uninfo_r2)

  truelabel = NULL
  for (i in 1:length(index)){
    truelabel = c(truelabel, rep(i,length(index[[i]])))
  }
  return(list(data1 = data1, data2 = data2, weight1 = weight1, weight2 = weight2,truelabel = truelabel) )
}


get_data = function(size, index, data_divide, dist, eff_size,sigma ){
  data1 = matrix(0,size, 1000)
  divide = lapply(X = str_split(str_split(data_divide,pattern = "/")[[1]], pattern = ""),
                  FUN = as.numeric
  )
  ind1 = NULL
  for(i in 1:length(divide[[1]])){
    ind1 = c(ind1, index[[divide[[1]][i]]])
  }
  ind2 = NULL
  for(i in 1:length(divide[[2]])){
    ind2 = c(ind2, index[[divide[[2]][i]]])
  }
  if (dist == "normal"){
    data1[ind1,1:50] = rnorm(length(ind1)*50, eff_size[1], sqrt(sigma))
    data1[ind2,1:50] = rnorm(length(ind2)*50, -eff_size[1], sqrt(sigma))
    data1[,51:1000] = rnorm(size*950, 0, sqrt(sigma))
  }else if(dist == "logit"){
    data1[ind1,1:50] = sapply(X = rnorm(length(ind1)*50, eff_size[1], sqrt(sigma)),
                              FUN = function(x) exp(x)/(1+exp(x)))
    data1[ind2,1:50] = sapply(X = rnorm(length(ind2)*50, -eff_size[1], sqrt(sigma)),
                              FUN = function(x) exp(x)/(1+exp(x)))
    data1[,51:1000] = sapply(X = rnorm(size*950, 0, sqrt(sigma)),
                              FUN = function(x) exp(x)/(1+exp(x)))

  }
  return(data1)
}
#n1 is # of correct weights for informative features
#n2 is # of correct weights for non-informative features
#r1, r2, uninfo_r1, uninfo_r2 define the magnitude of correct wights
get_weight = function(n1 = 50*0.6,n2 = 950*0.8,r1 = 1,r2 = 3,uninfo_r1 = 0,uninfo_r2 = 1){
  ge_weight=abs(runif(1000,0,1))
  s=sample(1:50,n1)
  ge_weight[s]=runif(n1,r1,r2)
  ss=sample(51:1000,n2)
  ge_weight[ss]=runif(n2,uninfo_r1,uninfo_r2)
  return(ge_weight)
}


















