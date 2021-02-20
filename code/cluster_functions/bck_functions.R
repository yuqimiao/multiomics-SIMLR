simulation_2 = function(x_size = 150,eff_size = c(5,5), sub_ratio = rep(1/3,3),dist = c("norm", "norm"),sigma = 2,n1_r = 0.6,n2_r = 0.8,r1 = 1,r2 = 3,uninfo_r1 = 0,uninfo_r2 = 1){
  ##GE
  if(dist[1] == "norm"){
    GE = normal_sim(x_size, sub_ratio =c(sub_ratio[1],sub_ratio[2]+sub_ratio[3]),eff_size[1],sigma = sigma)
    #c(sub_ratio[1]+sub_ratio[2], sub_ratio[3])
  }else if(dist[1] == "logit"){
    GE = logit_sim(x_size, sub_ratio =c(sub_ratio[1],sub_ratio[2]+sub_ratio[3]),eff_size[1],sigma = sigma)
  }
  ##MI
  if(dist[2] == "norm"){
    MI = normal_sim(x_size, sub_ratio =c(sub_ratio[1]+sub_ratio[2], sub_ratio[3]),eff_size[2],sigma = sigma)
    
  }else if(dist[2] == "logit"){
    MI = logit_sim(x_size, sub_ratio =c(sub_ratio[1]+sub_ratio[2], sub_ratio[3]),eff_size[2],sigma = sigma)
  }
  
  ## feature weight for GE  
  eff_index = c(50,50)
  ge_weight=abs(runif(1000,0,1))
  s=sample(1:eff_index[1],n1_r*eff_index[1])
  ge_weight[s]=runif(n1_r*eff_index[1],r1,r2)
  ss=sample((eff_index[1]+1):1000,n2_r*eff_index[1])
  ge_weight[ss]=runif(n2_r*eff_index[1],uninfo_r1,uninfo_r2)
  
  mi_weight=abs(runif(1000,0,1))
  s=sample(1:eff_index[2],n1_r*eff_index[2])
  mi_weight[s]=runif(n1_r*eff_index[2],r1,r2)
  ss=sample((eff_index[2]+1):1000,n2_r*eff_index[2])
  mi_weight[ss]=runif(n2_r*eff_index[2],uninfo_r1,uninfo_r2)
  a = sub_ratio[1]*x_size
  b = (sub_ratio[1]+sub_ratio[2])*x_size
  true_label = c(rep(1,a),rep(2,(b-a)),rep(3,x_size-b))
  
  return(list(GE = GE, MI = MI,ge_weight = ge_weight, mi_weight = mi_weight,true_label = true_label))
}

normal_sim = function(x_size = 150, sub_ratio = c(2/3,1/3), eff_size = 5,sigma = 2){
  data1 = matrix(0, x_size, 1000)
  a = sub_ratio[1]*x_size
  sub1_index = 1:a
  sub2_index = (a+1):x_size
  data1[sub1_index,1:50] = rnorm(a*50, -eff_size, sqrt(sigma))
  data1[sub2_index,1:50] = rnorm((x_size-a)*50, eff_size, sqrt(sigma))
  data1[,51:1000] = rnorm(x_size*950, 0, sqrt(sigma))
  return(data1)
}

beta_sim = function(x_size = 150, sub_ratio = c(2/3,1/3), par = list(c(1,5),c(5,1)),sigma = 2){
  data1 = matrix(0, x_size, 1000)
  a = sub_ratio[1]*x_size
  sub1_index = 1:a
  sub2_index = (a+1):x_size
  data1[sub1_index,1:50] = rbeta(a*50, shape1 = par[[1]][1], shape2 = par[[1]][2]) ## hypo
  data1[sub2_index,1:50] = rbeta((x_size-a)*50, shape1 = par[[2]][1], shape2 = par[[2]][2])
  data1[,51:1000] = rbeta(x_size*950, 5,5, sqrt(sigma))
  return(data1)
}

logit_sim = function(x_size = 150, sub_ratio = c(2/3,1/3), eff_size = 5,sigma = 2){
  data1 = matrix(0, x_size, 1000)
  a = sub_ratio[1]*x_size
  sub1_index = 1:a
  sub2_index = (a+1):x_size
  data1[sub1_index,1:50] = sapply(rnorm(a*50, -eff_size, sqrt(sigma)), FUN = function(x) exp(x)/(1+exp(x)))
  data1[sub2_index,1:50] = sapply(rnorm((x_size-a)*50, eff_size, sqrt(sigma)),FUN = function(x) exp(x)/(1+exp(x)))
  data1[,51:1000] = sapply(rnorm(x_size*950, 0, sqrt(sigma)),FUN = function(x) exp(x)/(1+exp(x)))
  return(data1)
}