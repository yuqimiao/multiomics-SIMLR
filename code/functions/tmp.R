a1 = 0.0001
a2 = 0.0005
b1 = seq(1,2,by = 0.1)
b2 = seq(1,3,by = 0.2)
w1_n = NULL
w2_n = NULL
for(i in 1:length(b1)){
  w1 = exp(a1+b1[i])
  w2 = exp(a2+b2[i])

  w1_n = c(w1_n,w1/(w1+w2))
  w2_n = c(w2_n,w2/(w1+w2))
}

a1 = 0.0001*100
a2 = 0.0005*100
b1 = seq(1,2,by = 0.1)
b2 = seq(1,3,by = 0.2)
w1_n1 = NULL
w2_n1 = NULL
for(i in 1:length(b1)){
  w1 = exp(a1+b1[i])
  w2 = exp(a2+b2[i])

  w1_n1 = c(w1_n1,w1/(w1+w2))
  w2_n1 = c(w2_n1,w2/(w1+w2))
}

## so the change of a1 and a2 will definetely change the basic value of the weight

data = simulation_3(size = par[[m]]$sample_size,
                    sigma = rep(par[[m]]$sigma,3),
                    eff_size = par[[m]]$eff_size)
truth = data$truelabel

W_g = GL_SIMLR_4(data_list = list(t(data$data1),t(data$data2),t(data$data3)),
                 c = 4,
                 beta = 0.1,
                 gamma = 1)
W_g1 = GL_SIMLR_4(data_list = list(t(data$data1),t(data$data2),t(data$data3)),
                  c = 4,
                  beta = 10,
                  gamma = 1)
W_g2 = GL_SIMLR_4(data_list = list(t(data$data1),t(data$data2),t(data$data3)),
                  c = 4,
                  beta = 100,
                  gamma = 1)

W_g3 = GL_SIMLR_4(data_list = list(t(data$data1),t(data$data2),t(data$data3)),
                 c = 4,
                 beta = 1,
                 gamma = 1,
                 F_norm = F,
                 standardize_term = T)
W_g4 = GL_SIMLR_4(data_list = list(t(data$data1),t(data$data2),t(data$data3)),
                  c = 4,
                  beta = 10,
                  gamma = 1,
                  F_norm = F,
                  standardize_term = T)

W_g5 = GL_SIMLR_4(data_list = list(t(data$data1),t(data$data2),t(data$data3)),
                  c = 4,
                  beta = 100,
                  gamma = 100,
                  F_norm = F,
                  standardize_term = T)
W_g3$w_list
W_g4$w_list
W_g5$w_list

compare(data$truelabel, W_g3$cluster, method = "nmi")
compare(data$truelabel, W_g4$cluster, method = "nmi")
compare(data$truelabel, W_g5$cluster, method = "nmi")

W_g_tunn$tunning_crit %>%
  pivot_longer(
    nmi:rand,
    names_to = "measure",
    values_to = "value"
  ) %>%
  ggplot(aes(x = rho, y = value, color = measure, group = measure))+
  geom_line()+
  facet_grid(beta~gamma)

W_g_tunn$tunning_crit %>%
  pivot_longer(
    w1:w2,
    names_to = "weight",
    values_to = "value"
  ) %>%
  ggplot(aes(x = rho, y = value, color = weight, group = weight))+
  geom_line()+
  facet_grid(beta~gamma)
