## tidy file
file_names = list.files("./data/simu_0216")
nmi_tune = NULL
for(i in file_names){
  load(paste("./data/simu_0216/",i, sep = ""))
  nmi_tune = rbind(nmi_tune,unlist(nmi_ls))
}
nmi_tunning = as_tibble(nmi_tune) %>%
  mutate(pars = str_replace(str_replace(file_names,pattern = "(nmi_)|(.Rdata)",
                                        replacement = ""),".Rdata",""),
         ID = seq_along(pars)) %>%
  separate(pars, into = c("beta","gamma","rho","alpha0","alpha"),sep = "_") %>%
  pivot_longer(
    cols = beta:alpha,
    names_to = "pars",
    values_to = "value"
  ) %>%
  mutate(value = as.numeric(str_extract(value, "\\d?\\.?\\d+"))) %>%
  pivot_longer(
    cols = scen7:scen9,
    names_to = "scenario",
    values_to = "nmi"
  ) %>%
  dplyr::select(ID, scenario, pars, value, nmi) %>%
  arrange(ID)

## performance and par summary
nmi_tunning_wide = nmi_tunning%>%
  pivot_wider(
    names_from = pars,
    values_from = value
  )

# see the maximum for each scenario
best_table = nmi_tunning_wide %>%
  group_by(scenario) %>%
  filter(nmi == max(nmi)) %>%
  summarise(max_nmi = max(nmi),
            mean_beta = mean(beta),
            mean_gamma = mean(gamma),
            mean_rho = mean(rho),
            mean_alpha0 = mean(alpha0),
            mean_alpha = mean(alpha))

# see the high performance 10% for each scenario
ten_perc_table = nmi_tunning_wide %>%
  group_by(scenario) %>%
  filter(nmi >= quantile(nmi,0.9)) %>%
  summarise(mean_nmi = mean(nmi),
            mean_beta = mean(beta),
            mean_gamma = mean(gamma),
            mean_rho = mean(rho),
            mean_alpha0 = mean(alpha0),
            mean_alpha = mean(alpha))

View(best_table)
View(ten_perc_table)

## best par extraction
best_par_finding = function(x){
  nmi_tunning %>%
    group_by(scenario) %>%
    filter(nmi >= quantile(nmi, 0.9)) %>%
    ungroup(scenario) %>%
    pivot_wider(
      names_from = pars,
      values_from = value
    ) %>%
    filter(scenario == x) ## scen1
}

best_par = lapply(paste("scen",7:9,sep = ""), best_par_finding)
names(best_par) = paste("scen",7:9,sep = "")
best_par1 = map_dfr(best_par,function(x){x[1,]})
save(best_par, file = "best_tuned_par_scen789.Rdata")

load("best_tuned_par.Rdata")
best_par_all[["scen7"]] = best_par$scen7
best_par_all[["scen8"]] = best_par$scen8
save(best_par_all,file = "best_par1-9.Rdata")
