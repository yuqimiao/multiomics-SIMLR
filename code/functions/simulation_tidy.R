## this function is to combine all the sub_jobs

job_combine = function(path, pattern, scenario_add = FALSE){
  file_name = list.files(path)[str_detect(list.files(path), pattern)]
  table_all = NULL
  for(i in 1:length(file_name)){
    table = get(load(paste(path, file_name[i], sep = "")))
    if(scenario_add){
      table = table %>%
        mutate(scenario = str_extract(file_name[[i]], pattern = paste("n[A-Za-z0-9].+?(?=_",pattern, ")",sep = ""))) %>%
        select(scenario,everything())
    }
    table_all = rbind(table_all, table)
  }
  return(table_all)
}
