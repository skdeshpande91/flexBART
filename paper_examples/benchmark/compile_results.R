load("../bakeoff_simulation/sim_settings.RData")

n_sim <- 50
methods <- c("flexBART", "flexDART",
             "BART1", "DART1", "targetBART", "targetDART")
tmp_results <-
  data.frame(matrix(nrow = n_sim,ncol = 6, dimnames = list(c(), methods)))

smse_train <- list()
for(d in data_list){
  smse_train[[d]] <- tmp_results
}
smse_test <- smse_train
cov_train <- smse_train
cov_test <- smse_train
timing <- smse_train


for(job_id in 1:nrow(sim_settings)){
  
  file_name <- paste0("bakeoff_results_", job_id, ".RData")
  if(file.exists(file_name)){
    load(file_name)
    tmp_results <- get(paste0("bakeoff_", job_id))
    sim_number <- tmp_results$sim_number
    d <- tmp_results$dataset_name
    
    smse_train[[d]][sim_number,] <- tmp_results$smse_train
    smse_test[[d]][sim_number,] <- tmp_results$smse_test
    cov_train[[d]][sim_number,] <- tmp_results$cov_train
    cov_test[[d]][sim_number,] <- tmp_results$cov_test
    timing[[d]][sim_number,] <- tmp_results$timing
    
    rm(list = paste0("bakeoff_", job_id))
    rm(tmp_results)
  }
}
save(smse_train, smse_test, cov_train, cov_test, timing, file = "bakeoff_results.RData")
t_test_bart <- rep(NA, times = length(data_list))
names(t_test_bart) <- data_list
t_test_dart <- rep(NA, times = length(data_list))
names(t_test_dart) <- data_list

improve_bart <- rep(NA, times = length(data_list))
names(improve_bart) <- data_list

improve_dart <- rep(NA, times = length(data_list))
names(improve_dart) <- data_list

for(d in data_list){
  t_test_bart[d] <- t.test(smse_test[[d]][,"flexBART"],
                           smse_test[[d]][,"BART1"],
                           paired = TRUE, alternative = "less")$p.value
  t_test_dart[d] <- t.test(smse_test[[d]][,"flexDART"],
                           smse_test[[d]][,"DART1"],
                           paired = TRUE, alternative = "less")$p.value
  improve_bart[d] <- sum(smse_test[[d]][,"flexBART"] < smse_test[[d]][,"BART1"])
  improve_dart[d] <- sum(smse_test[[d]][,"flexDART"] < smse_test[[d]][,"DART1"])
  
}



