load("sim_settings.RData")
year_list <- unique(sim_settings$year)
method_list <- unique(sim_settings$method)

results_mat <- 
  data.frame(matrix(NA, nrow = 10, ncol = length(method_list),
                    dimnames = list(c(), method_list)))

brier_train <- list()
for(year in year_list) brier_train[[paste0("y.",year)]] <- results_mat
brier_test <- brier_train

logloss_train <- brier_train
logloss_test <- brier_train
misclass_train <- brier_train
misclass_test <- brier_train
timing <- brier_train


for(job_id in 1:nrow(sim_settings)){
  file_name <- paste0("pitchFraming_", job_id, ".RData")
  if(file.exists(file_name)){
    load(file_name)
    tmp_results <- get(paste0("pitchFraming_", job_id))
    method <- sim_settings[job_id,"method"]
    year <- tmp_results$year
    sim_number <- tmp_results$sim_number
    
    brier_train[[paste0("y.",year)]][sim_number, method] <- tmp_results$brier["train"]
    brier_test[[paste0("y.",year)]][sim_number, method] <- tmp_results$brier["test"]
    
    logloss_train[[paste0("y.",year)]][sim_number, method] <- tmp_results$logloss["train"]
    logloss_test[[paste0("y.",year)]][sim_number, method] <- tmp_results$logloss["test"]

    misclass_train[[paste0("y.",year)]][sim_number, method] <- tmp_results$misclass["train"]
    misclass_test[[paste0("y.",year)]][sim_number, method] <- tmp_results$misclass["test"]
    
    timing[[paste0("y.",year)]][sim_number, method] <- tmp_results$timing

  }
  
}
save(misclass_train, misclass_test,
     logloss_train, logloss_test,
     brier_train, brier_test,
     timing, file = "pitchFraming_results.RData")
all_misclass_test <- dplyr::bind_rows(misclass_test)
all_brier_test <- dplyr::bind_rows(brier_test)
all_timing <- dplyr::bind_rows(timing)
