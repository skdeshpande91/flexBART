load("sim_settings.RData")
n <- 100

rmse_train <- data.frame(matrix(NA, nrow = n, ncol = length(method_list),
                                dimnames = list(c(), method_list)))


rmse_test <- rmse_train
coverage_train <- rmse_train
coverage_test <- rmse_train
timing <- rmse_train

for(job_id in 1:nrow(sim_settings)){
  file_name <- paste0("philly_obs_cv_results", job_id, ".RData")
  if(file.exists(file_name)){
    load(file_name)
    results <- get(paste0("philly_obs_cv_", job_id))
    
    method <- as.character(results$method)
    sim_number <- results$sim_number
    rmse_train[sim_number, method] <- results$rmse_train
    rmse_test[sim_number, method] <- results$rmse_test
    coverage_train[sim_number, method] <- results$cov_train
    coverage_test[sim_number, method] <- results$cov_test
    timing[sim_number, method] <- results$timing["elapsed"]
    rm(list = c("results", paste0("philly_obs_cv_", job_id)))
  } else{
    print(paste("Missing", job_id))
  }
}
save(rmse_train, rmse_test,
     coverage_train, coverage_test,
     timing, file = "philly_obs_cv_result.RData")

boxplot(rmse_test)


sum(rmse_test[,"networkBART"] < rmse_test[,"BART2"])
