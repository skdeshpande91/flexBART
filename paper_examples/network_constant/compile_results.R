load("sim_settings.RData")

get_name <- function(results){
  method_name <- NULL
  if(results$alg == "networkBART"){
    method_name <- paste0("networkBART_", results$graph_cut_type)
  } else if(results$alg == "flexBART"){
    method_name <- "flexBART"
  } else if(results$alg == "BART"){
    if(is.na(results$embed)){
      if(results$sparse) method_name <- "DART"
      else method_name <- "BART"
    } else{
      if(results$sparse){
        method_name <- paste0("DART_ase", results$embed_dim)
      } else{
        method_name <- paste0("BART_ase", results$embed_dim)
      }
    }
  }
  return(method_name)
}
all_names <- c(paste0("networkBART_", 1:7),
               "flexBART",
               "BART", "DART",
               paste0("BART_ase", c(1,3,5)),
               paste0("DART_ase", c(3,5)))

n_methods <- length(all_names)
n_sim <- 50
tmp_results <- 
  data.frame(matrix(nrow = n_sim, ncol = n_methods, dimnames = list(c(), all_names)))

fit_mse_train <- list(n_obs10 = tmp_results,
                      n_obs100 = tmp_results)
fit_mse_test <- fit_mse_train
fit_cov_train <- fit_mse_train
fit_cov_test <- fit_mse_train
timing <- fit_mse_train

for(job_id in 1:nrow(sim_settings)){
  if(job_id %% floor(nrow(sim_settings)/10) == 0) print(paste0("Job", job_id, "of", nrow(sim_settings)))
  file_name <- paste0("network_constant1_results_", job_id, ".RData")
  if(file.exists(file_name)){
    load(file_name)
    results <- get(paste0("network_constant1_", job_id))
    method_name <- get_name(results)
    n_obs_name <- paste0("n_obs", results$n_obs)
    sim_number <- results$sim_number
    if(!is.null(method_name) & method_name %in% all_names){
      timing[[n_obs_name]][sim_number, method_name] <- results$timing
      fit_mse_train[[n_obs_name]][sim_number, method_name] <- results$fit_mse_train
      fit_mse_test[[n_obs_name]][sim_number, method_name] <- results$fit_mse_test
      
      fit_cov_train[[n_obs_name]][sim_number, method_name] <- results$fit_cov_train
      fit_cov_test[[n_obs_name]][sim_number, method_name] <- results$fit_cov_test
    } else{
      print(paste0("Could not resolve name for", job_id))
    }
    rm(list = paste0("network_constant1_", job_id))
    rm(results)
  } else{
    print(paste("Job ID =", job_id, "not found"))
  }
}

save(fit_mse_train, fit_mse_test,
     fit_cov_train, fit_cov_test,
     timing, file = "network_constant1_results.RData")

  


