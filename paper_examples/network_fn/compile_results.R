load("sim_settings.RData")
get_name <- function(results){
  method_name <- NULL
  if(results$alg == "networkBART"){
    if(results$sparse){
      method_name <- paste0("networkDART_", results$graph_cut_type)
    } else{
      method_name <- paste0("networkBART_", results$graph_cut_type)
    }
  } else if(results$alg == "flexBART"){
    if(results$sparse) method_name <- "flexDART"
    else method_name <- "flexBART"
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
               paste0("networkDART_", 1:7),
               "flexBART", "flexDART",
               "BART", "DART",
               paste0("BART_ase", c(1,3,5)),
               paste0("DART_ase", c(1,3,5)))

n_methods <- length(all_names)
n_sim <- 50
tmp_results <- 
  data.frame(matrix(nrow = n_sim, ncol = n_methods, dimnames = list(c(), all_names)))

mse_train <- list(n_obs10 = tmp_results,
                      n_obs100 = tmp_results)
mse_test1 <- mse_train
mse_test2 <- mse_train
cov_train <- mse_train
cov_test1 <- mse_train
cov_test2 <- mse_train
timing <- mse_train

for(job_id in 1:nrow(sim_settings)){
  if(job_id %% floor(nrow(sim_settings)/10) == 0) print(paste0("Job", job_id, "of", nrow(sim_settings)))
  file_name <- paste0("network_fn_results_", job_id, ".RData")
  if(file.exists(file_name)){
    load(file_name)
    results <- get(paste0("network_fn_", job_id))
    method_name <- get_name(results)
    n_obs_name <- paste0("n_obs",sim_settings[job_id, "n_obs"])
    #t_name <- paste0("t", results$t)
    sim_number <- results$sim_number
    if(!is.null(method_name) & method_name %in% all_names){
      timing[[n_obs_name]][sim_number, method_name] <- results$timing
      mse_train[[n_obs_name]][sim_number, method_name] <- results$mse_train
      mse_test1[[n_obs_name]][sim_number, method_name] <- results$mse_test1
      mse_test2[[n_obs_name]][sim_number, method_name] <- results$mse_test2
      cov_train[[n_obs_name]][sim_number, method_name] <- results$cov_train
      cov_test1[[n_obs_name]][sim_number, method_name] <- results$cov_test1
      cov_test2[[n_obs_name]][sim_number, method_name] <- results$cov_test2
    } else{
      print(paste0("Could not resolve name for", job_id))
    }
    rm(list = paste0("network_fn_", job_id))
    rm(results)
  } else{
    print(paste("Job ID =", job_id, "not found"))
  }
}

save(mse_train, mse_test1, mse_test2,
     cov_train, cov_test1, cov_test2,
     timing, file = "network_fn_results.RData")

  


