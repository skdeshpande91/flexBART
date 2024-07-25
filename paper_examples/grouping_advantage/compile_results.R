load("sim_settings.RData")

n_sim <- 50
methods <- c("flexBART", "flexDART",
             "BART", "DART",
             "oracleBART", "oracleDART",
             "targetBART", "targetDART")
n_methods <- length(methods)

tmp_results <-
  data.frame(matrix(nrow = n_sim,ncol = n_methods, dimnames = list(c(), methods)))
tmp_results_grp <-
  array(dim = c(n_sim, n_methods, 10),
        dimnames = list(c(), methods, c()))

mse_train_unif <- list()
mse_train_grp_unif <- list()
for(dgp_ix in 1:4){
  mse_train_unif[[dgp_ix]] <- 
    list(n1000 = tmp_results, 
         n5000 = tmp_results, n10000 = tmp_results)

  mse_train_grp_unif[[dgp_ix]] <- 
    list(n1000 = tmp_results_grp, 
         n5000 = tmp_results_grp, n10000 = tmp_results_grp)
}

mse_test_unif <- mse_train_unif
cov_train_unif <- mse_train_unif
cov_test_unif <- mse_train_unif
timing_unif <- mse_train_unif

mse_test_grp_unif <- mse_train_grp_unif
cov_train_grp_unif <- mse_train_grp_unif
cov_test_grp_unif <- mse_train_grp_unif

# Set up containers for imbalanced categories
mse_train_imbal <- mse_train_unif
mse_test_imbal <- mse_train_unif
cov_train_imbal <- mse_train_unif
cov_test_imbal <- mse_train_unif
timing_imbal <- mse_train_unif
mse_train_grp_imbal <- mse_train_grp_unif
mse_test_grp_imbal<- mse_train_grp_unif
cov_train_grp_imbal <- mse_train_grp_unif
cov_test_grp_imbal <- mse_train_grp_unif


for(job_id in 1:nrow(sim_settings)){
  
  file_name <- paste0("grouping_advantage_results_", job_id, ".RData")
  if(file.exists(file_name)){
    load(file_name)
    results <- get(paste0("grouping_advantage_", job_id))
    sim_number <- results$sim_number
    dgp_ix <- results$dgp_ix
    n_name <- paste0("n", results$n_all)
    
    if(results$cat_unif){
      mse_train_unif[[dgp_ix]][[n_name]][sim_number,] <- results$mse_train
      mse_test_unif[[dgp_ix]][[n_name]][sim_number,] <- results$mse_test
      cov_train_unif[[dgp_ix]][[n_name]][sim_number,] <- results$cov_train
      cov_test_unif[[dgp_ix]][[n_name]][sim_number,] <- results$cov_test
      timing_unif[[dgp_ix]][[n_name]][sim_number,] <- results$timing
      
      mse_train_grp_unif[[dgp_ix]][[n_name]][sim_number,,] <- results$mse_train_grp
      mse_test_grp_unif[[dgp_ix]][[n_name]][sim_number,,] <- results$mse_test_grp
      cov_train_grp_unif[[dgp_ix]][[n_name]][sim_number,,] <- results$cov_train_grp
      cov_test_grp_unif[[dgp_ix]][[n_name]][sim_number,,] <- results$cov_test_grp
    } else{
      mse_train_imbal[[dgp_ix]][[n_name]][sim_number,] <- results$mse_train
      mse_test_imbal[[dgp_ix]][[n_name]][sim_number,] <- results$mse_test
      cov_train_imbal[[dgp_ix]][[n_name]][sim_number,] <- results$cov_train
      cov_test_imbal[[dgp_ix]][[n_name]][sim_number,] <- results$cov_test
      timing_imbal[[dgp_ix]][[n_name]][sim_number,] <- results$timing
      
      mse_train_grp_imbal[[dgp_ix]][[n_name]][sim_number,,] <- results$mse_train_grp
      mse_test_grp_imbal[[dgp_ix]][[n_name]][sim_number,,] <- results$mse_test_grp
      cov_train_grp_imbal[[dgp_ix]][[n_name]][sim_number,,] <- results$cov_train_grp
      cov_test_grp_imbal[[dgp_ix]][[n_name]][sim_number,,] <- results$cov_test_grp
    }
    
    rm(list = paste0("grouping_advantage_", job_id))
    rm(results)
  }
}
save(mse_train_unif, mse_test_unif, cov_train_unif, cov_test_unif, timing_unif, 
     mse_train_grp_unif, mse_test_grp_unif, cov_train_grp_unif, cov_test_grp_unif,
     mse_train_imbal, mse_test_imbal, cov_train_imbal, cov_test_imbal, timing_imbal, 
     mse_train_grp_imbal, mse_test_grp_imbal, cov_train_grp_imbal, cov_test_grp_imbal,
     file = "grouping_advantage_results.RData")