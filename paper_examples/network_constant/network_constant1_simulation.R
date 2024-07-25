source("helpers.R")
load("sim_settings.RData")
load("network_constant1_data.RData")

args <- commandArgs(TRUE)
block_id <- as.numeric(args[1])

for(job_id in block_starts[block_id]:block_ends[block_id]){
  alg <- sim_settings[job_id, "alg"]
  graph_cut_type <- sim_settings[job_id, "graph_cut_type"]
  sparse <- sim_settings[job_id, "sparse"]
  embed <- sim_settings[job_id, "embed"]
  embed_dim <- sim_settings[job_id, "embed_dim"]
  sim_number <- sim_settings[job_id, "sim_number"]
  n_obs <- sim_settings[job_id, "n_obs"]
  
  sigma <- 1
  set.seed(618*n_obs + sim_number)
  nd <- 1000
  burn <- 1000
  ################################################################################
  # Generate the data and get training/testing split
  ################################################################################
  source("generate_network_constant1_data.R")
  if(alg == "networkBART"){
    timing <- 
      system.time(
        fit <- 
          flexBART::network_BART(Y_train = Y_train,
                                 vertex_id_train = vertex_id_train,
                                 vertex_id_test = vertex_id_test,
                                 A = A,
                                 sparse = sparse,
                                 nd = nd, burn = burn,
                                 graph_split = TRUE,
                                 graph_cut_type = graph_cut_type,
                                 save_trees = FALSE, verbose = TRUE))[["elapsed"]]
  } else if(alg == "flexBART"){
    timing <- 
      system.time(
        fit <- 
          flexBART::flexBART(Y_train = Y_train,
                             X_cat_train = X_cat_train,
                             X_cat_test = X_cat_test,
                             cat_levels_list = list(0:n_vertex-1),
                             sparse = sparse, 
                             nd = nd, burn = burn,
                             save_trees = FALSE, verbose = FALSE))[["elapsed"]]
  } else if(alg == "BART"){
    if(is.na(embed)){
      timing <- system.time(
        fit <-
          BART::wbart(x.train = bart_df_train,
                      y.train = Y_train,
                      x.test = bart_df_test, 
                      sparse = sparse,
                      sigest = sd(Y_train),
                      ndpost = nd, nskip = burn,
                      printevery = nd + burn + 1))[["elapsed"]]
    } else if(embed == "ase"){
      timing <- system.time(
        fit <-
          BART::wbart(x.train = ase_mm_train,
                      y.train = Y_train,
                      x.test = ase_mm_test, 
                      sparse = sparse,
                      sigest = sd(Y_train),
                      ndpost = nd, nskip = burn,
                      printevery = nd + burn + 1))[["elapsed"]]
    }
  }
  fit_sum <- summarize(mu_samples = fit$yhat.test)
  fit_sum_train <- fit_sum[train_vertices,]
  fit_sum_test <- fit_sum[test_vertices,]
  
  mu_train <- mu[train_vertices]
  mu_test <- mu[test_vertices]
  
  fit_mse_train <- mean( (mu_train - fit_sum_train[,"MEAN"])^2 )
  fit_mse_test <- mean( (mu[test_vertices] - fit_sum_test[,"MEAN"])^2 )
  
  fit_cov_train <- 
    coverage(mu_train,
             lower = fit_sum_train[,"L95"],
             upper = fit_sum_train[,"U95"])
  fit_cov_test <- 
    coverage(mu_test,
             lower = fit_sum_test[,"L95"],
             upper = fit_sum_test[,"U95"])
  assign(paste0("network_constant1_", job_id),
         list(alg = alg, graph_cut_type = graph_cut_type,
              sparse = sparse, embed = embed, embed_dim = embed_dim, 
              n_obs = n_obs,
              sim_number = sim_number,
              fit_mse_train = fit_mse_train,
              fit_mse_test = fit_mse_test,
              fit_cov_train = fit_cov_train,
              fit_cov_test = fit_cov_test,
              timing = timing))
  save(list = paste0("network_constant1_", job_id),
       file = paste0("network_constant1_results_", job_id, ".RData"))  
}