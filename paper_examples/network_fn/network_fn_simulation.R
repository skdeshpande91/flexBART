source("helpers.R")
load("sim_settings.RData")
load("network_fn_data.RData")

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
  p_cont <- 10
  set.seed(618*n_obs + sim_number)
  nd <- 1000
  burn <- 1000
  ################################################################################
  # Generate the data and get training/testing split
  ################################################################################
  source("generate_network_fn_data.R")
  if(alg == "networkBART"){
    timing <- 
      system.time(
        fit <- 
          flexBART::network_BART(Y_train = Y_train,
                                 X_cont_train = X_cont_train,
                                 vertex_id_train = vertex_id_train,
                                 X_cont_test = X_cont_test,
                                 vertex_id_test = vertex_id_test,
                                 A = A,
                                 sparse = sparse,
                                 nd = nd, burn = burn,
                                 graph_split = TRUE,
                                 graph_cut_type = graph_cut_type,
                                 save_trees = FALSE, verbose = FALSE))[["elapsed"]]
  } else if(alg == "flexBART"){
    timing <- 
      system.time(
        fit <- 
          flexBART::flexBART(Y_train = Y_train,
                             X_cont_train = X_cont_train,
                             X_cat_train = X_cat_train,
                             X_cont_test = X_cont_test,
                             X_cat_test = X_cat_test,
                             cat_levels_list = list(0:n_vertex-1),
                             sparse = sparse, 
                             nd = nd, burn = burn,
                             save_trees = FALSE, verbose = FALSE))[["elapsed"]]
  } else if(alg == "BART"){
    if(is.na(embed)){
      timing <- system.time(
        fit <- BART::wbart(x.train = bart_df_train,
                      y.train = Y_train,
                      x.test = bart_df_test, 
                      sparse = sparse,
                      ndpost = nd, nskip = burn,
                      printevery = nd + burn + 1))[["elapsed"]]
    } else if(embed == "ase"){
      timing <- system.time(
        fit <-
          BART::wbart(x.train = ase_mm_train,
                      y.train = Y_train,
                      x.test = ase_mm_test, 
                      sparse = sparse,
                      ndpost = nd, nskip = burn,
                      printevery = nd + burn + 1))[["elapsed"]]
    }
  }
  
  train_sum <- summarize(mu_samples = fit$yhat.train)
  
  test_sum <- summarize(mu_samples = fit$yhat.test)
  test_sum1 <- test_sum[test_index_1,]
  test_sum2 <- test_sum[test_index_2,]

  
  mse_train <- mean( (mu_train - train_sum[,"MEAN"])^2 )
  mse_test1 <- mean( (mu_test1 - test_sum1[,"MEAN"])^2 )
  mse_test2 <- mean( (mu_test2 - test_sum2[,"MEAN"])^2 )
  
  cov_train <-
    coverage(mu_train,
             lower = train_sum[,"L95"],
             upper = train_sum[,"U95"])
  
  cov_test1 <-
    coverage(mu_test1,
             lower = test_sum1[,"L95"],
             upper = test_sum1[,"U95"])
  cov_test2 <-
    coverage(mu_test2,
             lower = test_sum2[,"L95"],
             upper = test_sum2[,"U95"])

  assign(paste0("network_fn_", job_id),
         list(alg = alg, graph_cut_type = graph_cut_type,
              sparse = sparse, embed = embed, embed_dim = embed_dim, 
              n_obs = n_obs,
              sim_number = sim_number,
              mse_train = mse_train,
              mse_test1 = mse_test1,
              mse_test2 = mse_test2,
              cov_train = cov_train,
              cov_test1 = cov_test1,
              cov_test2 = cov_test2,
              timing = timing))
  save(list = paste0("network_fn_", job_id),
       file = paste0("network_fn_results_", job_id, ".RData"))  
}