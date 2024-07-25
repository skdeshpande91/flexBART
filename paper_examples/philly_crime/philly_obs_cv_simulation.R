load("philly_crime_data.RData")
load("sim_settings.RData")
source("helpers.R")
source("targetBART.R")

args <- commandArgs(TRUE)
block_id <- as.numeric(args[1])

for(job_id in block_starts[block_id]:block_ends[block_id]){
  method <- sim_settings[job_id, "method"]
  sim_number <- sim_settings[job_id, "sim_number"]
  

  set.seed(405+sim_number)
  test_index <- sort(sample(1:N, size = floor(0.25 * N), replace = FALSE))
  train_index <- (1:N)[-test_index]
  
  
  X_cont_train <- matrix(X_cont_all[train_index,], ncol = 1)
  X_cont_test <- matrix(X_cont_all[test_index,], ncol = 1)
  
  X_cat_train <- matrix(X_cat_all[train_index,], ncol = 1)
  X_cat_test <- matrix(X_cat_all[test_index,], ncol = 1)
  
  vertex_id_train <- vertex_id_all[train_index]
  vertex_id_test <- vertex_id_all[test_index]
  Y_train <- Y_all[train_index]
  Y_test <- Y_all[test_index]
  
  bart_df_train <- bart_df_all[train_index,]
  bart_df_test <- bart_df_all[test_index,]
  
  set.seed(724+sim_number)
  
  nd <- 1000
  burn <- 1000
  
  
  if(method %in% paste0("networkBART_", 1:7)){
    graph_cut_type <- as.numeric(strsplit(method, split = "_")[[1]][2])
    fit_time <- system.time(
      fit <- 
        flexBART::network_BART(Y_train = Y_train,
                               vertex_id_train = vertex_id_train,
                               X_cont_train = X_cont_train,
                               vertex_id_test = vertex_id_test,
                               X_cont_test = X_cont_test,
                               unif_cuts = unif_cuts,
                               cutpoints_list = cutpoints_list,
                               A = A_tract, 
                               graph_cut_type = graph_cut_type,
                               save_trees = FALSE,
                               nd = nd, burn = burn,
                               verbose = FALSE))
  } else if(method == "flexBART"){
    fit_time <- system.time(
      fit <- 
        flexBART::flexBART(Y_train = Y_train,
                           X_cont_train = X_cont_train,
                           X_cat_train = X_cat_train,
                           X_cont_test = X_cont_test,
                           X_cat_test = X_cat_test,
                           unif_cuts = unif_cuts,
                           cutpoints_list = cutpoints_list,
                           cat_levels_list = cat_levels_list,
                           nd = nd, burn = burn,
                           save_trees = FALSE,
                           verbose = FALSE))
  } else if(method == "BART"){
    fit_time <- system.time(
      fit <- 
        BART::wbart(x.train = bart_df_train,
                    y.train = Y_train, x.test = bart_df_test, 
                    sparse = FALSE, sigest = sd(Y_train),
                    ndpost = nd, nskip = burn, printevery = nd + burn + 1))
  } else if(method == "DART"){
    fit_time <- system.time(
      fit <- 
        BART::wbart(x.train = bart_df_train,
                    y.train = Y_train, x.test = bart_df_test, 
                    sparse = TRUE, sigest = sd(Y_train),
                    ndpost = nd, nskip = burn, printevery = nd + burn + 1))
  } else if(method %in% c("BART_ase1", "BART_ase3", "BART_ase5")){
    embed_dim <- as.numeric(strsplit(method, split = "_ase")[[1]][2])
    ase_df_all <- data.frame(Time = X_cont_all, X_ase_all[,1:embed_dim])
    ase_df_train <- ase_df_all[train_index,]
    ase_df_test <- ase_df_all[test_index,]
    fit_time <- system.time(
      fit <-
        BART::wbart(x.train = ase_df_train,
                    y.train = Y_train,
                    x.test = ase_df_test,
                    sparse = FALSE, sigest = sd(Y_train),
                    ndpost = nd, nskip = burn, printevery = nd+burn+1))["elapsed"]
  } else if (method == "BART_latlon"){
    latlon_df_train <- bart_latlon_df_all[train_index,]
    latlon_df_test <- bart_latlon_df_all[test_index,]
    fit_time <- system.time(
      fit <-
        BART::wbart(x.train = latlon_df_train,
                    y.train = Y_train,
                    x.test = latlon_df_test,
                    sparse = FALSE, sigest = sd(Y_train),
                    ndpost = nd, nskip = burn, printevery = nd+burn+1))["elapsed"]
  } else if(method == "targetBART"){
    fit_time <-
      system.time(
        fit <-
          targetBART(Y_train, 
                     df_train = bart_df_train,
                     df_test = bart_df_test,
                     p_cont = 1, p_cat = 1,
                     cat_levels_list = cat_levels_list,
                     nd = nd, burn = burn, sparse = FALSE))["elapsed"]
  } else{
    stop("Bad method")
  }
  
  sum_train <- summarize_post_pred(fit$yhat.train, fit$sigma[-(1:burn)])
  sum_test <- summarize_post_pred(fit$yhat.test, fit$sigma[-(1:burn)])
  
  assign(paste0("philly_obs_cv_", job_id),
         list(method = method, sim_number = sim_number,
              rmse_train = sqrt(mean( (Y_train - sum_train[,"MEAN"])^2)),
              cov_train = coverage(Y_train, sum_train[,"L95"], sum_train[,"U95"]),
              rmse_test = sqrt(mean( (Y_test - sum_test[,"MEAN"])^2)),
              cov_test = coverage(Y_test, sum_test[,"L95"], sum_test[,"U95"]),
              timing = fit_time))
  save(list = paste0("philly_obs_cv_", job_id),
       file  = paste0("philly_obs_cv_results", job_id, ".RData"))
}




