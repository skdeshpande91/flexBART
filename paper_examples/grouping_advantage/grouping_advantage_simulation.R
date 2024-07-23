source("base_functions.R")
source("true_mu.R")
source("bart_oracle.R")
source("bart_target.R")
source("helpers.R")
load("sim_settings.RData")

args <- commandArgs(TRUE)
block_id <- as.numeric(args[1])

n_list <- unique(sim_settings$n_all)
methods_list <- c("flexBART", "flexDART",
                  "BART", "DART",
                  "oracleBART", "oracleDART",
                  "targetBART", "targetDART")

for(job_id in block_starts[block_id]:block_ends[block_id]){
  
  n_all <- sim_settings[job_id, "n_all"]
  dgp_ix <- sim_settings[job_id, "dgp"]
  cat_unif <- sim_settings[job_id, "cat_unif"]
  sim_number <- sim_settings[job_id, "sim_number"]
  n_ix <- which(n_list == n_all)
  set.seed(620*n_ix + dgp_ix * 100 + sim_number)
  sigma <- 1
  source("generate_data.R")
  
  nd <- 1000
  burn <- 1000
  
  tmp_results <- rep(NA, times = length(methods_list))
  names(tmp_results) <- methods_list
  
  tmp_results_grp <- matrix(NA, nrow = length(methods_list), ncol = 10, 
                            dimnames = list(methods_list, c()))
  
  mse_train <- tmp_results
  mse_test <- tmp_results
  cov_train <- tmp_results
  cov_test <- tmp_results
  timing <- tmp_results
  mse_train_grp <- tmp_results_grp
  mse_test_grp <- tmp_results_grp
  cov_train_grp <- tmp_results_grp
  cov_test_grp <- tmp_results_grp
  
  
  flexBART_time <- 
    system.time(
      flexBART_fit <-
        flexBART::flexBART(Y_train = Y_train,
                           X_cont_train = X_cont_train,
                           X_cat_train = X_cat_train,
                           X_cont_test = X_cont_test,
                           X_cat_test = X_cat_test,
                           unif_cuts = unif_cuts,
                           cat_levels_list = cat_levels_list,
                           sparse = FALSE,
                           nd = nd, burn = burn, verbose = FALSE))["elapsed"]
  
  flexDART_time <- 
    system.time(
      flexDART_fit <-
        flexBART::flexBART(Y_train = Y_train,
                           X_cont_train = X_cont_train,
                           X_cat_train = X_cat_train,
                           X_cont_test = X_cont_test,
                           X_cat_test = X_cat_test,
                           unif_cuts = unif_cuts,
                           cat_levels_list = cat_levels_list,
                           sparse = TRUE,
                           nd = nd, burn = burn, verbose = FALSE))["elapsed"]
  
  BART_time <-
    system.time(
      BART_fit <-
        BART::wbart(x.train = bart_df_train,
                    y.train = Y_train,
                    x.test = bart_df_test,
                    ndpost = nd, nskip = burn,
                    sparse = FALSE,
                    sigest = sd(Y_train),
                    printevery = 2001))[[1]]
  
  DART_time <-
    system.time(
      DART_fit <-
        BART::wbart(x.train = bart_df_train,
                    y.train = Y_train,
                    x.test = bart_df_test,
                    ndpost = nd, nskip = burn,
                    sparse = TRUE,
                    sigest = sd(Y_train),
                    printevery = 2001))[[1]]
  
  oracleBART_time <-
    system.time(
      oracleBART_fit <-
        bart_oracle(Y_train = Y_train,
                    df_train = bart_df_train,
                    df_test = bart_df_test,
                    ix_list = oracle_ix_list,
                    sparse = FALSE))[[1]]
  oracleDART_time <-
    system.time(
      oracleDART_fit <-
        bart_oracle(Y_train = Y_train,
                    df_train = bart_df_train,
                    df_test = bart_df_test,
                    ix_list = oracle_ix_list,
                    sparse = TRUE))[[1]]
  
  targetBART_time <-
    system.time(
      targetBART_fit <-
        bart_target(Y_train = Y_train,
                    df_train = bart_df_train,
                    df_test = bart_df_test,
                    sparse = FALSE))["elapsed"]
  targetDART_time <-
    system.time(
      targetDART_fit <-
        bart_target(Y_train = Y_train,
                    df_train = bart_df_train,
                    df_test = bart_df_test,
                    sparse = TRUE))["elapsed"]
  
  for(m in methods_list){
    fit <- get(paste0(m, "_fit"))
    timing[m] <- get(paste0(m, "_time"))
    train_sum <- summarize(fit$yhat.train)
    test_sum <- summarize(fit$yhat.test)
    
    mse_train[m] <- 
      square_error(mu_train, train_sum[,"MEAN"])
    cov_train[m] <-
      coverage(mu_train, train_sum[,"L95"], train_sum[,"U95"])
    mse_test[m] <- 
      square_error(mu_test, test_sum[,"MEAN"])
    cov_test[m] <-
      coverage(mu_test, test_sum[,"L95"], test_sum[,"U95"])

    for(l in 0:9){
      train_ix <- which(X_cat_train[,1] == l)
      test_ix <- which(X_cat_test[,1] == l)
      
      mse_train_grp[m,l+1] <-
        square_error(mu_train[train_ix], train_sum[train_ix,"MEAN"])
      cov_train_grp[m,l+1] <-
        coverage(mu_train[train_ix], train_sum[train_ix,"L95"], train_sum[train_ix,"U95"])
      mse_test_grp[m,l+1] <-
        square_error(mu_test[test_ix], test_sum[test_ix,"MEAN"])
      cov_test_grp[m,l+1] <-
        coverage(mu_test[test_ix], test_sum[test_ix,"L95"], test_sum[test_ix,"U95"])
    }
  }
  assign(paste0("grouping_advantage_", job_id),
         list(job_id = job_id,
              n_all = n_all,
              sim_number = sim_number,
              dgp_ix = dgp_ix,
              cat_unif = cat_unif,
              mse_train = mse_train,
              mse_test = mse_test,
              cov_train = cov_train,
              cov_test = cov_test,
              mse_train_grp = mse_train_grp,
              mse_test_grp = mse_test_grp,
              cov_train_grp = cov_train_grp,
              cov_test_grp = cov_test_grp,
              timing = timing))
  save(list = paste0("grouping_advantage_", job_id),
       file = paste0("grouping_advantage_results_", job_id, ".RData"))
}
