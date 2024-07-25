load("sim_settings.RData")
source("helpers.R")
source("targetBART.R")

args <- commandArgs(TRUE)
block_id <- as.numeric(args[1])

for(job_id in block_starts[block_id]:block_ends[block_id]){
  dataset_name <- sim_settings[job_id, "dataset"]
  data_ix <- which(data_list == dataset_name)
  
  sim_number <- sim_settings[job_id, "sim_number"]
  load(paste0("data/", dataset_name, ".RData"))
  tmp_data <- get(paste0(dataset_name, "_data"))
  n_all <- tmp_data$n_all
  p_cont <- tmp_data$p_cont
  p_cat <- tmp_data$p_cat
  X_cont_all <- tmp_data$X_cont_all
  X_cat_all <- tmp_data$X_cat_all
  bart_df_all <- tmp_data$bart_df_all
  bart_mm_all <- tmp_data$bart_mm_all
  Y_all <- tmp_data$Y_all
  
  unif_cuts <- tmp_data$unif_cuts
  cutpoints_list <- tmp_data$cutpoints_list
  cat_levels_list <- tmp_data$cat_levels_list
  rm(tmp_data)
  
  ############################
  # Create test/train split
  ############################
  set.seed(619*data_ix + sim_number)
  train_index <- sort(sample(1:n_all, floor(0.75*n_all), replace = FALSE))
  test_index <- sort( (1:n_all)[-train_index])
  
  X_cont_train <- X_cont_all[train_index,]
  X_cont_test <- X_cont_all[test_index,]

  if(p_cat > 1){
    X_cat_train <- X_cat_all[train_index,]
    X_cat_test <- X_cat_all[test_index,]
  } else{
    X_cat_train <- matrix(X_cat_all[train_index,], ncol = 1)
    X_cat_test <- matrix(X_cat_all[test_index,], ncol = 1)
  }
  bart_mm_train <- bart_mm_all[train_index,]
  bart_mm_test <- bart_mm_all[test_index,]
  
  bart_df_train <- bart_df_all[train_index,]
  bart_df_test <- bart_df_all[test_index,]
  
  
  Y_train <- Y_all[train_index]
  Y_test <- Y_all[test_index]
  nd <- 1000
  burn <- 1000
  
  # Output containers
  # Get the results
  smse_train <- rep(NA, times = 6)
  names(smse_train) <- c("flexBART", "flexDART",
                         "BART1","DART1",
                         "targetBART", "targetDART")
  smse_test <- smse_train
  cov_train <- smse_train
  cov_test <- smse_train
  timing <- smse_train
  
  flexbart_time <- 
    system.time(
      flexbart_fit <-
        flexBART::flexBART(Y_train = Y_train,
                           X_cont_train = X_cont_train,
                           X_cat_train = X_cat_train,
                           X_cont_test = X_cont_test,
                           X_cat_test = X_cat_test,
                           unif_cuts = unif_cuts,
                           cutpoints_list = cutpoints_list,
                           cat_levels_list = cat_levels_list,
                           sparse = FALSE,
                           verbose = FALSE))["elapsed"]
  flexbart_train_sum <- summarize(flexbart_fit$yhat.train,
                                  flexbart_fit$sigma[-(1:burn)])
  flexbart_test_sum <- summarize(flexbart_fit$yhat.test,
                                 flexbart_fit$sigma[-(1:burn)])
  
  smse_train["flexBART"] <- 
    smse(truth = Y_train, pred = flexbart_train_sum[,"MEAN"], baseline = mean(Y_train))
  smse_test["flexBART"] <- 
    smse(truth = Y_test, pred = flexbart_test_sum[,"MEAN"], baseline = mean(Y_train))
  
  cov_train["flexBART"] <-
    coverage(truth = Y_train, lower = flexbart_train_sum[,"L95"], upper = flexbart_train_sum[,"U95"])
  cov_test["flexBART"] <-
    coverage(truth = Y_test, lower = flexbart_test_sum[,"L95"], upper = flexbart_test_sum[,"U95"])
  
  timing["flexBART"] <- flexbart_time
  
  flexdart_time <- 
    system.time(
      flexdart_fit <-
        flexBART::flexBART(Y_train = Y_train,
                           X_cont_train = X_cont_train,
                           X_cat_train = X_cat_train,
                           X_cont_test = X_cont_test,
                           X_cat_test = X_cat_test,
                           unif_cuts = unif_cuts,
                           cutpoints_list = cutpoints_list,
                           cat_levels_list = cat_levels_list,
                           sparse = TRUE,
                           verbose = FALSE))["elapsed"]
  flexdart_train_sum <- summarize(flexdart_fit$yhat.train,
                                  flexdart_fit$sigma[-(1:burn)])
  flexdart_test_sum <- summarize(flexdart_fit$yhat.test,
                                 flexdart_fit$sigma[-(1:burn)])
  
  smse_train["flexDART"] <- 
    smse(truth = Y_train, pred = flexdart_train_sum[,"MEAN"], baseline = mean(Y_train))
  smse_test["flexDART"] <- 
    smse(truth = Y_test, pred = flexdart_test_sum[,"MEAN"], baseline = mean(Y_train))
  cov_train["flexDART"] <-
    coverage(truth = Y_train, lower = flexdart_train_sum[,"L95"], upper = flexdart_train_sum[,"U95"])
  cov_test["flexDART"] <-
    coverage(truth = Y_test, lower = flexdart_test_sum[,"L95"], upper = flexdart_test_sum[,"U95"])
  
  timing["flexDART"] <- flexdart_time
  
  bart1_time <-
    system.time(
      bart1_fit <-
        BART::wbart(x.train = bart_mm_train,
                    y.train = Y_train,
                    x.test = bart_mm_test,
                    rm.const = TRUE,
                    sparse = FALSE,
                    sigest = sd(Y_train),
                    ndpost = nd, nskip = burn,
                    printevery = 1+burn+nd))["elapsed"]
  
  bart1_train_sum <- summarize(bart1_fit$yhat.train,
                               bart1_fit$sigma[-(1:burn)])
  bart1_test_sum <- summarize(bart1_fit$yhat.test,
                              bart1_fit$sigma[-(1:burn)])
  smse_train["BART1"] <- 
    smse(truth = Y_train, pred = bart1_train_sum[,"MEAN"], baseline = mean(Y_train))
  smse_test["BART1"] <- 
    smse(truth = Y_test, pred = bart1_test_sum[,"MEAN"], baseline = mean(Y_train))
  cov_train["BART1"] <-
    coverage(truth = Y_train, lower = bart1_train_sum[,"L95"], upper = bart1_train_sum[,"U95"])
  cov_test["flexBART"] <-
    coverage(truth = Y_test, lower = bart1_test_sum[,"L95"], upper = bart1_test_sum[,"U95"])
  timing["BART1"] <- bart1_time
  

  dart1_time <-
    system.time(
      dart1_fit <-
        BART::wbart(x.train = bart_mm_train,
                    y.train = Y_train,
                    x.test = bart_mm_test,
                    rm.const = TRUE,
                    sparse = TRUE,
                    sigest = sd(Y_train),
                    ndpost = nd, nskip = burn,
                    printevery = 1 + nd + burn))["elapsed"]
  dart1_train_sum <- summarize(dart1_fit$yhat.train,
                               dart1_fit$sigma[-(1:burn)])
  dart1_test_sum <- summarize(dart1_fit$yhat.test,
                              dart1_fit$sigma[-(1:burn)])
  smse_train["DART1"] <- 
    smse(truth = Y_train, pred = dart1_train_sum[,"MEAN"], baseline = mean(Y_train))
  smse_test["DART1"] <- 
    smse(truth = Y_test, pred = dart1_test_sum[,"MEAN"], baseline = mean(Y_train))
  cov_train["DART1"] <-
    coverage(truth = Y_train, lower = dart1_train_sum[,"L95"], upper = dart1_train_sum[,"U95"])
  cov_test["DART1"] <-
    coverage(truth = Y_test, lower = dart1_test_sum[,"L95"], upper = dart1_test_sum[,"U95"])
  timing["DART1"] <- dart1_time

  
  target_bart_time <-
    system.time(
      target_bart_fit <-
        targetBART(Y_train, 
                   df_train = bart_df_train,
                   df_test = bart_df_test,
                   p_cont = p_cont, p_cat = p_cat,
                   cat_levels_list = cat_levels_list,
                   nd = nd, burn = burn, sparse = FALSE))["elapsed"]
  target_bart_train_sum <- summarize(target_bart_fit$yhat.train,
                                target_bart_fit$sigma[-(1:burn)])
  target_bart_test_sum <- summarize(target_bart_fit$yhat.test,
                               target_bart_fit$sigma[-(1:burn)])
  smse_train["targetBART"] <- 
    smse(truth = Y_train, pred = target_bart_train_sum[,"MEAN"], baseline = mean(Y_train))
  smse_test["targetBART"] <- 
    smse(truth = Y_test, pred = target_bart_test_sum[,"MEAN"], baseline = mean(Y_train))
  cov_train["targetBART"] <-
    coverage(truth = Y_train, lower = target_bart_train_sum[,"L95"], upper = target_bart_train_sum[,"U95"])
  cov_test["targetBART"] <-
    coverage(truth = Y_test, lower = target_bart_test_sum[,"L95"], upper = target_bart_test_sum[,"U95"])
  timing["targetBART"] <- target_bart_time
  
  target_dart_time <-
    system.time(
      target_dart_fit <-
        targetBART(Y_train, 
                   df_train = bart_df_train,
                   df_test = bart_df_test,
                   p_cont = p_cont, p_cat = p_cat,
                   cat_levels_list = cat_levels_list,
                   nd = nd, burn = burn, sparse = TRUE))["elapsed"]
  target_dart_train_sum <- summarize(target_dart_fit$yhat.train,
                                target_dart_fit$sigma[-(1:burn)])
  target_dart_test_sum <- summarize(target_dart_fit$yhat.test,
                               target_dart_fit$sigma[-(1:burn)])
  smse_train["targetDART"] <- 
    smse(truth = Y_train, pred = target_dart_train_sum[,"MEAN"], baseline = mean(Y_train))
  smse_test["targetDART"] <- 
    smse(truth = Y_test, pred = target_dart_test_sum[,"MEAN"], baseline = mean(Y_train))
  cov_train["targetDART"] <-
    coverage(truth = Y_train, lower = target_dart_train_sum[,"L95"], upper = target_dart_train_sum[,"U95"])
  cov_test["targetDART"] <-
    coverage(truth = Y_test, lower = target_dart_test_sum[,"L95"], upper = target_dart_test_sum[,"U95"])
  timing["targetDART"] <- target_dart_time
  
  
  assign(paste0("bakeoff_", job_id),
         list(dataset_name = dataset_name,
              sim_number = sim_number,
              smse_train = smse_train,
              smse_test = smse_test,
              cov_train = cov_train,
              cov_test = cov_test,
              timing = timing))
  save(list = paste0("bakeoff_", job_id),
       file = paste0("bakeoff_results_", job_id, ".RData"))
}