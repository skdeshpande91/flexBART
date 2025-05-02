library(Rcpp)
library(RcppArmadillo)
sourceCpp("../flexBART/src/detect_nesting.cpp")
source("../flexBART/R/get_z_col_id.R")
source("../flexBART/R/flexBART_structures.R")
source("../flexBART/R/get_covariate_info.R")
source("../flexBART/R/parse_formula.R")
source("../flexBART/R/prepare_data.R")
source("../flexBART/R/parse_controls.R")
source("../flexBART/R/parse_hyper.R")
source("../flexBART/R/get_sigma.R")

sourceCpp("../flexBART/src/single_ensm_unnested_fit.cpp")
sourceCpp("../flexBART/src/multiple_ensemble_unnested_fit.cpp")


rmse_train <- c(single = NA, multi = NA, dbarts = NA, bart = NA, flex = NA)
rmse_test <- c(single = NA, multi = NA, dbarts = NA, bart = NA, flex = NA)
timing <- c(single = NA, multi = NA, dbarts = NA, bart = NA, flex = NA)



source("generate_friedman_data.R")

train_data <- friedman_train
test_data <- friedman_test

frmla <- Y~bart(.)
tmp_form <- parse_formula(frmla, train_data)
outcome_name <- tmp_form$outcome_name
cov_ensm <- tmp_form$cov_ensm

tmp_data <- prepare_data(train_data = train_data,
                         outcome_name = outcome_name,
                         cov_ensm = cov_ensm,
                         test_data = test_data)

R <- ncol(cov_ensm)
y_range <- 
  max(tmp_data$training_info$std_Y) - min(tmp_data$training_info$std_Y)
sigest <- 
  get_sigma(tmp_data$training_info)
hyper <- 
  parse_hyper(R = R,
              y_range = y_range,
              sigest = sigest, M = 200, sparse = FALSE)


control <- parse_controls()

y_mean <- tmp_data$training_info$y_mean
y_sd <- tmp_data$training_info$y_sd


###############################
# New single ensemble fit

single_time <-
  system.time(
    single_fit <- 
      .single_unnested_fit(Y_train = tmp_data$training_info$std_Y,
                           tZ_train = t(tmp_data$training_info$Z),
                           tX_cont_train = t(tmp_data$training_info$X_cont),
                           tX_cat_train = t(tmp_data$training_info$X_cat),
                           tZ_test = t(tmp_data$testing_info$Z),
                           tX_cont_test = t(tmp_data$testing_info$X_cont),
                           tX_cat_test = t(tmp_data$testing_info$X_cat),
                           cutpoints_list = tmp_data$training_info$cutpoints,
                           cat_levels_list = tmp_data$training_info$cat_levels_list,
                           edge_mat_list = tmp_data$training_info$edge_mat_list,
                           graph_cut_type = hyper$graph_cut_type,
                           sparse = hyper$sparse, a_u = hyper$a_u, b_u = hyper$b_u,
                           mu0 = hyper$mu0_vec[1], tau = hyper$tau_vec[1],
                           sigest = hyper$sigest, lambda = hyper$lambda, nu = hyper$nu,
                           M = hyper$M_vec[1],
                           nd = control$nd, burn = control$burn, thin = control$thin,
                           save_samples = control$save_samples, save_trees = control$save_trees,
                           verbose = control$verbose, print_every = control$print_every))

muhat_single_train <- y_mean + y_sd * single_fit$fit_train_mean
muhat_single_test <- y_mean + y_sd * single_fit$fit_test_mean

rmse_train["single"] <- sqrt( mean( (mu_train - muhat_single_train)^2 ))
rmse_test["single"] <- sqrt( mean( (mu_test - muhat_single_test)^2 ))
timing["single"] <- single_time["elapsed"]


multi_time <-
  system.time(
    multi_fit <- 
      ._multi_unnested_fit(Y_train = tmp_data$training_info$std_Y,
                           tZ_train = t(tmp_data$training_info$Z),
                           tX_cont_train = t(tmp_data$training_info$X_cont),
                           tX_cat_train = t(tmp_data$training_info$X_cat),
                           sigest = hyper$sigest,
                           cov_ensm = cov_ensm,
                           cutpoints_list = tmp_data$training_info$cutpoints,
                           cat_levels_list = tmp_data$training_info$cat_levels_list,
                           edge_mat_list = tmp_data$training_info$edge_mat_list,
                           tZ_test = t(tmp_data$testing_info$Z),
                           tX_cont_test = t(tmp_data$testing_info$X_cont),
                           tX_cat_test = t(tmp_data$testing_info$X_cat),
                           M_vec = hyper$M_vec,
                           alpha_vec = hyper$alpha_vec, beta_vec = hyper$beta_vec,
                           mu0_vec = hyper$mu0_vec, tau_vec = hyper$tau_vec,
                           graph_cut_type = hyper$graph_cut_type,
                           sparse = hyper$sparse, a_u = hyper$a_u, b_u = hyper$b_u,
                           nu = hyper$nu,lambda = hyper$lambda, 
                           nd = control$nd, burn = control$burn, thin = control$thin,
                           save_samples = control$save_samples, save_trees = control$save_trees,
                           verbose = control$verbose, print_every = control$print_every))

muhat_multi_train <- y_mean + y_sd * multi_fit$fit_train_mean
muhat_multi_test <- y_mean + y_sd * multi_fit$fit_test_mean


rmse_train["multi"] <- sqrt( mean( (mu_train - muhat_multi_train)^2 ))
rmse_test["multi"] <- sqrt( mean( (mu_test - muhat_multi_test)^2 ))
timing["multi"] <- multi_time["elapsed"]

dbart_time <-
  system.time(
    dbart_fit <- dbarts::bart(x.train = friedman_train[,colnames(friedman_train) != "Y"],
                              y.train = friedman_train[,"Y"],
                              x.test = friedman_test[,colnames(friedman_test) != "Y"],
                              numcut = 100,
                              ndpost = 1000, nskip = 1000, keeptrees = TRUE))
rmse_train["dbarts"] <- sqrt( mean( (mu_train - dbart_fit$yhat.train.mean)^2 ))
rmse_test["dbarts"] <- sqrt( mean( (mu_test - dbart_fit$yhat.test.mean)^2 ))
timing["dbarts"] <- dbart_time["elapsed"]

bart_time <-
  system.time(
    bart_fit <- BART::wbart(x.train = friedman_train[,colnames(friedman_train) != "Y"],
                            y.train = friedman_train[,"Y"],
                            x.test = friedman_test[,colnames(friedman_test) != "Y"],
                            ndpost = 1000, nskip = 1000))

rmse_train["bart"] <- sqrt( mean( (mu_train - bart_fit$yhat.train.mean)^2 ))
rmse_test["bart"] <- sqrt( mean( (mu_test - bart_fit$yhat.test.mean)^2 ))
timing["bart"] <- bart_time["elapsed"]

flex_time <-
  system.time(
    flex_fit <- flexBART::flexBART(Y_train =friedman_train[,"Y"],
                                   X_cont_train = tmp_data$training_info$X_cont,
                                   X_cat_train = tmp_data$training_info$X_cat,
                                   X_cont_test = tmp_data$testing_info$X_cont,
                                   X_cat_test = tmp_data$testing_info$X_cat,
                                   cat_levels_list = tmp_data$training_info$cat_levels_list,
                                   sparse = hyper$sparse))

rmse_train["flex"] <- sqrt( mean( (mu_train - flex_fit$yhat.train.mean)^2 ))
rmse_test["flex"] <- sqrt( mean( (mu_test - flex_fit$yhat.test.mean)^2 ))
timing["flex"] <- flex_time["elapsed"]    
