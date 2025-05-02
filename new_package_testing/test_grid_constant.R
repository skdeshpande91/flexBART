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

sourceCpp("../flexBART/src/single_ensm_fit.cpp")
sourceCpp("../flexBART/src/multi_ensm_fit.cpp")

source("generate_grid_data.R")

frmla <- Y ~ bart(X1)
tmp_form <- parse_formula(frmla, train_data)

outcome_name <- tmp_form$outcome_name
cov_ensm <- tmp_form$cov_ensm
adjacency_list <- list(X1 = A)


tmp_data <- 
  prepare_data(train_data = train_data,
               outcome_name = outcome_name, 
               cov_ensm = cov_ensm, 
               test_data = test_data,
               adjacency_list = adjacency_list)


R <- ncol(tmp_data$training_info$Z)
y_range <- 
  max(tmp_data$training_info$std_Y) - min(tmp_data$training_info$std_Y)

sigest <- 
  get_sigma(tmp_data$training_info, 
            tmp_data$data_info)
hyper <- 
  parse_hyper(R = R,
              y_range = y_range,
              sigest = sigest, M = 200, sparse = FALSE)

control <- parse_controls(save_samples = FALSE, save_trees = FALSE)


y_mean <- tmp_data$training_info$y_mean
y_sd <- tmp_data$training_info$y_sd

rmse_train <- c(single = NA, multi = NA, dbarts = NA, bart = NA)
rmse_test <- c(single = NA, multi = NA, dbarts = NA, bart = NA)
timing <- c(single = NA, multi = NA, dbarts = NA, bart = NA)

###############################
# New single ensemble fit

single_time <-
  system.time(
    single_fit <- 
      .single_fit(Y_train = tmp_data$training_info$std_Y,
                  cov_ensm = cov_ensm,
                  tX_cont_train = t(tmp_data$training_info$X_cont),
                  tX_cat_train = t(tmp_data$training_info$X_cat),
                  tX_cont_test = t(tmp_data$testing_info$X_cont),
                  tX_cat_test = t(tmp_data$testing_info$X_cat),
                  cutpoints_list = tmp_data$training_info$cutpoints,
                  cat_levels_list = tmp_data$training_info$cat_levels_list,
                  edge_mat_list = tmp_data$training_info$edge_mat_list,
                  nest_list = tmp_data$training_info$nest_list,
                  graph_cut_type = hyper$graph_cut_type,
                  sparse = hyper$sparse, 
                  a_u = hyper$a_u, b_u = hyper$b_u,
                  nest_v = hyper$nest_v,
                  nest_v_option = hyper$nest_v_option,
                  nest_c = hyper$nest_c,
                  mu0 = hyper$mu0_vec[1],
                  tau = hyper$tau_vec[1],
                  M = hyper$M_vec[1],
                  alpha = hyper$alpha_vec[1],
                  beta = hyper$beta_vec[1],
                  sigest = hyper$sigest,
                  nu = hyper$nu,
                  lambda = hyper$lambda,
                  nd = control$nd, burn = control$burn, thin = control$thin,
                  save_samples = control$save_samples, save_trees = control$save_trees,
                  verbose = control$verbose, print_every = control$print_every))

muhat_single <- y_mean + y_sd * single_fit$fit_test_mean

rmse_train["single"] <- sqrt( mean( (mu[train_vertices] - muhat_single[train_vertices])^2 ))
rmse_test["single"] <- sqrt( mean( (mu[test_vertices] - muhat_single[test_vertices])^2 ))
timing["single"] <- single_time["elapsed"]

multi_time <-
  system.time(
    multi_fit <- 
      ._multi_fit(Y_train = tmp_data$training_info$std_Y,
                  cov_ensm = cov_ensm,
                  tZ_train = t(tmp_data$training_info$Z),
                  tX_cont_train = t(tmp_data$training_info$X_cont),
                  tX_cat_train = t(tmp_data$training_info$X_cat),
                  tZ_test = t(tmp_data$testing_info$Z),
                  tX_cont_test = t(tmp_data$testing_info$X_cont),
                  tX_cat_test = t(tmp_data$testing_info$X_cat),
                  cutpoints_list = tmp_data$training_info$cutpoints,
                  cat_levels_list = tmp_data$training_info$cat_levels_list,
                  edge_mat_list = tmp_data$training_info$edge_mat_list,
                  nest_list = tmp_data$training_info$nest_list,
                  graph_cut_type = hyper$graph_cut_type,
                  sparse = hyper$sparse, 
                  a_u = hyper$a_u, b_u = hyper$b_u,
                  nest_v = hyper$nest_v,
                  nest_v_option = hyper$nest_v_option,
                  nest_c = hyper$nest_c,
                  M_vec = hyper$M_vec,
                  alpha_vec = hyper$alpha_vec, beta_vec = hyper$beta_vec,
                  mu0_vec = hyper$mu0_vec, tau_vec = hyper$tau_vec,
                  sigest = hyper$sigest,
                  nu = hyper$nu,lambda = hyper$lambda, 
                  nd = control$nd, burn = control$burn, thin = control$thin,
                  save_samples = control$save_samples, save_trees = control$save_trees,
                  verbose = control$verbose, print_every = control$print_every))

muhat_multi <- y_mean + y_sd * multi_fit$fit_test_mean

rmse_train["multi"] <- 
  sqrt( mean( (mu[train_vertices] - muhat_multi[train_vertices])^2 ))
rmse_test["multi"] <- 
  sqrt( mean( (mu[test_vertices] - muhat_multi[test_vertices])^2 ))
timing["multi"] <- multi_time["elapsed"]


dbart_time <-
  system.time(
    dbart_fit <- dbarts::bart(x.train = train_model_mat,
                              y.train = train_data[,"Y"],
                              x.test = test_model_mat,
                              ndpost = 1000, nskip = 1000, keeptrees = TRUE))
rmse_train["dbarts"] <- 
  sqrt( mean( (mu[train_vertices] - dbart_fit$yhat.test.mean[train_vertices])^2 ))
rmse_test["dbarts"] <- 
  sqrt( mean( (mu[test_vertices] - dbart_fit$yhat.test.mean[test_vertices])^2 ))
timing["dbarts"] <- dbart_time["elapsed"]

bart_time <-
  system.time(
    bart_fit <- BART::wbart(x.train = train_data[,colnames(train_data) != "Y"],
                            y.train = train_data[,"Y"],
                            x.test = test_data[,colnames(test_data) != "Y"],
                            sparse = TRUE,
                            ndpost = 1000, nskip = 1000))

rmse_train["bart"] <- 
  sqrt( mean( (mu[train_vertices] - bart_fit$yhat.test.mean[train_vertices])^2 ))
rmse_test["bart"] <- 
  sqrt( mean( (mu[test_vertices] - bart_fit$yhat.test.mean[test_vertices])^2 ))
timing["bart"] <- bart_time["elapsed"]


