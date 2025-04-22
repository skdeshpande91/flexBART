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
sourceCpp("../flexBART/src/flexBART_fit.cpp")
#sourceCpp("reBART/src/rescale_beta.cpp")

source("generate_grid_data.R")

frmla <- Y ~ bart(X1)
tmp_form <- parse_formula(frmla, colnames(train_data))
outcome_name <- tmp_form$outcome_name
cov_ensm <- tmp_form$cov_ensm
adjacency_list <- list(X1 = A)


tmp_data <- 
  prepare_data(train_data = train_data,
               outcome_name = outcome_name, 
               cov_ensm = cov_ensm, 
               test_data = test_data,
               adjacency_list = adjacency_list)


R <- 1
y_range <- 
  max(tmp_data$training_info$std_Y) - min(tmp_data$training_info$std_Y)

hyper <- parse_hyper(R = R, y_range = y_range,
                     M = 50, graph_cut_type = 2)
control <- parse_controls(save_samples = FALSE, save_trees = FALSE)

######################

fit <-
  .flexBART_fit(Y_train = tmp_data$training_info$std_Y,
                tX_cont_train = t(tmp_data$training_info$X_cont),
                tX_cat_train = t(tmp_data$training_info$X_cat),
                tX_cont_test = t(tmp_data$testing_info$X_cont),
                tX_cat_test = t(tmp_data$testing_info$X_cat),
                cutpoints_list = tmp_data$training_info$cutpoints,
                cat_levels_list = tmp_data$training_info$cat_levels_list,
                edge_mat_list = tmp_data$training_info$edge_mat_list,
                graph_cut_type = hyper$graph_cut_type,
                sparse = FALSE, a_u = hyper$a_u, b_u = hyper$b_u,
                sigest = hyper$sigest, lambda = hyper$lambda, nu = hyper$nu,
                M = hyper$M_vec[1],
                mu0 = hyper$mu0_vec[1], tau = hyper$tau_vec[1],
                nd = control$nd, burn = control$burn, thin = control$thin,
                save_samples = control$save_samples, save_trees = control$save_trees,
                verbose = control$verbose, print_every = control$print_every)


fit2 <- 
  .flexBART_fit(Y_train = tmp_data$training_info$std_Y,
                tX_cont_train = t(tmp_data$training_info$X_cont),
                tX_cat_train = t(tmp_data$training_info$X_cat),
                tX_cont_test = t(tmp_data$testing_info$X_cont),
                tX_cat_test = t(tmp_data$testing_info$X_cat),
                cutpoints_list = tmp_data$training_info$cutpoints,
                cat_levels_list = tmp_data$training_info$cat_levels_list,
                edge_mat_list = NULL,
                graph_cut_type = hyper$graph_cut_type,
                sparse = FALSE, a_u = hyper$a_u, b_u = hyper$b_u,
                sigest = hyper$sigest, lambda = hyper$lambda, nu = hyper$nu,
                M = hyper$M_vec[1],
                mu0 = hyper$mu0_vec[1], tau = hyper$tau_vec[1],
                nd = control$nd, burn = control$burn, thin = control$thin,
                save_samples = control$save_samples, save_trees = control$save_trees,
                verbose = control$verbose, print_every = control$print_every)
