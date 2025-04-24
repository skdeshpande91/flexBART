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
#sourceCpp("../flexBART/src/flexBART_fit.cpp")
sourceCpp("../flexBART/src/vcbart_fit.cpp")

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


R <- ncol(tmp_data$training_info$Z)
y_range <- 
  max(tmp_data$training_info$std_Y) - min(tmp_data$training_info$std_Y)

hyper <- parse_hyper(R = R, y_range = y_range,
                     M = 50, graph_cut_type = 2)
control <- parse_controls(save_samples = FALSE, save_trees = FALSE)

######################
# flexBART single ensemble fit
flex_fit1 <-
  .flexBART_fit(Y_train = tmp_data$training_info$std_Y,
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
                sparse = FALSE, a_u = hyper$a_u, b_u = hyper$b_u,
                sigest = hyper$sigest, lambda = hyper$lambda, nu = hyper$nu,
                M = hyper$M_vec[1],
                mu0 = hyper$mu0_vec[1], tau = hyper$tau_vec[1],
                nd = control$nd, burn = control$burn, thin = control$thin,
                save_samples = control$save_samples, save_trees = control$save_trees,
                verbose = control$verbose, print_every = control$print_every)
muhat_flex1 <- flex_fit1$fit_test_mean * sd(train_data$Y) + mean(train_data$Y)

flex_fit2 <- 
  .flexBART_fit(Y_train = tmp_data$training_info$std_Y,
                tZ_train = t(tmp_data$training_info$Z),
                tX_cont_train = t(tmp_data$training_info$X_cont),
                tX_cat_train = t(tmp_data$training_info$X_cat),
                tZ_test = t(tmp_data$testing_info$Z),
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
muhat_flex2 <- flex_fit2$fit_test_mean * sd(train_data$Y) + mean(train_data$Y)


vc_fit1 <- 
  ._vcbart_fit(Y_train = tmp_data$training_info$std_Y,
               tZ_train = t(tmp_data$training_info$Z),
               tX_cont_train = t(tmp_data$training_info$X_cont),
               tX_cat_train = t(tmp_data$training_info$X_cat),
               sigest = hyper$sigest, 
               cov_ensm = cov_ensm,
               cutpoints_list = tmp_data$training_info$cutpoints,
               cat_levels_list = tmp_data$training_info$cat_levels_list,
               edge_mat_list = tmp_data$training_info$edge_mat_list,
               nest_list = tmp_data$training_info$nest_list,
               tZ_test = t(tmp_data$testing_info$Z),
               tX_cont_test = t(tmp_data$testing_info$X_cont),
               tX_cat_test = t(tmp_data$testing_info$X_cat),
               M_vec = hyper$M_vec,
               alpha_vec = hyper$alpha_vec, beta_vec = hyper$beta_vec,
               mu0 = hyper$mu0_vec, tau = hyper$tau_vec,
               graph_cut_type = hyper$graph_cut_type,
               nest_v = hyper$nest_v, nest_v_option = hyper$nest_v_option,
               nest_c = hyper$nest_c,
               sparse = FALSE, a_u = hyper$a_u, b_u = hyper$b_u,
               nu = hyper$nu,lambda = hyper$lambda, 
               nd = control$nd, burn = control$burn, thin = control$thin,
               save_samples = control$save_samples, save_trees = control$save_trees,
               verbose = control$verbose, print_every = control$print_every)
muhat_vc1 <- vc_fit1$fit_test_mean * sd(train_data$Y) + mean(train_data$Y)


vc_fit2 <- 
  ._vcbart_fit(Y_train = tmp_data$training_info$std_Y,
               tZ_train = t(tmp_data$training_info$Z),
               tX_cont_train = t(tmp_data$training_info$X_cont),
               tX_cat_train = t(tmp_data$training_info$X_cat),
               sigest = hyper$sigest, 
               cov_ensm = cov_ensm,
               cutpoints_list = tmp_data$training_info$cutpoints,
               cat_levels_list = tmp_data$training_info$cat_levels_list,
               edge_mat_list = NULL,
               nest_list = tmp_data$training_info$nest_list,
               tZ_test = t(tmp_data$testing_info$Z),
               tX_cont_test = t(tmp_data$testing_info$X_cont),
               tX_cat_test = t(tmp_data$testing_info$X_cat),
               M_vec = hyper$M_vec,
               alpha_vec = hyper$alpha_vec, beta_vec = hyper$beta_vec,
               mu0 = hyper$mu0_vec, tau = hyper$tau_vec,
               graph_cut_type = hyper$graph_cut_type,
               nest_v = hyper$nest_v, nest_v_option = hyper$nest_v_option,
               nest_c = hyper$nest_c,
               sparse = FALSE, a_u = hyper$a_u, b_u = hyper$b_u,
               nu = hyper$nu,lambda = hyper$lambda, 
               nd = control$nd, burn = control$burn, thin = control$thin,
               save_samples = control$save_samples, save_trees = control$save_trees,
               verbose = control$verbose, print_every = control$print_every)
muhat_vc2 <- vc_fit2$fit_test_mean * sd(train_data$Y) + mean(train_data$Y)


wbart_fit <-
  BART::wbart(x.train = train_data[,"X1"],
              y.train = train_data[,"Y"],
              x.test = test_data)
muhat_wbart <- wbart_fit$yhat.test.mean




sqrt(mean( (mu[train_vertices] - muhat_flex1[train_vertices])^2 ))
sqrt(mean( (mu[train_vertices] - muhat_vc1[train_vertices])^2 ))

sqrt(mean( (mu[train_vertices] - muhat_flex2[train_vertices])^2 ))
sqrt(mean( (mu[train_vertices] - muhat_vc2[train_vertices])^2 ))
sqrt(mean( (mu[train_vertices] - muhat_wbart[train_vertices])^2 ))

sqrt(mean( (mu[test_vertices] - muhat_flex1[test_vertices])^2 ))
sqrt(mean( (mu[test_vertices] - muhat_vc1[test_vertices])^2 ))

sqrt(mean( (mu[test_vertices] - muhat_flex2[test_vertices])^2 ))
sqrt(mean( (mu[test_vertices] - muhat_vc2[test_vertices])^2 ))

sqrt(mean( (mu[test_vertices] - muhat_wbart[test_vertices])^2 ))


######################
# Plot the differences


mu_lim <- 
  c(-1,1) * 
  max(abs(c(muhat1, muhat2, muhat_wbart)))

g_flex1 <- g
g_flex2 <- g
g_bart <- g
g_heldout <- g

scaled_mu <- scales::rescale(mu, to = c(0,1), from = mu_lim)
scaled_flex1 <- scales::rescale(muhat1, to = c(0,1), from = mu_lim)
scaled_flex2 <- scales::rescale(muhat2, to = c(0,1), from = mu_lim)
scaled_bart <- scales::rescale(muhat_wbart, to = c(0,1), from = mu_lim)

V(g)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_mu)/255)
V(g_flex1)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_flex1)/255)
V(g_flex2)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_flex2)/255)
V(g_bart)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_bart)/255)

par(mar = c(1,1,1,1), mgp = c(1.8, 0.5, 0), mfrow = c(2,2))
plot(g, layout = layout_on_grid, vertex.label = NA, main = "Truth")
plot(g_flex1, layout = layout_on_grid, vertex.label = NA,main = "flexBART-adj")
plot(g_flex2, layout = layout_on_grid, vertex.label = NA,main = "flexBART-no-adj")
plot(g_bart, layout = layout_on_grid, vertex.label = NA,main = "BART")

