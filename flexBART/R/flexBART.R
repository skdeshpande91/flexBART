flexBART <- function(Y_train,
                     X_cont_train,
                     X_cat_train,
                     X_cont_test,
                     X_cat_test,
                     unif_cuts,
                     cutpoints_list,
                     cat_levels_list,
                     edge_mat_list,
                     graph_split,
                     graph_cut_type,
                     perc_rounds, perc_threshold,
                     sparse,
                     M,
                     nd, burn, thin,
                     save_trees, verbose, print_every)
{
  y_mean <- mean(Y_train)
  y_sd <- sd(Y_train)
  std_Y_train <- (Y_train - y_mean)/y_sd # standardize the output
  tau <- (max(std_Y_train) - min(std_Y_train))/(2 * 2 * sqrt(M)) # CGM10 prior sd on all leaf parameters
  nu <- 3
  lambda <- qchisq(0.1, df = nu)/nu
  
  
  fit <- .flexBART_fit(Y_train = std_Y_train,
                       tX_cont_train = t(X_cont_train),
                       tX_cat_train = t(X_cat_train),
                       tX_cont_test = t(X_cont_test),
                       tX_cat_test = t(X_cat_test),
                       unif_cuts = unif_cuts,
                       cutpoints_list = cutpoints_list,
                       cat_levels_list = cat_levels_list,
                       edge_mat_list = edge_mat_list,
                       graph_split = graph_split,
                       graph_cut_type = graph_cut_type,
                       perc_rounds = perc_rounds, perc_threshold = perc_threshold,
                       rc_split = FALSE, prob_rc = 0, a_rc = 1, b_rc = 1,
                       sparse = sparse, a_u = 0.5, b_u = 1,
                       mu0 = 0, tau = tau, 
                       lambda = lambda, nu = nu,
                       M = M, nd = nd, burn = burn, thin = thin,
                       save_trees = save_trees, verbose = verbose, 
                       print_every = print_every)

  yhat_train_samples <- y_mean + y_sd * fit$fit_train
  yhat_test_samples <- y_mean + y_sd * fit$fit_test
  sigma_samples <- y_sd * fit$sigma
  
  results <- list()
  results[["y_mean"]] <- y_mean
  results[["y_sd"]] <- y_sd
  results[["yhat_train"]] <- y_mean + y_sd * fit$fit_train
  results[["yhat_test"]] <- y_mean + y_sd * fit$fit_test
  results[["sigma"]] <- y_sd * fit$fit_test
  results[["varcounts"]] <- fit$var_count
  return(results)
}



