flexBART <- function(Y_train,
                     X_cont_train = matrix(0, nrow = 1, ncol = 1),
                     X_cat_train = matrix(0, nrow = 1, ncol = 1),
                     X_cont_test = matrix(0, nrow = 1, ncol = 1),
                     X_cat_test = matrix(0L, nrow = 1, ncol = 1),
                     unif_cuts = rep(TRUE, times = ncol(X_cont_train)),
                     cutpoints_list = NULL,
                     cat_levels_list = NULL,
                     a_cat = 0, b_cat = 0,
                     sparse = FALSE,
                     M = 200,
                     nd = 1000, burn = 1000, thin = 1,
                     save_samples = TRUE,
                     save_trees = TRUE, verbose = TRUE, print_every = floor( (nd*thin + burn))/10)
{
  y_mean <- mean(Y_train)
  y_sd <- stats::sd(Y_train)
  std_Y_train <- (Y_train - y_mean)/y_sd # standardize the output
  tau <- (max(std_Y_train) - min(std_Y_train))/(2 * 2 * sqrt(M)) # CGM10 prior sd on all leaf parameters
  nu <- 3
  lambda <- stats::qchisq(0.1, df = nu)/nu
  
  fit <- .flexBART_fit(Y_train = std_Y_train,
                       tX_cont_train = t(X_cont_train),
                       tX_cat_train = t(X_cat_train),
                       tX_cont_test = t(X_cont_test),
                       tX_cat_test = t(X_cat_test),
                       unif_cuts = unif_cuts,
                       cutpoints_list = cutpoints_list,
                       cat_levels_list = cat_levels_list,
                       edge_mat_list = NULL,
                       graph_split = rep(FALSE, times = ncol(X_cat_train)),
                       graph_cut_type = 0,
                       a_cat = a_cat, b_cat = b_cat,
                       perc_rounds = 0, perc_threshold = 0,
                       rc_split = FALSE, prob_rc = 0, a_rc = 1, b_rc = 1,
                       sparse = sparse, a_u = 0.5, b_u = 1,
                       mu0 = 0, tau = tau, 
                       lambda = lambda, nu = nu,
                       M = M, nd = nd, burn = burn, thin = thin,
                       save_samples = save_samples,
                       save_trees = save_trees, verbose = verbose, 
                       print_every = print_every)
  
  yhat_train_mean <- y_mean + y_sd * fit$fit_train_mean
  if(save_samples){
    yhat_train_samples <- y_mean + y_sd * fit$fit_train
  }
  if(!is.null(fit$fit_test_mean)){
    yhat_test_mean <- y_mean + y_sd * fit$fit_test_mean
    if(save_samples){
      yhat_test_samples <- y_mean + y_sd * fit$fit_test
    }
  }
  sigma_samples <- y_sd * fit$sigma
  
  results <- list()
  results[["y_mean"]] <- y_mean
  results[["y_sd"]] <- y_sd
  results[["yhat.train.mean"]] <- yhat_train_mean
  if(save_samples) results[["yhat.train"]] <- yhat_train_samples
  if(!is.null(fit$fit_test_mean)){
    results[["yhat.test.mean"]] <- yhat_test_mean
    if(save_samples) results[["yhat.test"]] <- yhat_test_samples
  }
  results[["sigma"]] <- y_sd * fit$sigma
  results[["varcounts"]] <- fit$var_count
  if(save_trees) results[["trees"]] <- fit$trees
  return(results)
}