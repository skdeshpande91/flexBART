probit_flexBART <- function(Y_train,
                            X_cont_train = matrix(0, nrow = 1, ncol = 1),
                            X_cat_train = matrix(0, nrow = 1, ncol = 1),
                            X_cont_test = matrix(0, nrow = 1, ncol = 1),
                            X_cat_test = matrix(0L, nrow = 1, ncol = 1),
                            unif_cuts = rep(TRUE, times = ncol(X_cont_train)),
                            cutpoints_list = NULL,
                            cat_levels_list,
                            sparse = FALSE,
                            M = 200, mu0 = stats::qnorm(mean(Y_train)), tau = 1/sqrt(M),
                            nd = 1000, burn = 1000, thin = 1,
                            save_trees = TRUE, verbose = TRUE, print_every = floor( (nd*thin + burn))/10)
{
  if(!is.integer(Y_train)) stop("Y_train must be an integer vector")
  if(!all(Y_train %in% c(0,1))) stop("All elements of Y_train must be 0 or 1")
  
  
  
  fit <- .probit_flexBART_fit(Y_train = Y_train,
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
                              perc_rounds = 0, perc_threshold = 0,
                              rc_split = FALSE, prob_rc = 0, a_rc = 1, b_rc = 1,
                              sparse = sparse, a_u = 0.5, b_u = 1,
                              mu0 = mu0, tau = tau, 
                              M = M, nd = nd, burn = burn, thin = thin,
                              save_trees = save_trees, verbose = verbose, 
                              print_every = print_every)
  
  results <- list()
  results[["yhat_train"]] <- fit$fit_train
  results[["yhat_test"]] <- fit$fit_test
  results[["varcounts"]] <- fit$var_count
  if(save_trees) results[["trees"]] <- fit$trees
  return(results)
}