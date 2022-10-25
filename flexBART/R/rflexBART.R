rflexBART <- function(M, 
                      X_cont = matrix(0, nrow = 1, ncol = 1),
                      X_cat = matrix(0L, nrow = 1, ncol = 1),
                      unif_cuts = rep(TRUE, times = ncol(X_cont)),
                      cutpoints_list = NULL,
                      cat_levels_list = NULL,
                      a_cat = 0, b_cat = 0,
                      alpha = 0.95, beta = 2, mu0 = 0, tau = 1/sqrt(M),
                      verbose = TRUE, print_every = floor(M/10))
{
  tree_draws <- .drawEnsemble(tX_cont = t(X_cont),
                              tX_cat = t(X_cat),
                              unif_cuts = unif_cuts,
                              cutpoints_list = cutpoints_list,
                              cat_levels_list = cat_levels_list,
                              edge_mat_list = NULL,
                              graph_split = rep(FALSE, times = ncol(X_cat)),
                              graph_cut_type = 0,
                              a_cat = a_cat, b_cat = b_cat,
                              perc_rounds = 0, perc_threshold = 0.5,
                              rc_split = FALSE, prob_rc = 0, a_rc = 1, b_rc = 1,
                              alpha = alpha, beta = beta,
                              mu0 = mu0, tau = tau, M = M,
                              verbose = verbose, print_every = print_every)
  return(tree_draws)
}