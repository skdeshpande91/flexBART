rflexBART <- function(M, 
                      X_cont,
                      X_cat,
                      unif_cuts,
                      cutpoints_list,
                      cat_levels_list,
                      edge_mat_list,
                      graph_split,
                      graph_cut_type,
                      perc_rounds, perc_threshold,
                      alpha, beta, mu0, tau,
                      verbose = verbose, print_every = print_every)
{
  tree_draws <- .drawEnsemble(tX_cont = t(X_cont),
                              tX_cat = t(X_cat),
                              unif_cuts = unif_cuts,
                              cutpoints_list = cutpoints_list,
                              cat_levels_list = cat_levels_list,
                              edge_mat_list = edge_mat_list,
                              graph_split = graph_split,
                              graph_cut_type = graph_cut_type,
                              perc_rounds = perc_rounds, perc_threshold = perc_threshold,
                              rc_split = FALSE, prob_rc = 0, a_rc = 1, b_rc = 1,
                              alpha = alpha, beta = beta,
                              mu0 = mu0, tau = tau, M = M,
                              verbose = verbose, print_every = print_every)
  return(tree_draws)
}