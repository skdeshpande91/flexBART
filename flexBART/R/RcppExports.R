# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

.flexBART_fit <- function(Y_train, tX_cont_train, tX_cat_train, tX_cont_test, tX_cat_test, unif_cuts, cutpoints_list, cat_levels_list, edge_mat_list, graph_split, graph_cut_type, sparse, a_u, b_u, mu0, tau, lambda, nu, M, nd, burn, thin, save_samples, save_trees, verbose, print_every) {
    .Call('_flexBART_flexBART_fit', PACKAGE = 'flexBART', Y_train, tX_cont_train, tX_cat_train, tX_cont_test, tX_cat_test, unif_cuts, cutpoints_list, cat_levels_list, edge_mat_list, graph_split, graph_cut_type, sparse, a_u, b_u, mu0, tau, lambda, nu, M, nd, burn, thin, save_samples, save_trees, verbose, print_every)
}

.predict_flexBART <- function(tree_draws, tX_cont, tX_cat, probit = FALSE, verbose = TRUE, print_every = 50L) {
    .Call('_flexBART_predict_flexBART', PACKAGE = 'flexBART', tree_draws, tX_cont, tX_cat, probit, verbose, print_every)
}

.probit_flexBART_fit <- function(Y_train, tX_cont_train, tX_cat_train, tX_cont_test, tX_cat_test, unif_cuts, cutpoints_list, cat_levels_list, edge_mat_list, graph_split, graph_cut_type, sparse, a_u, b_u, mu0, tau, M, nd, burn, thin, save_samples, save_trees, verbose, print_every) {
    .Call('_flexBART_probit_flexBART_fit', PACKAGE = 'flexBART', Y_train, tX_cont_train, tX_cat_train, tX_cont_test, tX_cat_test, unif_cuts, cutpoints_list, cat_levels_list, edge_mat_list, graph_split, graph_cut_type, sparse, a_u, b_u, mu0, tau, M, nd, burn, thin, save_samples, save_trees, verbose, print_every)
}

.drawEnsemble <- function(tX_cont, tX_cat, unif_cuts, cutpoints_list, cat_levels_list, edge_mat_list, graph_split, graph_cut_type, alpha, beta, mu0, tau, M, verbose, print_every) {
    .Call('_flexBART_drawEnsemble', PACKAGE = 'flexBART', tX_cont, tX_cat, unif_cuts, cutpoints_list, cat_levels_list, edge_mat_list, graph_split, graph_cut_type, alpha, beta, mu0, tau, M, verbose, print_every)
}

summarize_post_pred <- function(fit_samples, sigma_samples) {
    .Call('_flexBART_summarize_post_pred', PACKAGE = 'flexBART', fit_samples, sigma_samples)
}

