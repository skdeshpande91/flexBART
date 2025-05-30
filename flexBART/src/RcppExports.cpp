// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// flexBART_fit
Rcpp::List flexBART_fit(Rcpp::NumericVector Y_train, Rcpp::NumericMatrix tX_cont_train, Rcpp::IntegerMatrix tX_cat_train, Rcpp::NumericMatrix tX_cont_test, Rcpp::IntegerMatrix tX_cat_test, Rcpp::LogicalVector unif_cuts, Rcpp::Nullable<Rcpp::List> cutpoints_list, Rcpp::Nullable<Rcpp::List> cat_levels_list, Rcpp::Nullable<Rcpp::List> edge_mat_list, Rcpp::LogicalVector graph_split, int graph_cut_type, bool sparse, double a_u, double b_u, double mu0, double tau, double lambda, double nu, int M, int nd, int burn, int thin, bool save_samples, bool save_trees, bool verbose, int print_every);
RcppExport SEXP _flexBART_flexBART_fit(SEXP Y_trainSEXP, SEXP tX_cont_trainSEXP, SEXP tX_cat_trainSEXP, SEXP tX_cont_testSEXP, SEXP tX_cat_testSEXP, SEXP unif_cutsSEXP, SEXP cutpoints_listSEXP, SEXP cat_levels_listSEXP, SEXP edge_mat_listSEXP, SEXP graph_splitSEXP, SEXP graph_cut_typeSEXP, SEXP sparseSEXP, SEXP a_uSEXP, SEXP b_uSEXP, SEXP mu0SEXP, SEXP tauSEXP, SEXP lambdaSEXP, SEXP nuSEXP, SEXP MSEXP, SEXP ndSEXP, SEXP burnSEXP, SEXP thinSEXP, SEXP save_samplesSEXP, SEXP save_treesSEXP, SEXP verboseSEXP, SEXP print_everySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Y_train(Y_trainSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type tX_cont_train(tX_cont_trainSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type tX_cat_train(tX_cat_trainSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type tX_cont_test(tX_cont_testSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type tX_cat_test(tX_cat_testSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type unif_cuts(unif_cutsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type cutpoints_list(cutpoints_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type cat_levels_list(cat_levels_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type edge_mat_list(edge_mat_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type graph_split(graph_splitSEXP);
    Rcpp::traits::input_parameter< int >::type graph_cut_type(graph_cut_typeSEXP);
    Rcpp::traits::input_parameter< bool >::type sparse(sparseSEXP);
    Rcpp::traits::input_parameter< double >::type a_u(a_uSEXP);
    Rcpp::traits::input_parameter< double >::type b_u(b_uSEXP);
    Rcpp::traits::input_parameter< double >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type nd(ndSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< bool >::type save_samples(save_samplesSEXP);
    Rcpp::traits::input_parameter< bool >::type save_trees(save_treesSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type print_every(print_everySEXP);
    rcpp_result_gen = Rcpp::wrap(flexBART_fit(Y_train, tX_cont_train, tX_cat_train, tX_cont_test, tX_cat_test, unif_cuts, cutpoints_list, cat_levels_list, edge_mat_list, graph_split, graph_cut_type, sparse, a_u, b_u, mu0, tau, lambda, nu, M, nd, burn, thin, save_samples, save_trees, verbose, print_every));
    return rcpp_result_gen;
END_RCPP
}
// predict_flexBART
Rcpp::NumericMatrix predict_flexBART(Rcpp::List tree_draws, Rcpp::NumericMatrix tX_cont, Rcpp::IntegerMatrix tX_cat, bool probit, bool verbose, int print_every);
RcppExport SEXP _flexBART_predict_flexBART(SEXP tree_drawsSEXP, SEXP tX_contSEXP, SEXP tX_catSEXP, SEXP probitSEXP, SEXP verboseSEXP, SEXP print_everySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type tree_draws(tree_drawsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type tX_cont(tX_contSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type tX_cat(tX_catSEXP);
    Rcpp::traits::input_parameter< bool >::type probit(probitSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type print_every(print_everySEXP);
    rcpp_result_gen = Rcpp::wrap(predict_flexBART(tree_draws, tX_cont, tX_cat, probit, verbose, print_every));
    return rcpp_result_gen;
END_RCPP
}
// probit_flexBART_fit
Rcpp::List probit_flexBART_fit(Rcpp::IntegerVector Y_train, Rcpp::NumericMatrix tX_cont_train, Rcpp::IntegerMatrix tX_cat_train, Rcpp::NumericMatrix tX_cont_test, Rcpp::IntegerMatrix tX_cat_test, Rcpp::LogicalVector unif_cuts, Rcpp::Nullable<Rcpp::List> cutpoints_list, Rcpp::Nullable<Rcpp::List> cat_levels_list, Rcpp::Nullable<Rcpp::List> edge_mat_list, Rcpp::LogicalVector graph_split, int graph_cut_type, bool sparse, double a_u, double b_u, double mu0, double tau, int M, int nd, int burn, int thin, bool save_samples, bool save_trees, bool verbose, int print_every);
RcppExport SEXP _flexBART_probit_flexBART_fit(SEXP Y_trainSEXP, SEXP tX_cont_trainSEXP, SEXP tX_cat_trainSEXP, SEXP tX_cont_testSEXP, SEXP tX_cat_testSEXP, SEXP unif_cutsSEXP, SEXP cutpoints_listSEXP, SEXP cat_levels_listSEXP, SEXP edge_mat_listSEXP, SEXP graph_splitSEXP, SEXP graph_cut_typeSEXP, SEXP sparseSEXP, SEXP a_uSEXP, SEXP b_uSEXP, SEXP mu0SEXP, SEXP tauSEXP, SEXP MSEXP, SEXP ndSEXP, SEXP burnSEXP, SEXP thinSEXP, SEXP save_samplesSEXP, SEXP save_treesSEXP, SEXP verboseSEXP, SEXP print_everySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type Y_train(Y_trainSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type tX_cont_train(tX_cont_trainSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type tX_cat_train(tX_cat_trainSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type tX_cont_test(tX_cont_testSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type tX_cat_test(tX_cat_testSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type unif_cuts(unif_cutsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type cutpoints_list(cutpoints_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type cat_levels_list(cat_levels_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type edge_mat_list(edge_mat_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type graph_split(graph_splitSEXP);
    Rcpp::traits::input_parameter< int >::type graph_cut_type(graph_cut_typeSEXP);
    Rcpp::traits::input_parameter< bool >::type sparse(sparseSEXP);
    Rcpp::traits::input_parameter< double >::type a_u(a_uSEXP);
    Rcpp::traits::input_parameter< double >::type b_u(b_uSEXP);
    Rcpp::traits::input_parameter< double >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type nd(ndSEXP);
    Rcpp::traits::input_parameter< int >::type burn(burnSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< bool >::type save_samples(save_samplesSEXP);
    Rcpp::traits::input_parameter< bool >::type save_trees(save_treesSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type print_every(print_everySEXP);
    rcpp_result_gen = Rcpp::wrap(probit_flexBART_fit(Y_train, tX_cont_train, tX_cat_train, tX_cont_test, tX_cat_test, unif_cuts, cutpoints_list, cat_levels_list, edge_mat_list, graph_split, graph_cut_type, sparse, a_u, b_u, mu0, tau, M, nd, burn, thin, save_samples, save_trees, verbose, print_every));
    return rcpp_result_gen;
END_RCPP
}
// drawEnsemble
Rcpp::List drawEnsemble(Rcpp::NumericMatrix tX_cont, Rcpp::IntegerMatrix tX_cat, Rcpp::LogicalVector unif_cuts, Rcpp::Nullable<Rcpp::List> cutpoints_list, Rcpp::Nullable<Rcpp::List> cat_levels_list, Rcpp::Nullable<Rcpp::List> edge_mat_list, Rcpp::LogicalVector graph_split, int graph_cut_type, double alpha, double beta, double mu0, double tau, int M, bool verbose, int print_every);
RcppExport SEXP _flexBART_drawEnsemble(SEXP tX_contSEXP, SEXP tX_catSEXP, SEXP unif_cutsSEXP, SEXP cutpoints_listSEXP, SEXP cat_levels_listSEXP, SEXP edge_mat_listSEXP, SEXP graph_splitSEXP, SEXP graph_cut_typeSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP mu0SEXP, SEXP tauSEXP, SEXP MSEXP, SEXP verboseSEXP, SEXP print_everySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type tX_cont(tX_contSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type tX_cat(tX_catSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type unif_cuts(unif_cutsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type cutpoints_list(cutpoints_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type cat_levels_list(cat_levels_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::List> >::type edge_mat_list(edge_mat_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::LogicalVector >::type graph_split(graph_splitSEXP);
    Rcpp::traits::input_parameter< int >::type graph_cut_type(graph_cut_typeSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type mu0(mu0SEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< int >::type print_every(print_everySEXP);
    rcpp_result_gen = Rcpp::wrap(drawEnsemble(tX_cont, tX_cat, unif_cuts, cutpoints_list, cat_levels_list, edge_mat_list, graph_split, graph_cut_type, alpha, beta, mu0, tau, M, verbose, print_every));
    return rcpp_result_gen;
END_RCPP
}
// summarize_post_pred
Rcpp::NumericMatrix summarize_post_pred(Rcpp::NumericMatrix fit_samples, Rcpp::NumericVector sigma_samples);
RcppExport SEXP _flexBART_summarize_post_pred(SEXP fit_samplesSEXP, SEXP sigma_samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type fit_samples(fit_samplesSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type sigma_samples(sigma_samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(summarize_post_pred(fit_samples, sigma_samples));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_flexBART_flexBART_fit", (DL_FUNC) &_flexBART_flexBART_fit, 26},
    {"_flexBART_predict_flexBART", (DL_FUNC) &_flexBART_predict_flexBART, 6},
    {"_flexBART_probit_flexBART_fit", (DL_FUNC) &_flexBART_probit_flexBART_fit, 24},
    {"_flexBART_drawEnsemble", (DL_FUNC) &_flexBART_drawEnsemble, 15},
    {"_flexBART_summarize_post_pred", (DL_FUNC) &_flexBART_summarize_post_pred, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_flexBART(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
