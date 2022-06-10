#include "funs.h"
// TO DO: rewrite to work with the new graph stuff

// [[Rcpp::export(".predict_flexBART")]]
Rcpp::NumericMatrix predict_flexBART(Rcpp::List tree_draws,
                                     Rcpp::NumericMatrix tX_cont,
                                     Rcpp::IntegerMatrix tX_cat,
                                     bool probit = false,
                                     bool verbose = true, int print_every = 50)
{
  set_str_conversion set_str; // for converting sets of integers into strings

  int n = 0;
  int p_cont = 0;
  int p_cat = 0;
  
  parse_training_data(n, p_cont, p_cat, tX_cont, tX_cat);
  int p = p_cont + p_cat;
  data_info di;
  di.n = n;
  di.p_cont = p_cont;
  di.p_cat = p_cat;
  di.p = p;
  if(p_cont > 0) di.x_cont = tX_cont.begin();
  if(p_cat > 0) di.x_cat = tX_cat.begin();
  
  int nd = tree_draws.size();
  Rcpp::CharacterVector first_tree_vec = tree_draws[0];
  int M = first_tree_vec.size();
  
  Rcpp::Rcout << "nd = " << nd << "M = " << M;
  Rcpp::Rcout << " n = " << n << " p_cont = " << p_cont << " p_cat = " << p_cat << std::endl;
  
  std::vector<double> allfit(n);
  Rcpp::NumericMatrix pred_out(n,nd);
  
  for(int iter = 0; iter < nd; iter++){
    if( (iter%print_every == 0)){
      Rcpp::Rcout << "  Iteration: " << iter << " of " << nd <<std::endl;
      Rcpp::checkUserInterrupt();
    }
    Rcpp::CharacterVector tmp_string_vec = tree_draws[iter];
    if(tmp_string_vec.size() != M){
      // did we somehow not record enough tree strings?
      // this should really never be hit
      // essentially we're mimicing the R code all(sapply(tree_draws, FUN = length) = length(tree_draws[1]))
      Rcpp::Rcout << "iter = " << iter << " # tree strings = " << tmp_string_vec.size() << std::endl;
      Rcpp::stop("Unexpected number of tree strings!");
    } else{
      std::vector<tree> t_vec(M);
      for(int m = 0; m < M; m++){
        // tmp_string_vec is an Rcpp::CharacterVector
        // let's extract a single element from the CharacterVector and turn it into a std::string
        // that can be passed to read_tree
        std::string tmp_string = Rcpp::as<std::string>(tmp_string_vec[m]); // convert content of the
        read_tree(t_vec[m], tmp_string, set_str);
      }
      fit_ensemble(allfit, t_vec, di);
      if(probit){
        for(int i = 0; i < n; i++) pred_out(i,iter) = R::pnorm(allfit[i], 0.0, 1.0, true, false);
      } else{
        for(int i = 0; i < n; i++) pred_out(i,iter) = allfit[i];
      } 
    } // closes if/else checking that we have M strings for the draw of the ensemble
  } // closes loop over all draws of the ensemble
 
  return pred_out;
}
