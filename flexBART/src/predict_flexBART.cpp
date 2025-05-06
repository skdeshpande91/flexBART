#include "funs.h"

// [[Rcpp::export(".single_ensm_predict")]]
Rcpp::NumericMatrix predict_flexBART(Rcpp::List tree_draws,
                                     Rcpp::NumericMatrix tX_cont,
                                     Rcpp::IntegerMatrix tX_cat,
                                     bool probit,
                                     bool verbose, int print_every)
{
  set_str_conversion set_str; // for converting sets of integers into strings
  int n = 0;
  int p_cont = 0;
  int p_cat = 0;
  
  if(tX_cont.size() > 1) p_cont = tX_cont.rows();
  if(tX_cat.size() > 1) p_cat = tX_cat.rows();
  
  if(p_cont > 0 && p_cat == 0) n = tX_cont.cols();
  else if(p_cont == 0 && p_cat > 0) n = tX_cat.cols();
  else if(p_cont > 0 && p_cat > 0) n = tX_cont.cols();
  else Rcpp::stop("Need to supply tX_cont or tX_cat!");
  
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

  std::vector<double> allfit(n);
  Rcpp::NumericMatrix pred_out(nd,n);
  
  if(probit){
    for(int iter = 0; iter < nd; ++iter){
      if(iter == 0 || iter == nd-1) Rcpp::Rcout << "  Iteration: " << iter+1 << " of " << nd << std::endl;
      else if(iter%print_every == 0){
        Rcpp::Rcout << "  Iteration: " << iter << " of " << nd << std::endl;
        Rcpp::checkUserInterrupt();
      }
      Rcpp::CharacterVector tmp_string_vec = tree_draws[iter];
      if(tmp_string_vec.size() != M){
        Rcpp::Rcout << "iter = " << iter << " # tree strings = " << tmp_string_vec.size() << std::endl;
        Rcpp::stop("Unexpected number of tree strings!");
      } else{
        std::vector<tree> t_vec(M);
        for(int m = 0; m < M; ++m){
          // tmp_string_vec is an Rcpp::CharacterVector
          // let's extract a single element from the CharacterVector and turn it into a std::string
          // that can be passed to read_tree
          std::string tmp_string = Rcpp::as<std::string>(tmp_string_vec[m]); // convert content of the
          read_tree(t_vec[m], tmp_string, set_str);
        } // closes loop populating vector of trees
        fit_ensemble(allfit, t_vec, di);
        for(int i = 0; i < n; ++i) pred_out(iter,i) = R::pnorm(allfit[i], 0.0, 1.0, true, false);
      } // closes if/else checking that we have string for every tree
    } // closes loop over tree samples
  } else{
    for(int iter = 0; iter < nd; ++iter){
      if(iter == 0 || iter == nd-1) Rcpp::Rcout << "  Iteration: " << iter+1 << " of " << nd << std::endl;
      else if(iter%print_every == 0){
        Rcpp::Rcout << "  Iteration: " << iter << " of " << nd << std::endl;
        Rcpp::checkUserInterrupt();
      }
      Rcpp::CharacterVector tmp_string_vec = tree_draws[iter];
      if(tmp_string_vec.size() != M){
        Rcpp::Rcout << "iter = " << iter << " # tree strings = " << tmp_string_vec.size() << std::endl;
        Rcpp::stop("Unexpected number of tree strings!");
      } else{
        std::vector<tree> t_vec(M);
        for(int m = 0; m < M; ++m){
          // tmp_string_vec is an Rcpp::CharacterVector
          // let's extract a single element from the CharacterVector and turn it into a std::string
          // that can be passed to read_tree
          std::string tmp_string = Rcpp::as<std::string>(tmp_string_vec[m]); // convert content of the
          read_tree(t_vec[m], tmp_string, set_str);
        } // closes loop populating vector of trees
        fit_ensemble(allfit, t_vec, di);
        for(int i = 0; i < n; ++i) pred_out(iter,i) = allfit[i];
      } // closes if/else checking that we have string for every tree
    } // closes loop over tree samples
  } // closes if/else checking whether we're doing probit
  return pred_out;
}
