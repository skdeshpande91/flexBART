//
//  rflexbart.cpp
//  
//
//  Created by Sameer Deshpande on 1/20/22.
//

#include "draw_tree.h"

// [[Rcpp::export]]
Rcpp::List rflex_bart(int M,
                      int p_cont,
                      int p_cat,
                      double prob_aa, double prob_rc,
                      Rcpp::Nullable<Rcpp::List> cat_levels_list,
                      Rcpp::Nullable<Rcpp::List> adj_support_list)
{
  
  Rcpp::RNGScope scope;
  RNG gen;
  
  set_str_conversion set_str; // for converting sets of integers into string
  
  int p = p_cont + p_cat;
  
  // format the categorical levels
  std::vector<std::set<int>> cat_levels;
  std::vector<int> K; // number of levels for the different categorical variables
  std::vector<std::vector<unsigned int>> adj_support;
  
  if(p_cat > 0){
    if(cat_levels_list.isNotNull() && adj_support_list.isNotNull()){
      Rcpp::List tmp_cat_levels = Rcpp::List(cat_levels_list);
      Rcpp::List tmp_adj_support = Rcpp::List(adj_support_list);
      /*
      if( (tmp_cat_levels.size() == p_cat) && (tmp_adj_support.size() == p_cat) ){
        for(int j = 0; j < p_cat; j++){
          Rcpp::IntegerVector levels_vec = Rcpp::as<Rcpp::IntegerVector>(tmp_cat_levels[j]);
          std::set<int> levels_set;
          for(int l = 0; l < levels_vec.size(); l++) levels_set.insert(levels_vec[l]);
          cat_levels.push_back(levels_set);
          K.push_back(levels_set.size());
          
          Rcpp::IntegerVector adj_rvec = Rcpp::as<Rcpp::IntegerVector>(tmp_adj_support[j]);
          std::vector<unsigned int> adj_uvec;
          for(int l = 0; l < adj_rvec.size(); l++) adj_uvec.push_back( (unsigned int) adj_rvec[l]);
          adj_support.push_back(adj_uvec);
        }
      } else{
        Rcpp::Rcout << "p_cat = " << p_cat;
        Rcpp::Rcout << "cat_levels_list.size() = " << tmp_cat_levels.size();
        Rcpp::Rcout << "adj_levels_list.size() = " << tmp_adj_support.size() << std::endl;
        Rcpp::stop("cat_levels_list adj_levels_list must both have size equal to p_cat!");
      }
      */
      parse_categorical(cat_levels, adj_support, K, p_cat, tmp_cat_levels, tmp_adj_support);
    }
    
    
    // sanity check:
    //for(int j = 0; j < p_cat; j++){
    //  Rcpp::Rcout << " cat. var " << j << " has " << K[j] << " levels " << " and " << adj_support[j].size() << " edges" << std::endl;
    //}
    
  }

  
  data_info di;
  di.p_cont = p_cont;
  di.p_cat = p_cat;
  di.p = p;
  di.K = &K;
  di.cat_levels = &cat_levels;
  di.adj_support = &adj_support;
  
  tree_prior_info tree_pi;
  if(p_cont == 0){
    // no continuous variables so no axis aligned or random combination rules
    tree_pi.prob_aa = 0.0;
    tree_pi.prob_rc = 0.0;
  } else{
    if(p_cat == 0){
      // no categorical variables, ensure that prob_aa + prob_rc = 1
      tree_pi.prob_aa = prob_aa/(prob_aa + prob_rc);
      tree_pi.prob_rc = prob_rc/(prob_aa + prob_rc);
    } else{
      tree_pi.prob_aa = prob_aa;
      tree_pi.prob_rc = prob_rc;
    }
  }
  
  std::vector<double> theta_aa;
  double theta_rc = 0.0;
  std::vector<double> theta_cat;
  if(p_cont > 0){
    theta_aa.resize(p_cont, 1.0/( (double) p_cont));
    theta_rc = 2.0/( (double) p_cont);
    tree_pi.theta_aa = &theta_aa;
    tree_pi.theta_rc = &theta_rc;
    
  }
  if(p_cat > 0){
    theta_cat.resize(p_cat, 1.0/ ( (double) p_cat) );
    tree_pi.theta_cat = &theta_cat;
  }
  
  std::vector<tree> t_vec(M);
  Rcpp::CharacterVector tree_string_vec(M);
  for(int m = 0; m < M; m++){
    draw_tree(t_vec[m], di, tree_pi, true, gen);
    t_vec[m].print();
    tree_string_vec[m] = write_tree(t_vec[m], di, set_str);
  }
  
  
  Rcpp::List results;
  results["trees"] = tree_string_vec;
  return results;
  
}
