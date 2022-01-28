//
//  rflexbart.cpp
//  
//
//  Created by Sameer Deshpande on 1/20/22.
//

#include "draw_tree.h"


// [[Rcpp::export(".draw_tree")]]
Rcpp::List drawTree(Rcpp::NumericMatrix tX_cont,
                    Rcpp::IntegerMatrix tX_cat,
                    Rcpp::Nullable<Rcpp::List> cat_levels_list,
                    Rcpp::Nullable<Rcpp::List> adj_support_list,
                    bool mst_split, bool mst_reweight,
                    double alpha, double beta,
                    double mu0, double tau,
                    double prob_aa, double prob_rc)
{
  Rcpp::RNGScope scope;
  RNG gen;
  set_str_conversion set_str; // for converting sets of integers into string
  
  int n = 0;
  int p_cont = 0;
  int p_cat = 0;
  parse_training_data(n,p_cont, p_cat, tX_cont, tX_cat);

  int p = p_cont + p_cat;
  
  // format the categorical levels
  std::vector<std::set<int>> cat_levels;
  std::vector<int> K; // number of levels for the different categorical variables
  std::vector<std::vector<unsigned int>> adj_support;
  
  if(p_cat > 0){
    if(cat_levels_list.isNotNull() && adj_support_list.isNotNull()){
      Rcpp::List tmp_cat_levels = Rcpp::List(cat_levels_list);
      Rcpp::List tmp_adj_support = Rcpp::List(adj_support_list);
      parse_categorical(cat_levels, adj_support, K, p_cat, tmp_cat_levels, tmp_adj_support);
    }
  }

  data_info di;
  di.n = n;
  di.p_cont = p_cont;
  di.p_cat = p_cat;
  di.p = p;
  if(p_cont > 0) di.x_cont = tX_cont.begin();
  if(p_cat > 0){
    di.x_cat = tX_cat.begin();
    di.cat_levels = &cat_levels;
    di.K = &K;
    di.adj_support = &adj_support;
  }
  
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
  
  
  tree_pi.mst_split = mst_split;
  tree_pi.mst_reweight = mst_reweight;
  tree_pi.alpha = alpha;
  tree_pi.beta = beta;
  tree_pi.mu0 = mu0;
  tree_pi.tau = tau;
  
  tree t;
  suff_stat ss;
  double tmp_mu;
  Rcpp::NumericVector fit(n);
  
  draw_tree(t, di, tree_pi, gen);
  t.print();
  
  tree_traversal(ss, t, di);
  for(suff_stat_it ss_it = ss.begin(); ss_it != ss.end(); ++ss_it){
    tmp_mu = t.get_ptr(ss_it->first)->get_mu(); // get the value of mu in the leaf
    for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
      fit[*it] = tmp_mu;
    }
  }
  
  Rcpp::CharacterVector tree_string(1);
  tree_string[0] = write_tree(t, di, set_str);
  
  
  Rcpp::List results;
  results["fit"] = fit;
  results["trees"] = tree_string;
  return results;
  
}


// [[Rcpp::export(".draw_ensemble")]]
Rcpp::List drawEnsemble(Rcpp::NumericMatrix tX_cont,
                        Rcpp::IntegerMatrix tX_cat,
                        Rcpp::Nullable<Rcpp::List> cat_levels_list,
                        Rcpp::Nullable<Rcpp::List> adj_support_list,
                        bool mst_split, bool mst_reweight,
                        double alpha, double beta,
                        double mu0, double tau,
                        double prob_aa, double prob_rc, int M)
{
  Rcpp::RNGScope scope;
  RNG gen;
  set_str_conversion set_str; // for converting sets of integers into string
  
  int n = 0;
  int p_cont = 0;
  int p_cat = 0;
  parse_training_data(n,p_cont, p_cat, tX_cont, tX_cat);

  int p = p_cont + p_cat;
  
  // format the categorical levels
  std::vector<std::set<int>> cat_levels;
  std::vector<int> K; // number of levels for the different categorical variables
  std::vector<std::vector<unsigned int>> adj_support;
  
  if(p_cat > 0){
    if(cat_levels_list.isNotNull() && adj_support_list.isNotNull()){
      Rcpp::List tmp_cat_levels = Rcpp::List(cat_levels_list);
      Rcpp::List tmp_adj_support = Rcpp::List(adj_support_list);
      parse_categorical(cat_levels, adj_support, K, p_cat, tmp_cat_levels, tmp_adj_support);
    }
  }

  data_info di;
  di.n = n;
  di.p_cont = p_cont;
  di.p_cat = p_cat;
  di.p = p;
  if(p_cont > 0) di.x_cont = tX_cont.begin();
  if(p_cat > 0){
    di.x_cat = tX_cat.begin();
    di.cat_levels = &cat_levels;
    di.K = &K;
    di.adj_support = &adj_support;
  }
  
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
  
  
  tree_pi.mst_split = mst_split;
  tree_pi.mst_reweight = mst_reweight;
  tree_pi.alpha = alpha;
  tree_pi.beta = beta;
  tree_pi.mu0 = mu0;
  tree_pi.tau = tau;
  
  
  Rcpp::NumericMatrix tree_fits(n,M);
  Rcpp::IntegerMatrix leaf_id(n,M);
  Rcpp::NumericVector fit(n);
  for(int i = 0; i < n; i++) fit[i] = 0.0;
  Rcpp::CharacterVector tree_strings(M);
  
  
  tree t;
  suff_stat ss;
  double tmp_mu;
  
  for(int m = 0; m < M; m++){
    if(m % 10 == 0) Rcpp::Rcout << "Drawing tree " << m+1 << " of " << M << std::endl; 
    t.to_null();
    draw_tree(t, di, tree_pi, gen);
    ss.clear();
    tree_traversal(ss,t,di);
    for(suff_stat_it ss_it = ss.begin(); ss_it != ss.end(); ++ss_it){
      tmp_mu = t.get_ptr(ss_it->first)->get_mu();
      for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
        tree_fits(*it,m) = tmp_mu;
        fit[*it] += tmp_mu;
        leaf_id(*it,m) = ss_it->first; // id of the leaf
      }
    }
    tree_strings[m] = write_tree(t, di, set_str);
  }
  
  Rcpp::List results;
  results["fit"] = fit;
  results["trees"] = tree_strings;
  results["tree_fits"] = tree_fits;
  results["leaf"] = leaf_id;
  return results;
}

/*
// [[Rcpp::export(".rflexBART")]]
Rcpp::List rflexBART(Rcpp::NumericMatrix tX_cont,
                     Rcpp::IntegerMatrix tX_cat,
                     Rcpp::Nullable<Rcpp::List> cat_levels_list,
                     Rcpp::Nullable<Rcpp::List> adj_support_list,
                     bool mst_split, bool mst_reweight,
                     double mu0, double tau,
                     double prob_aa, double prob_rc,
                     int M, int nd, bool save_trees)
{
  
  Rcpp::RNGScope scope;
  RNG gen;
  set_str_conversion set_str; // for converting sets of integers into string
  
  int n = 0;
  int p_cont = 0;
  int p_cat = 0;
  parse_training_data(n,p_cont, p_cat, tX_cont, tX_cat);

  int p = p_cont + p_cat;
  
  // format the categorical levels
  std::vector<std::set<int>> cat_levels;
  std::vector<int> K; // number of levels for the different categorical variables
  std::vector<std::vector<unsigned int>> adj_support;
  
  if(p_cat > 0){
    if(cat_levels_list.isNotNull() && adj_support_list.isNotNull()){
      Rcpp::List tmp_cat_levels = Rcpp::List(cat_levels_list);
      Rcpp::List tmp_adj_support = Rcpp::List(adj_support_list);
      parse_categorical(cat_levels, adj_support, K, p_cat, tmp_cat_levels, tmp_adj_support);
    }
  }

  data_info di;
  di.n = n;
  di.p_cont = p_cont;
  di.p_cat = p_cat;
  di.p = p;
  if(p_cont > 0) di.x_cont = tX_cont.begin();
  if(p_cat > 0){
    di.x_cat = tX_cat.begin();
    di.cat_levels = &cat_levels;
    di.K = &K;
    di.adj_support = &adj_support;
  }
  
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
  
  
  tree_pi.mst_split = mst_split;
  tree_pi.mst_reweight = mst_reweight;
  tree_pi.mu0 = mu0;
  tree_pi.tau = tau;
  
  // now we are ready for to do the trees
  Rcpp::List tree_draws(nd);
  
  
  for(int iter = 0; iter < nd; iter++){
    for(int m = 0; m < M; m++){
      t_vec[m].to_null(); // kill off the current value of the m-th tree
      draw_tree(t_vec[m], di, tree_pi, true, gen);
    }
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
*/
