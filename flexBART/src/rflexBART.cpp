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
                    Rcpp::LogicalVector unif_cuts,
                    Rcpp::Nullable<Rcpp::List> cutpoints_list,
                    Rcpp::Nullable<Rcpp::List> cat_levels_list,
                    Rcpp::LogicalVector graph_split, int graph_cut_type,
                    Rcpp::Nullable<Rcpp::List> adj_support_list,
                    bool rc_split, double prob_rc,
                    double alpha, double beta,
                    double mu0, double tau)
{
  Rcpp::RNGScope scope;
  RNG gen;
  set_str_conversion set_str; // for converting sets of integers into string
  
  int n = 0;
  int p_cont = 0;
  int p_cat = 0;
  parse_training_data(n,p_cont, p_cat, tX_cont, tX_cat);
  int p = p_cont + p_cat;
  //Rcpp::Rcout << "n = " << n << " p_cont = " << p_cont << " p_cat" << std::endl;
  
  // format the categorical levels
  std::vector<std::set<double>> cutpoints;
  std::vector<std::set<int>> cat_levels;
  std::vector<int> K; // number of levels for the different categorical variables
  std::vector<std::vector<unsigned int>> adj_support;
  
  if(p_cont > 0){
    if(cutpoints_list.isNotNull()){
      Rcpp::List tmp_cutpoints = Rcpp::List(cutpoints_list);
      parse_cutpoints(cutpoints, p_cont, tmp_cutpoints, unif_cuts);
    }
  }
  
  if(p_cat > 0){
    if(cat_levels_list.isNotNull()){
      Rcpp::List tmp_cat_levels = Rcpp::List(cat_levels_list);
      parse_cat_levels(cat_levels, K, p_cat, tmp_cat_levels);
    } else{
      Rcpp::stop("Must provide categorical levels.");
    }
    if(adj_support_list.isNotNull()){
      Rcpp::List tmp_adj_support = Rcpp::List(adj_support_list);
      parse_cat_adj(adj_support, p_cat, tmp_adj_support, graph_split);
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
  
  
  std::vector<double> theta(p, 1.0/ (double) p);
  double theta_rc = 0.0;
  if(p_cont >= 2 && rc_split){
    theta_rc = 2.0/( (double) p_cont);
  }
  
  tree_prior_info tree_pi;
  tree_pi.theta = &theta;
  tree_pi.unif_cuts = unif_cuts.begin();
  tree_pi.cutpoints = &cutpoints;
  tree_pi.graph_split = graph_split.begin();
  tree_pi.graph_cut_type = graph_cut_type;
  tree_pi.prob_rc = &prob_rc;
  tree_pi.theta_rc = &theta_rc;
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
  //results["trees"] = tree_string;
  return results;
  
}


// [[Rcpp::export(".draw_ensemble")]]
Rcpp::List drawEnsemble(Rcpp::NumericMatrix tX_cont,
                        Rcpp::IntegerMatrix tX_cat,
                        Rcpp::LogicalVector unif_cuts,
                        Rcpp::Nullable<Rcpp::List> cutpoints_list,
                        Rcpp::Nullable<Rcpp::List> cat_levels_list,
                        Rcpp::LogicalVector graph_split, int graph_cut_type,
                        Rcpp::Nullable<Rcpp::List> adj_support_list,
                        bool rc_split, double prob_rc,
                        double alpha, double beta,
                        double mu0, double tau, int M,
                        bool verbose, int print_every)
{
  Rcpp::RNGScope scope;
  RNG gen;
  set_str_conversion set_str; // for converting sets of integers into string
  
  int n = 0;
  int p_cont = 0;
  int p_cat = 0;
  parse_training_data(n,p_cont, p_cat, tX_cont, tX_cat);
  int p = p_cont + p_cat;
  //Rcpp::Rcout << "n = " << n << " p_cont = " << p_cont << " p_cat" << std::endl;
  
  // format the categorical levels
  std::vector<std::set<double>> cutpoints;
  std::vector<std::set<int>> cat_levels;
  std::vector<int> K; // number of levels for the different categorical variables
  std::vector<std::vector<unsigned int>> adj_support;
  
  if(p_cont > 0){
    if(cutpoints_list.isNotNull()){
      Rcpp::List tmp_cutpoints = Rcpp::List(cutpoints_list);
      parse_cutpoints(cutpoints, p_cont, tmp_cutpoints, unif_cuts);
    }
  }
  
  if(p_cat > 0){
    if(cat_levels_list.isNotNull()){
      Rcpp::List tmp_cat_levels = Rcpp::List(cat_levels_list);
      parse_cat_levels(cat_levels, K, p_cat, tmp_cat_levels);
    } else{
      Rcpp::stop("Must provide categorical levels.");
    }
    if(adj_support_list.isNotNull()){
      Rcpp::List tmp_adj_support = Rcpp::List(adj_support_list);
      parse_cat_adj(adj_support, p_cat, tmp_adj_support, graph_split);
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
  
  
  std::vector<double> theta(p, 1.0/ (double) p);
  double theta_rc = 0.0;
  if(p_cont >= 2 && rc_split){
    theta_rc = 2.0/( (double) p_cont);
  }
  
  tree_prior_info tree_pi;
  tree_pi.theta = &theta;
  tree_pi.unif_cuts = unif_cuts.begin();
  tree_pi.cutpoints = &cutpoints;
  tree_pi.graph_split = graph_split.begin();
  tree_pi.graph_cut_type = graph_cut_type;
  tree_pi.prob_rc = &prob_rc;
  tree_pi.theta_rc = &theta_rc;
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
    if(verbose && m % print_every == 0) Rcpp::Rcout << "Drawing tree " << m+1 << " of " << M << std::endl; 
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
