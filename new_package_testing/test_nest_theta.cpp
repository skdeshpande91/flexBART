// X0,X1: continuous
// X2-X4: classroom, school, district
// X5-X6: blkgrp, tract
// X7: race

#include "../flexBART/src/rng.h"
#include "../flexBART/src/tree.h"
#include "../flexBART/src/graph_funs.h"
#include "../flexBART/src/data_parsing_funs.h"
#include "../flexBART/src/rule_funs.h"

// [[Rcpp::export]]
void test_nest_theta(Rcpp::IntegerMatrix cov_ensm,
                     Rcpp::NumericMatrix tX_cont_train,
                     Rcpp::IntegerMatrix tX_cat_train,
                     Rcpp::Nullable<Rcpp::List> cutpoints_list,
                     Rcpp::Nullable<Rcpp::List> cat_levels_list,
                     Rcpp::Nullable<Rcpp::List> edge_mat_list,
                     Rcpp::Nullable<Rcpp::List> nest_list)
{
  Rcpp::RNGScope scope;
  RNG gen;
  
  int R = cov_ensm.cols(); // how many ensembles
  int p_cont = 0;
  int p_cat = 0;
  if(tX_cont_train.size() > 1){
    p_cont = tX_cont_train.rows(); // how many continuous covariates
  }
  if(tX_cat_train.size() > 1){
    p_cat = tX_cat_train.rows(); // how many categorical covariates
  }
  int p = p_cont + p_cat;
  
  // BEGIN: set cutpoints & categorical levels + parse network structure
  std::vector<std::set<double>> cutpoints;
  if(p_cont > 0) parse_cutpoints(cutpoints, p_cont, cutpoints_list);
  
  std::vector<std::set<int>> cat_levels;
  std::vector<std::vector<edge>> edges;
  if(p_cat > 0){
    parse_cat_levels(cat_levels, p_cat, cat_levels_list);
    parse_graphs(edges, p_cat, edge_mat_list);
  }
  // END: set cutpoints & categorical levels + parse network structure

  // BEGIN: build graph encoding nesting relationships b/w categorical predictors
  std::vector<hi_lo_map> nesting;
  std::vector<edge_map> nest_graph_in;
  std::vector<edge_map> nest_graph_out;
  std::vector<std::map<int, std::set<int>>> nest_graph_components;
  parse_nesting(nesting, nest_graph_in, nest_graph_out, nest_graph_components, p_cont, cov_ensm, cat_levels, nest_list);
  // END: build graph encoding nesting relationships b/w categorical predictors

  
  tree_prior_info tree_pi;
  if(p_cont > 0) tree_pi.cutpoints = &cutpoints;
  if(p_cat > 0){
    tree_pi.cat_levels = &cat_levels;
    tree_pi.edges = &edges;
    tree_pi.nesting = &nesting;
    tree_pi.nest_in = &(nest_graph_in[0]);
    tree_pi.nest_out = &(nest_graph_out[0]);
    tree_pi.nest_components = &(nest_graph_components[0]);
  }
  tree_pi.nest_v = true;
  tree_pi.nest_c = false; // for the purposes of computed nest_theta, just use random cutset
  
  tree t;
  std::vector<double> nest_theta(p, 0.0);
  int nid = 1;
  Rcpp::Rcout << "At node " << nid << std::endl;
  for(int nvo = 0; nvo < 3; ++nvo){
    nest_theta.clear();
    nest_theta.resize(p, 0.0);
    tree_pi.nest_v_option = nvo;
    compute_nested_theta(nest_theta, t, nid, p_cont, p_cat, tree_pi);
    Rcpp::Rcout << "  " << tree_pi.nest_v_option << ":";
    for(std::vector<double>::iterator it = nest_theta.begin(); it != nest_theta.end(); ++it){
      Rcpp::Rcout << " " << *it;
    }
    Rcpp::Rcout << std::endl;
  }
  
  rule_t rule;
  rule.is_cat = true;
  rule.v_cat = 1; // school is Xcat1
  std::set<int> avail_levels;
  t.get_ptr(nid)->get_rg_cat(avail_levels, rule.v_cat);
  if(avail_levels.size() <= 1) avail_levels = tree_pi.cat_levels->at(rule.v_cat);
  partition_levels(rule, avail_levels, tree_pi, gen);
  t.birth(nid, rule);
  
  nid = 2;
  Rcpp::Rcout << "At node " << nid << std::endl;
  for(int nvo = 0; nvo < 3; ++nvo){
    nest_theta.clear();
    nest_theta.resize(p, 0.0);
    tree_pi.nest_v_option = nvo;
    compute_nested_theta(nest_theta, t, nid, p_cont, p_cat, tree_pi);
    Rcpp::Rcout << "  " << tree_pi.nest_v_option << ":";
    for(std::vector<double>::iterator it = nest_theta.begin(); it != nest_theta.end(); ++it){
      Rcpp::Rcout << " " << *it;
    }
    Rcpp::Rcout << std::endl;
  }
  
  
}
