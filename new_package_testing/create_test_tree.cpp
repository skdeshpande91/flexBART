//
// Xcont 0, Xcont 1: continuous variables
// Xcat0-Xcat2: classroom school district
// Xcat3-Xcat4: blkgrp, tract
// Xcat5: race

#include "../flexBART/src/rng.h"
#include "../flexBART/src/tree.h"
#include "../flexBART/src/graph_funs.h"
#include "../flexBART/src/data_parsing_funs.h"
#include "../flexBART/src/rule_funs.h"

void create_tree(int r, tree &t, tree_prior_info &tree_pi, RNG &gen)
{
  t.to_null();
  rule_t rule;
  std::set<int> avail_levels;
  
  int nx_nid = 0;
  tree::tree_p nx = t.get_ptr(1);
  
  if(r == 0){
    // ensemble on gets classroom, school, district
    // Node 1: split on school (Xcat1)
    nx_nid = 1;
    nx = t.get_ptr(nx_nid);
    rule.clear();
    rule.is_cat = true;
    rule.v_cat = 1;
    avail_levels.clear();
    nx->get_rg_nested_cat(avail_levels, rule.v_cat, tree_pi);
    if(avail_levels.size() <= 1) avail_levels = tree_pi.cat_levels->at(rule.v_cat); // trivial split
    partition_levels(rule, avail_levels, tree_pi, gen);
    t.birth(nx_nid, rule);
    
    // node 2: split on school
    nx_nid = 2;
    nx = t.get_ptr(nx_nid);
    rule.clear();
    rule.is_cat = true;
    rule.v_cat = 1;
    avail_levels.clear();
    nx->get_rg_nested_cat(avail_levels, rule.v_cat, tree_pi);
    if(avail_levels.size() <= 1) avail_levels = tree_pi.cat_levels->at(rule.v_cat); // trivial split
    partition_levels(rule, avail_levels, tree_pi, gen);
    t.birth(nx_nid, rule);
    
    // node 3: split on district
    nx_nid = 3;
    nx = t.get_ptr(nx_nid);
    rule.clear();
    rule.is_cat = true;
    rule.v_cat = 2;
    avail_levels.clear();
    nx->get_rg_nested_cat(avail_levels, rule.v_cat, tree_pi);
    if(avail_levels.size() <= 1) avail_levels = tree_pi.cat_levels->at(rule.v_cat); // trivial split
    partition_levels(rule, avail_levels, tree_pi, gen);
    t.birth(nx_nid, rule);
    
    // node 4: split on classroom
    nx_nid = 4;
    nx = t.get_ptr(nx_nid);
    rule.clear();
    rule.is_cat = true;
    rule.v_cat = 0;
    avail_levels.clear();
    nx->get_rg_nested_cat(avail_levels, rule.v_cat, tree_pi);
    if(avail_levels.size() <= 1) avail_levels = tree_pi.cat_levels->at(rule.v_cat); // trivial split
    partition_levels(rule, avail_levels, tree_pi, gen);
    t.birth(nx_nid, rule);
  }
}

// [[Rcpp::export]]
void test_rg_nest(Rcpp::Nullable<Rcpp::List> cat_levels_list,
                  Rcpp::Nullable<Rcpp::List> edge_mat_list,
                  Rcpp::Nullable<Rcpp::List> nest_list,
                  Rcpp::IntegerMatrix cov_ensm,
                  int p_cont, int p_cat, int nest_v_option)
{
  Rcpp::RNGScope scope;
  RNG gen;
  
  int R = cov_ensm.cols();
  std::vector<std::set<int>> cat_levels;
  std::vector<std::vector<edge>> edges;
  
  parse_cat_levels(cat_levels, p_cat, cat_levels_list);
  parse_graphs(edges, p_cat, edge_mat_list);
  
  std::vector<hi_lo_map> nesting;
  std::vector<edge_map> nest_graph_in;
  std::vector<edge_map> nest_graph_out;
  std::vector<std::map<int, std::set<int>>> nest_graph_components;
  parse_nesting(nesting, nest_graph_in, nest_graph_out, nest_graph_components, p_cont, cov_ensm, cat_levels, nest_list);
  
  data_info di;
  di.R = R;
  di.p_cont = p_cont;
  di.p_cat = p_cat;
  di.p = p_cont + p_cat;
  
  std::vector<tree_prior_info> tree_pi_vec(R);
  for(int r = 0; r < R; ++r){
    tree_pi_vec[r].cat_levels = &cat_levels;
    tree_pi_vec[r].edges = &edges;
    tree_pi_vec[r].nesting = &nesting;
    tree_pi_vec[r].nest_v_option = nest_v_option;
    tree_pi_vec[r].nest_in = &(nest_graph_in[r]);
    tree_pi_vec[r].nest_out = &(nest_graph_out[r]);
    tree_pi_vec[r].nest_components = &(nest_graph_components[r]);
  }
  
  tree t;
  int r = 0;
  create_tree(r, t, tree_pi_vec[r], gen);
  t.print();
  
  
  // now go to each leaf and come up with nest_theta
  tree::npv bnv;
  t.get_bots(bnv);
  
  for(tree::npv_it l_it = bnv.begin(); l_it != bnv.end(); ++l_it){
    int nid = (*l_it)->get_nid();
    Rcpp::Rcout << "Leaf " << nid;
    
    /*
    Rcpp::Rcout << " ancestors split on categorical variables:";
    std::vector<int> anc_v;
    (*l_it)->get_anc_v_cat(anc_v);
    for(std::vector<int>::iterator it = anc_v.begin(); it != anc_v.end(); ++it) Rcpp::Rcout << " " << *it;
    Rcpp::Rcout << std::endl;
    */
    
    std::vector<double> nest_theta(di.p, 0.0);
    compute_nested_theta(nest_theta, t, nid, p_cont, p_cat, tree_pi_vec[0]);
    for(std::vector<double>::iterator it = nest_theta.begin(); it != nest_theta.end(); ++it) Rcpp::Rcout << " " << *it;
    Rcpp::Rcout << std::endl;
    
  }
  
  
  /*
  // what levels of classroom are available at node 2?
  int nx_nid = 4;
  tree::tree_p nx = t.get_ptr(nx_nid);
  int v_cat = 0;
  std::set<int> avail_levels;
  nx->get_rg_nested_cat(avail_levels, v_cat, tree_pi_vec[r]);
  
  Rcpp::Rcout << "avail levels for " << v_cat << " at " << nx_nid <<":";
  for(std::set<int>::iterator it = avail_levels.begin(); it != avail_levels.end(); ++it) Rcpp::Rcout << " " << *it;
  Rcpp::Rcout << std::endl;
  
  v_cat = 1;
  avail_levels.clear();
  nx->get_rg_nested_cat(avail_levels, v_cat, tree_pi_vec[r]);
  Rcpp::Rcout << "avail levels for " << v_cat << " at " << nx_nid <<":";
  for(std::set<int>::iterator it = avail_levels.begin(); it != avail_levels.end(); ++it) Rcpp::Rcout << " " << *it;
  Rcpp::Rcout << std::endl;
  
  v_cat = 2;
  avail_levels.clear();
  nx->get_rg_nested_cat(avail_levels, v_cat, tree_pi_vec[r]);
  Rcpp::Rcout << "avail levels for " << v_cat << " at " << nx_nid <<":";
  for(std::set<int>::iterator it = avail_levels.begin(); it != avail_levels.end(); ++it) Rcpp::Rcout << " " << *it;
  Rcpp::Rcout << std::endl;
  */
}
