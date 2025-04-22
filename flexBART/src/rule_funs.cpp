
#include "rule_funs.h"

void draw_rule(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen){
  rule.clear();
  int v_raw = gen.multinomial(di.p, tree_pi.theta);
  if(v_raw < di.p_cont){
    rule.is_cat = false;
    rule.v_aa = v_raw;
    draw_aa_rule(rule, t, nid, di, tree_pi, gen);
  } else{
    rule.is_cat = true;
    rule.v_cat = v_raw - di.p_cont;
    draw_cat_rule(rule, t, nid, di, tree_pi, gen);
  }
}

void draw_cat_rule(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen){
  tree::tree_p nx = t.get_ptr(nid); // at what node are we proposing this rule.
  int rule_counter = 0;
  std::set<int> avail_levels = tree_pi.cat_levels->at(rule.v_cat); // get the full set of levels for this variable
  nx->get_rg_cat(rule.v_cat, avail_levels); // determine the set of levels available at nx.
  // if there is only one level left for this variable at nx, we will just propose a trivial split
  // and will reset the value of avail_levels to be the full set of all levels for the variable
  if(avail_levels.size() <= 1) avail_levels = tree_pi.cat_levels->at(rule.v_cat);
  
  rule.l_vals.clear();
  rule.r_vals.clear();
  
  if(tree_pi.edges->at(rule.v_cat).size() > 0){
    // if we explicitly use the graph to split the variables
    graph_partition(rule.l_vals, rule.r_vals, tree_pi.edges->at(rule.v_cat), avail_levels, tree_pi.graph_cut_type, gen);
  } else{
    // otherwise we default to splitting the available levels uniformly at random: prob 0.5 to go to each child
    rule_counter = 0;
    while( ((rule.l_vals.size() == 0) || (rule.r_vals.size() == 0)) && rule_counter < 1000 ){
      rule.l_vals.clear();
      rule.r_vals.clear();
      for(set_it it = avail_levels.begin(); it != avail_levels.end(); ++it){
        if(gen.uniform() <= 0.5) rule.l_vals.insert(*it);
        else rule.r_vals.insert(*it);
      }
      ++(rule_counter);
    }
    if(rule_counter == 1000){
      Rcpp::stop("[draw_cat_rule]: failed to generate valid categorical split in 1000 attempts"); // this should almost surely not get triggered.
    }
  }
  if( (rule.l_vals.size() == 0) || (rule.r_vals.size() == 0) ){
    Rcpp::stop("[draw_cat_rule]: proposed an invalid categorical rule!");
  }
}

void draw_aa_rule(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen)
{
  double c_upper = 1.0;
  double c_lower = -1.0;
  tree::tree_p nx = t.get_ptr(nid); // at what node are we proposing this rule.
  if(tree_pi.cutpoints->at(rule.v_aa).size() > 0){
    // draw the cutpoint from the supplied cutpoints
    c_lower = *(tree_pi.cutpoints->at(rule.v_aa).begin()); // returns smallest element in set
    c_upper = *(tree_pi.cutpoints->at(rule.v_aa).rbegin()); // reverse iterator, returns largest value in set
    nx->get_rg_aa(rule.v_aa, c_lower, c_upper);
    if(c_lower >= c_upper){
      // this is a weird tree and we'll just propose a trivial split
      c_lower = *(tree_pi.cutpoints->at(rule.v_aa).begin());
      c_upper = *(tree_pi.cutpoints->at(rule.v_aa).rbegin());
    }
    std::vector<double> valid_cutpoints;
    if(tree_pi.cutpoints->at(rule.v_aa).count(c_lower) != 1 || tree_pi.cutpoints->at(rule.v_aa).count(c_upper) != 1){
      // c_lower and c_upper were not found in the set of available cutpoints
      Rcpp::Rcout << "[draw_rule]: attempting to select a cutpoint from given set" << std::endl;
      Rcpp::Rcout << "  lower bound is: " << c_lower << " count in set is " << tree_pi.cutpoints->at(rule.v_aa).count(c_lower) << std::endl;
      Rcpp::Rcout << "  upper bound is: " << c_upper << " count in set is " << tree_pi.cutpoints->at(rule.v_aa).count(c_upper) << std::endl;
      //Rcpp::Rcout << "  cutpoints are:";
      //for(std::set<double>::iterator it = tree_pi.cutpoints->at(rule.v_aa).begin(); it != tree_pi.cutpoints->at(rule.v_aa).end(); ++it) Rcpp::Rcout << " " << *it;
      //Rcpp::Rcout << std::endl;
      Rcpp::stop("we should never have a c that is outside the pre-defined set of cutpoints!");
    }
    // we want to draw from the cutpoints exclusive of c_lower & c_upper;
    // i.e. we want to start with the one just after c_lower and just before c_upper
    // std::set::lower_bound: iterator at first element that is not considered to come before
    // std::set::upper_bound: iterator at first element considered to come after
    // if value is not in set, lower_bound and upper_bound give same result
    // if value is in set: lower bound returns the value, upper bound returns the next value
    for(std::set<double>::iterator it = tree_pi.cutpoints->at(rule.v_aa).upper_bound(c_lower); it != tree_pi.cutpoints->at(rule.v_aa).lower_bound(c_upper); ++it){
      valid_cutpoints.push_back(*it);
    }
    int num_cutpoints = valid_cutpoints.size();
    if(num_cutpoints < 1){
      // no valid splits are available; we will just pick something, all of the observations will go to one child anyway...
      valid_cutpoints.clear();
      for(std::set<double>::iterator it = tree_pi.cutpoints->at(rule.v_aa).begin(); it != tree_pi.cutpoints->at(rule.v_aa).end(); ++it){
        valid_cutpoints.push_back(*it);
      }
      num_cutpoints = valid_cutpoints.size();
    }
    // at this point, valid cutpoints is a vector containing the available cutpoints at this node. we pick one uniformly.
    rule.c = valid_cutpoints[floor(gen.uniform() * num_cutpoints)];
    
  } else{
    // draw cutpoints uniformly
    nx->get_rg_aa(rule.v_aa, c_lower, c_upper);
    if(c_lower >= c_upper){
      c_lower = -1.0;
      c_upper = 1.0;
    }
    rule.c = gen.uniform(c_lower, c_upper);
  }
}
