#ifndef GUARD_rule_funs_h
#define GUARD_rule_funs_h

#include "tree.h"
#include "graph_funs.h"

void draw_aa_cutpoint(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen);
void partition_levels(rule_t &rule, std::set<int> &avail_levels, tree_prior_info &tree_pi, RNG &gen);
void compute_nested_theta(std::vector<double> &nest_theta, tree &t, int &nid, int &p_cont, int &p_cat, tree_prior_info &tree_pi);
void draw_rule(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen);
//void draw_rule_unnested(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen);
//void draw_rule_nested(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen);


//void draw_rule(rule_t & rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen);
/*
void draw_cat_rule(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen);
void draw_aa_rule(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen);
void draw_rule(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen);
*/

#endif
