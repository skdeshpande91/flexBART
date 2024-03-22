#ifndef GUARD_rule_funs_h
#define GUARD_rule_funs_h

#include "tree.h"
#include "graph_funs.h"
#include "polytope_funs.h"

void draw_cat_rule(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen);
//void draw_cont_rule(rule_t &rule, tree &t, int &nid, std::vector<double> &x0, data_info &di, tree_prior_info &tree_pi, RNG &gen);
void draw_aa_rule(rule-t &rule, tree &t, ind &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen);
void draw_rule(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen);


#endif
