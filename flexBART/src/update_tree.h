#ifndef GUARD_update_tree_h
#define GUARD_update_tree_h

#include "rule_funs.h"

//void compute_ss_grow(suff_stat &ss, int &nx_nid, rule_t &rule, data_info &di);
//void compute_ss_prune(suff_stat &ss, int &nxl_nid, int &nxr_nid, int& nx_nid, rule_t &rule, data_info &di);
void compute_suff_stat_grow(suff_stat &orig_suff_stat, suff_stat &new_suff_stat, int &nx_nid, rule_t &rule, data_info &di);
void compute_suff_stat_prune(suff_stat &orig_suff_stat, suff_stat &new_suff_stat, int &nl_nid, int &nr_nid, int &np_nid, data_info &di);

double compute_lil(suff_stat &ss, int &nid, int &r, double &sigma, data_info &di, tree_prior_info &tree_pi);

void draw_mu(tree &t, suff_stat &ss, int &r, double &sigma, data_info &di, tree_prior_info &tree_pi, RNG &gen);


void grow_tree_unnested(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, int &r, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);
void grow_tree_nested(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, int &r, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);

void prune_tree(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, int &r, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);
void update_tree(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, int &r, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);

void update_tree_unnested(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, int &r, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);
void update_tree_nested(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, int &r, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);

#endif /* update_tree_h */
