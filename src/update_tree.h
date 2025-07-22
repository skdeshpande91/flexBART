#ifndef GUARD_update_tree_h
#define GUARD_update_tree_h

#include "rule_funs.h"


void compute_jump_posterior_single(std::map<int, jump_post> &jp_map, suff_stat &ss, double &sigma, data_info &di, tree_prior_info &tree_pi);
void compute_jump_posterior_multi(std::map<int, jump_post> &jp_map, suff_stat &ss, int &r, double &sigma, data_info &di, tree_prior_info &tree_pi);
double compute_lil(int &nid, std::map<int, jump_post> &jp_map);
void draw_mu(tree &t, std::map<int, jump_post> &jp_map, RNG &gen);
void compute_ss_grow_single(suff_stat &ss, std::map<int, jump_post> &jp_map, int &nx_nid, rule_t &rule, double &sigma, data_info &di, tree_prior_info &tree_pi);
void compute_ss_grow_multi(suff_stat &ss, std::map<int, jump_post> &jp_map, int &nx_nid, rule_t &rule, int &r, double &sigma, data_info &di, tree_prior_info &tree_pi);
void compute_ss_grow(suff_stat &ss, int &nx_nid, rule_t &rule, data_info &di);

void compute_ss_prune_single(suff_stat &ss, std::map<int, jump_post> &jp_map, int &nxl_nid, int &nxr_nid, int &nx_nid, double &sigma, data_info &di, tree_prior_info &tree_pi);
void compute_ss_prune_multi(suff_stat &ss, std::map<int, jump_post> &jp_map, int &nxl_nid, int &nxr_nid, int &nx_nid, int &r, double &sigma, data_info &di, tree_prior_info &tree_pi);
void compute_ss_prune(suff_stat &ss, int &nxl_nid, int &nxr_nid, int &nx_nid, data_info &di);

void grow_tree_single(tree &t, suff_stat &ss_train, suff_stat &ss_test, std::map<int, jump_post> &jp_map, int &accept, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);
void grow_tree_multi(tree &t, suff_stat &ss_train, suff_stat &ss_test, std::map<int, jump_post> &jp_map, int &accept, int &r, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);

void prune_tree_single(tree &t, suff_stat &ss_train, suff_stat &ss_test, std::map<int, jump_post> &jp_map, int &accept, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);
void prune_tree_multi(tree &t, suff_stat &ss_train, suff_stat &ss_test, std::map<int, jump_post> &jp_map, int &accept, int &r, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);

void update_tree_single(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);
void update_tree_multi(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, int &r, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);

#endif /* update_tree_h */
