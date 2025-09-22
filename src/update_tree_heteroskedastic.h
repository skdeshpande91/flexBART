#ifndef GUARD_update_tree_heteroskedastic_h
#define GUARD_update_tree_heteroskedastic_h

#include "update_tree.h"

void compute_jump_posterior_heteroskedastic_single(std::map<int, jump_post> &jp_map, suff_stat &ss, data_info &di, tree_prior_info &tree_pi);
void compute_jump_posterior_heteroskedastic_multi(std::map<int, jump_post> &jp_map, suff_stat &ss, int &r, data_info &di, tree_prior_info &tree_pi);

void compute_ss_grow_heteroskedastic_single(suff_stat &ss, std::map<int, jump_post> &jp_map, int &nx_nid, rule_t &rule, data_info &di, tree_prior_info &tree_pi);
void compute_ss_grow_heteroskedaskic_multi(suff_stat &ss, std::map<int, jump_post> &jp_map, int &nx_nid, rule_t &rule, int &r, data_info &di, tree_prior_info &tree_pi);

void compute_ss_prune_heteroskedastic_single(suff_stat &ss, std::map<int, jump_post> &jp_map, int &nxl_nid, int &nxr_nid, int &nx_nid, data_info &di, tree_prior_info &tree_pi);
void compute_ss_prune_heteroskedastic_multi(suff_stat &ss, std::map<int, jump_post> &jp_map, int &nxl_nid, int &nxr_nid, int &nx_nid, int &r, data_info &di, tree_prior_info &tree_pi);

void grow_tree_heteroskedastic_single(tree &t, suff_stat &ss_train, suff_stat &ss_test, std::map<int, jump_post> &jp_map, int &accept, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);
void grow_tree_heteroskedastic_multi(tree &t, suff_stat &ss_train, suff_stat &ss_test, std::map<int, jump_post> &jp_map, int &accept, int &r, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);

void prune_tree_heteroskedastic_single(tree &t, suff_stat &ss_train, suff_stat &ss_test, std::map<int, jump_post> &jp_map, int &accept, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);
void prune_tree_heteroskedastic_multi(tree &t, suff_stat &ss_train, suff_stat &ss_test, std::map<int, jump_post> &jp_map, int &accept, int &r, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);

void update_tree_heteroskedastic_single(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);
void update_tree_heteroskedastic_multi(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, int &r, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);

#endif /* update_tree_heteroskedastic_h */