#ifndef GUARD_update_tree_gen_h
#define GUARD_update_tree_gen_h

#include "update_tree.h"
#include "gen_models.h"

void compute_laplace_approx_single(std::map<int, laplace_approx> &lap_map, suff_stat &ss, data_info &di, tree_prior_info &tree_pi, GenModel &gmp);
void compute_laplace_approx_multi(std::map<int, laplace_approx> &lap_map, suff_stat &ss, int &r, data_info &di, tree_prior_info &tree_pi, GenModel &gmp);
void draw_mu_gen(tree &t, std::map<int, laplace_approx> &lap_map, RNG &gen);
void compute_ss_grow_gen_single(suff_stat &ss, std::map<int, laplace_approx> &lap_map, int &nx_nid, rule_t &rule, data_info &di, tree_prior_info &tree_pi, GenModel &gmp);
void compute_ss_grow_gen_multi(suff_stat &ss, std::map<int, laplace_approx> &lap_map, int &nx_nid, rule_t &rule, int &r, data_info &di, tree_prior_info &tree_pi, GenModel &gmp);

void compute_ss_prune_gen_single(suff_stat &ss, std::map<int, laplace_approx> &lap_map, int &nxl_nid, int &nxr_nid, int &nx_nid, data_info &di, tree_prior_info &tree_pi, GenModel &gmp);
void compute_ss_prune_gen_multi(suff_stat &ss, std::map<int, laplace_approx> &lap_map, int &nxl_nid, int &nxr_nid, int &nx_nid, int &r, data_info &di, tree_prior_info &tree_pi, GenModel &gmp);

void grow_tree_gen_single(tree &t, suff_stat &ss_train, suff_stat &ss_test, std::map<int, laplace_approx> &lap_map, int &accept, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, GenModel &gmp, RNG &gen);
void grow_tree_gen_multi(tree &t, suff_stat &ss_train, suff_stat &ss_test, std::map<int, laplace_approx> &lap_map, int &accept, int &r, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, GenModel &gmp, RNG &gen);

void prune_tree_gen_single(tree &t, suff_stat &ss_train, suff_stat &ss_test, std::map<int, laplace_approx> &lap_map, int &accept, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, GenModel &gmp, RNG &gen);
void prune_tree_gen_multi(tree &t, suff_stat &ss_train, suff_stat &ss_test, std::map<int, laplace_approx> &lap_map, int &accept, int &r, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, GenModel &gmp, RNG &gen);

void update_tree_gen_single(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, GenModel &gmp, RNG &gen);
void update_tree_gen_multi(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, int &r, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, GenModel &gmp, RNG &gen);

#endif /* update_tree_gen_h */