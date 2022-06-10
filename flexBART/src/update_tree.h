#ifndef GUARD_update_tree_h
#define GUARD_update_tree_h

#include "funs.h"

void grow_tree(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);
void prune_tree(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);
void update_tree(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen);
#endif /* update_tree_h */
