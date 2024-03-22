

#ifndef GUARD_funs_h
#define GUARD_funs_h

#include "tree.h"


void tree_traversal(suff_stat &ss, tree &t, data_info &di);


//void fit_single_tree(double* ftemp, tree &t, data_info &di);
void fit_ensemble(std::vector<double> &fit, std::vector<tree> &t_vec, data_info &di);


std::string write_tree(tree &t, tree_prior_info &tree_pi, set_str_conversion &set_str);
void read_tree(tree &t, std::string &tree_string, set_str_conversion &set_str);

void update_theta_u(std::vector<double> &theta, double &u, std::vector<int> &var_count, int &p, double &a_u, double &b_u, RNG &gen);


#endif /* funs_h */
