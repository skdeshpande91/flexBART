

#ifndef GUARD_funs_h
#define GUARD_funs_h

#include "tree.h"


void tree_traversal(suff_stat &ss, tree &t, data_info &di);


//void fit_single_tree(double* ftemp, tree &t, data_info &di);
void fit_ensemble(std::vector<double> &fit, std::vector<tree> &t_vec, data_info &di);


void compute_suff_stat_grow(suff_stat &orig_suff_stat, suff_stat &new_suff_stat, int &nx_nid, rule_t &rule, tree &t, data_info &di);
void compute_suff_stat_prune(suff_stat &orig_suff_stat, suff_stat &new_suff_stat, int &nl_nid, int &nr_nid, int &np_nid, tree &t, data_info &di);

double compute_lil(suff_stat &ss, int &nid, double &sigma, data_info &di, tree_prior_info &tree_pi);
void draw_mu(tree &t, suff_stat &ss, double &sigma, data_info &di, tree_prior_info &tree_pi, RNG &gen);

std::string write_tree(tree &t, data_info &di, set_str_conversion &set_str);
void read_tree(tree &t, std::string &tree_string, data_info &di, set_str_conversion &set_str);

void draw_rule(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen);


void build_symmetric_edge_map(edge_map &emap, std::vector<edge> &edges);
std::vector<edge> get_induced_edges(std::vector<edge> &edges, std::set<int> &vertex_subset);
void dfs(int v, std::map<int, bool> &visited, std::vector<int> &comp, edge_map &emap);
void find_components(std::vector<std::vector<int> > &components, std::vector<edge> &edges, std::set<int> &vertices);
void get_unique_edges(std::vector<edge> &edges);

std::vector<edge> boruvka(std::vector<edge> &edges, std::set<int> &vertices);
std::vector<edge> wilson(std::vector<edge> &edges, std::set<int> &vertices, RNG &gen);

void graph_partition(std::set<int> &avail_levels, std::set<int> &l_vals, std::set<int> &r_vals, std::vector<edge> &edges, int &K, bool &use_wilson, int &cut_type, RNG &gen);
/*
void dfs(int i, std::vector<bool> &visited, std::vector<int> &comp, int &n, arma::mat &A);
void find_components(std::vector<std::vector<int> > &components, arma::mat &A);
std::pair<int,int> find_min_edge_weight(std::vector<int> &components, int &n, arma::mat &W);
arma::mat boruvka(arma::mat &W);
//void get_edge_probs(std::vector<double> &cut_ix_probs, const arma::mat &cut_A, const arma::uvec &mst_index, const int &n);
arma::uword get_cut_edge(const arma::mat &cut_A, const arma::mat &cut_W, const arma::uvec mst_edge_index, const int &mst_cut_type, RNG &gen);
void graph_partition(std::set<int> &vals, std::set<int> &l_vals, std::set<int> &r_vals, std::vector<unsigned int> &adj_support, int &K, int &mst_cut_type, RNG &gen);
*/
void update_theta_u(std::vector<double> &theta, double &u, std::vector<int> &var_count, int &p, double &a_u, double &b_u, RNG &gen);
void update_theta_rc(double& theta_rc, int &rc_var_count, int &rc_rule_count, double &a_rc, double &b_rc, int &p_cont, RNG &gen);
//void update_theta_cont(std::vector<double> &theta_cont, std::vector<int> &cont_var_count, int &cont_rule_count, double &a_cont, double &b_cont, int &p_cont, RNG &gen);
//void update_theta_u_cat(std::vector<double> &theta_cat, std::vector<int> &cat_var_count, double &u_cat, double& a_cat, double& b_cat, int &p_cat, RNG &gen);

#endif /* funs_h */
