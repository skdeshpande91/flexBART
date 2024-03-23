

#ifndef GUARD_graph_funs_h
#define GUARD_graph_funs_h

#include "structs.h"
#include "rng.h"
#include "tree.h"


// we will pass a list of Rcpp::NumericMatrices, which are computed using igraph::get_data_frame
// this function reads those matrices and builds a vector of edges
void parse_edge_mat(std::vector<edge> &edges, Rcpp::NumericMatrix &edge_mat, int &n_vertex);


// takes in the List of edge_mat's
void parse_graphs(std::vector<std::vector<edge>> &edges, int &p_cat, std::vector<int> &K, Rcpp::List &tmp_edge_mats, Rcpp::LogicalVector &graph_split);


void build_symmetric_edge_map(edge_map &emap, std::vector<edge> &edges, std::set<int> &vertices);
std::vector<edge> get_induced_edges(std::vector<edge> &edges, std::set<int> &vertex_subset);
void dfs(int v, std::map<int, bool> &visited, std::vector<int> &comp, edge_map &emap);
void find_components(std::vector<std::vector<int> > &components, std::vector<edge> &edges, std::set<int> &vertices);
void get_unique_edges(std::vector<edge> &edges);

arma::mat get_adjacency_matrix(edge_map emap);

arma::mat floydwarshall(std::vector<edge> &edges, std::set<int> &vertices);

void boruvka(std::vector<edge> &mst_edges, std::vector<edge> &edges, std::set<int> &vertices);
void wilson(std::vector<edge> &mst_edges, std::vector<edge> &edges, std::set<int> &vertices, RNG &gen);
void hotspot(std::set<int> &l_vals, std::set<int> &r_vals, std::vector<edge> &edges, std::set<int> &vertices, RNG &gen, bool debug = false);

void delete_unif_edge(std::set<int> &l_vals, std::set<int> &r_vals, std::vector<edge> &edges, std::set<int> &vertices, RNG &gen);
void signcheck_split(std::set<int> &l_vals, std::set<int> &r_vals, std::vector<edge> &edges, std::set<int> &vertices);

void graph_partition(std::set<int> &l_vals, std::set<int> &r_vals, std::vector<edge> &orig_edges, std::set<int> &avail_levels, int &cut_type, RNG &gen);


#endif /* graph_funs_h */
