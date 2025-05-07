#include "draw_tree.h"
#include "data_parsing_funs.h"
#include "funs.h"


// [[Rcpp::export("._draw_tree")]]
Rcpp::List drawTree(Rcpp::NumericMatrix tX_cont,
                    Rcpp::IntegerMatrix tX_cat,
                    Rcpp::IntegerMatrix cov_ensm,
                    Rcpp::Nullable<Rcpp::List> cutpoints_list,
                    Rcpp::Nullable<Rcpp::List> cat_levels_list,
                    Rcpp::Nullable<Rcpp::List> edge_mat_list,
                    Rcpp::Nullable<Rcpp::List> nest_list,
                    int graph_cut_type,
                    bool nest_v, int nest_v_option, bool nest_c,
                    double alpha, double beta,
                    int nd,
                    bool verbose, int print_every)
{
  Rcpp::RNGScope scope;
  RNG gen;
  set_str_conversion set_str; // for converting sets of integers into strings
  int p_cont = 0;
  int p_cat = 0;
  if(tX_cont.size() > 1) p_cont = tX_cont.rows();
  if(tX_cat.size() > 1) p_cat = tX_cat.rows();
  int p = p_cont + p_cat;
  
  int n = 0;
  if(p_cont > 0) n = tX_cont.cols();
  else if(p_cat > 0) n = tX_cat.cols();
  
  
  int R = 1;
  // BEGIN: set cutpoints & categorical levels + parse network structure
  std::vector<std::set<double>> cutpoints;
  if(p_cont > 0) parse_cutpoints(cutpoints, p_cont, cutpoints_list);
  
  std::vector<std::set<int>> cat_levels;
  std::vector<std::vector<edge>> edges;
  if(p_cat > 0){
    parse_cat_levels(cat_levels, p_cat, cat_levels_list);
    parse_graphs(edges, p_cat, edge_mat_list);
  }
  // END: set cutpoints & categorical levels + parse network structure
  
  // BEGIN: build graph encoding nesting relationships b/w categorical predictors
  std::vector<hi_lo_map> nesting;
  std::vector<edge_map> nest_graph_in;
  std::vector<edge_map> nest_graph_out;
  std::vector<std::map<int, std::set<int>>> nest_graph_components;
  parse_nesting(nesting, nest_graph_in, nest_graph_out, nest_graph_components, p_cont, cov_ensm, cat_levels, nest_list);
  // END: build graph encoding nesting relationships b/w categorical predictors

  // BEGIN: create splitting probabilities
  // declare stuff for variable selection
  std::vector<double> theta(p, 1.0/ (double) p);
  std::vector<int> var_count(p, 0); // count how many times a variable has been used in a splitting rule
  int rule_count = 0; // how many total decision rules are there in the ensemble
  // END: create splitting probabilities
  
  data_info di;
  di.n = n;
  di.p_cont = p_cont;
  di.p_cat = p_cat;
  di.p = p;
  di.R = R;
  if(p_cont > 0) di.x_cont = tX_cont.begin();
  if(p_cat > 0) di.x_cat = tX_cat.begin();
  
  // BEGIN: create tree prior info object
  tree_prior_info tree_pi;
  tree_pi.theta = &theta;
  tree_pi.var_count = &var_count;
  tree_pi.rule_count = &rule_count;
  tree_pi.nest_v = nest_v;
  tree_pi.nest_v_option = nest_v_option;
  tree_pi.nest_c = nest_c;
  
  if(p_cont > 0) tree_pi.cutpoints = &cutpoints;
  
  if(p_cat > 0){
    tree_pi.cat_levels = &cat_levels;
    tree_pi.edges = &edges;
    tree_pi.graph_cut_type = graph_cut_type;
  }
  tree_pi.nesting = &nesting;
  tree_pi.nest_in = &(nest_graph_in[0]);
  tree_pi.nest_out = &(nest_graph_out[0]);
  tree_pi.nest_components = &(nest_graph_components[0]);
  tree_pi.alpha = alpha;
  tree_pi.beta = beta;
  // END: create tree prior info object
  
  //
  Rcpp::IntegerVector num_clusters(nd);
  Rcpp::IntegerVector num_singletons(nd);
  Rcpp::IntegerVector num_empty(nd);
  Rcpp::IntegerVector max_cluster_size(nd);
  Rcpp::IntegerVector min_cluster_size(nd);
  arma::mat kernel = arma::zeros<arma::mat>(n,n); // kernel(i,ii) counts #times obs i & j in same leaf
    
  for(int iter = 0; iter < nd; ++iter){
    num_clusters[iter] = 0;
    num_singletons[iter] = 0;
    num_empty[iter] = 0;
    max_cluster_size[iter] = 0;
    min_cluster_size[iter] = 0;
  }
  
  tree t;
  suff_stat ss;
  
  for(int iter = 0; iter < nd; ++iter){
    if(verbose && iter % print_every == 0) Rcpp::Rcout << "  Drawing tree " << iter+1 << " of " << nd << std::endl;
    t.to_null();
    draw_tree(t, di, tree_pi, gen);
    ss.clear();
    tree_traversal(ss,t,di);
    num_clusters[iter] = ss.size();
    int singleton = 0;
    int empty = 0;
    int max_size = 0;
    int min_size = n;
    
    for(suff_stat_it ss_it = ss.begin(); ss_it != ss.end(); ++ss_it){
      
      int cluster_size = ss_it->second.size();
      if(cluster_size == 1) ++singleton;
      if(cluster_size == 0) ++empty;
      if(cluster_size > max_size) max_size = cluster_size;
      if(cluster_size < min_size) min_size = cluster_size;
      
      if(cluster_size > 1){
        for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
          for(int_it iit = it; iit != ss_it->second.end(); ++iit){
            if(*it != *iit){
              kernel(*it, *iit) += 1.0;
              kernel(*iit, *it) += 1.0;
            } else{
              kernel(*it, *iit) += 1.0;
            }
          } // closes second loop over obs in leaf
        } // closes loop over obs in leaf
      }
      num_singletons[iter] = singleton;
      num_empty[iter] = empty;
      max_cluster_size[iter] = max_size;
      min_cluster_size[iter] = min_size;
    } // closes loop over leafs
  }
  kernel /= (double) nd;
  
  Rcpp::List results;
  results["num_leafs"] = num_clusters;
  results["num_singletons"] = num_singletons;
  results["num_empty"] = num_empty;
  results["max_leaf_size"] = max_cluster_size;
  results["min_leaf_size"] = min_cluster_size;
  results["kernel"] = kernel;
  return results;
}
