#include "../flexBART/src/data_parsing_funs.h"
#include "../flexBART/src/rule_funs.h"

// [[Rcpp::export]]
void test_nest_graph(Rcpp::Nullable<Rcpp::List> cat_levels_list,
                     Rcpp::Nullable<Rcpp::List> nest_list,
                     Rcpp::IntegerMatrix cov_ensm,
                     int p_cont, int p_cat)
{
  int R = cov_ensm.cols();
  
  std::vector<std::set<int>> cat_levels;
  parse_cat_levels(cat_levels, p_cat, cat_levels_list);
  
  std::vector<hi_lo_map> nesting;
  std::vector<edge_map> nest_graph_in;
  std::vector<edge_map> nest_graph_out;
  std::vector<std::map<int, std::set<int>>> nest_graph_components;
  parse_nesting(nesting, nest_graph_in, nest_graph_out, nest_graph_components, p_cont, cov_ensm, cat_levels, nest_list);
  
  for(int r = 0; r < R; ++r){
    Rcpp::Rcout << "Ensemble " << r << " nest graph has components: " << std::endl;
    for(std::map<int, std::set<int>>::iterator c_it = nest_graph_components[r].begin(); c_it != nest_graph_components[r].end(); ++c_it){
      Rcpp::Rcout << "    Component " << c_it->first << ":";
      for(std::set<int>::iterator it = c_it->second.begin(); it != c_it->second.end(); ++it) Rcpp::Rcout << " " << *it;
      Rcpp::Rcout << std::endl;
    }
    Rcpp::Rcout << "    In edges:" << std::endl;
    for(edge_map_it e_it = nest_graph_in[r].begin(); e_it != nest_graph_in[r].end(); ++e_it){
      Rcpp::Rcout << "      " << e_it->first << ":";
      for(std::vector<edge>::iterator it = e_it->second.begin(); it != e_it->second.end(); ++it){
        Rcpp::Rcout << " " << it->source << "-" << it->sink;
      }
      Rcpp::Rcout << std::endl;
    }
    Rcpp::Rcout << "    Out edges:" << std::endl;
    for(edge_map_it e_it = nest_graph_out[r].begin(); e_it != nest_graph_out[r].end(); ++e_it){
      Rcpp::Rcout << "      " << e_it->first << ":";
      for(std::vector<edge>::iterator it = e_it->second.begin(); it != e_it->second.end(); ++it){
        Rcpp::Rcout << " " << it->source << "-" << it->sink;
      }
      Rcpp::Rcout << std::endl;
    }
  }
}
