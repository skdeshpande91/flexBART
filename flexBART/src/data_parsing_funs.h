//
//  data_parsing_funs.hpp
//  
//
//  Created by Sameer Deshpande on 11/19/23.
//

#ifndef data_parsing_funs_h
#define data_parsing_funs_h
#include "structs.h"

void parse_cutpoints(std::vector<std::set<double>> &cutpoints, int p_cont, Rcpp::Nullable<Rcpp::List> &cutpoints_list);
void parse_cat_levels(std::vector<std::set<int>> &cat_levels, int &p_cat, Rcpp::Nullable<Rcpp::List> &cat_levels_list);
void parse_nesting(std::vector<hi_lo_map> &nesting, std::vector<edge_map> &nest_graph_in,std::vector<edge_map> &nest_graph_out,
                   int &p_cont, Rcpp::IntegerMatrix &cov_ensm, std::vector<std::set<int>> &cat_levels, Rcpp::Nullable<Rcpp::List> &nest_list);

//void parse_cutpoints(std::vector<std::set<double>> &cutpoints, int p_cont, Rcpp::List &tmp_cutpoints, Rcpp::LogicalVector &unif_cuts);
//void parse_cat_levels(std::vector<std::set<int>> &cat_levels, std::vector<int> &K, int &p_cat, Rcpp::List &tmp_cat_levels);
//void parse_training_data(int &n_train, int &p_cont, int &p_cat, Rcpp::NumericMatrix &tX_cont_train, Rcpp::IntegerMatrix &tX_cat_train);
//void parse_testing_data(int &n_test, Rcpp::NumericMatrix &tX_cont_test, Rcpp::IntegerMatrix &tX_cat_test, const int &p_cat, const int &p_cont);

#endif /* data_parsing_funs_hpp */
