//
//  data_parsing_funs.cpp
//  
//
//  Created by Sameer Deshpande on 11/19/23.
//

#include "data_parsing_funs.h"


void parse_cutpoints(std::vector<std::set<double>> &cutpoints, int p_cont, Rcpp::Nullable<Rcpp::List> &cutpoints_list)
{
  cutpoints.clear();
  cutpoints.resize(p_cont, std::set<double>());
  if(cutpoints_list.isNotNull()){
    // cutpoints provided
    Rcpp::List tmp_cutpoints = Rcpp::List(cutpoints_list);
    if(tmp_cutpoints.size() == p_cont){
      for(int j = 0; j < p_cont; ++j){
        if(tmp_cutpoints[j] != R_NilValue){
          Rcpp::NumericVector cutpoints_vec = Rcpp::as<Rcpp::NumericVector>(tmp_cutpoints[j]);
          if(cutpoints_vec.size() <= 1){
            Rcpp::Rcout << "Only " << cutpoints_vec.size() << " cutpoints supplied for variable X_cont[," << j+1 << "]" << std::endl;
            Rcpp::stop("[parse_cutpoints]: Not enough cutpoints supplied!");
          } else{
            for(int l = 0; l < cutpoints_vec.size(); l++) cutpoints[j].insert(cutpoints_vec[l]); // save the cutpoints
          } // closes if/else checking that more than 2 cutpoints are provided for the variable
        } // closes if checking that cutpoints are provided for this variable
      } // closes loop over all continuous variables.
    } else{
      Rcpp::Rcout << "p_cont = " << p_cont;
      Rcpp::Rcout << "  cutpoints_list.size() = " << tmp_cutpoints.size() << std::endl;
      Rcpp::stop("cutpoints_list needs to have length p_cont!");
    } // closes if/else checking that the supplied cutpoints_list is the right size
  } // closes if checking that cutpoints were provided
}

/*
 // APRIL 22: stage for deprecation.
void parse_cutpoints(std::vector<std::set<double>> &cutpoints, int p_cont, Rcpp::List &tmp_cutpoints, Rcpp::LogicalVector &unif_cuts)
{
  cutpoints.clear();
  cutpoints.resize(p_cont, std::set<double>());
  if(tmp_cutpoints.size() == p_cont && unif_cuts.size() == p_cont){
    for(int j = 0; j < p_cont; j++){
      
      if(unif_cuts[j] == 0){
        Rcpp::NumericVector cutpoints_vec = Rcpp::as<Rcpp::NumericVector>(tmp_cutpoints[j]);
        if(cutpoints_vec.size() <= 1){
          Rcpp::Rcout << "Only " << cutpoints_vec.size() << " cutpoints supplied for variable X_cont[," << j+1 << "]" << std::endl;
          Rcpp::stop("[parse_cutpoints]: Not enough cutpoints supplied!");
        } else{
          for(int l = 0; l < cutpoints_vec.size(); l++) cutpoints[j].insert(cutpoints_vec[l]);
        }
      }
    }
  } else{
    Rcpp::Rcout << "p_cont = " << p_cont;
    Rcpp::Rcout << "  cutpoints_list.size() = " << tmp_cutpoints.size() << std::endl;
    Rcpp::Rcout << "  unif_cuts.size() = " << unif_cuts.size() << std::endl;
    Rcpp::stop("cutpoints_list & unif_cuts needs to have length p_cont!");
  }
}
*/


void parse_cat_levels(std::vector<std::set<int>> &cat_levels, int &p_cat, Rcpp::Nullable<Rcpp::List> &cat_levels_list)
{
  cat_levels.clear();
  cat_levels.resize(p_cat, std::set<int>());
  
  if(cat_levels_list.isNotNull()){
    Rcpp::List tmp_cat_levels = Rcpp::List(cat_levels_list);
    if(tmp_cat_levels.size() == p_cat){
      for(int j = 0; j < p_cat; ++j){
        if(tmp_cat_levels[j] != R_NilValue){
          Rcpp::IntegerVector levels_vec = Rcpp::as<Rcpp::IntegerVector>(tmp_cat_levels[j]);
          for(int l = 0; l < levels_vec.size(); l++) cat_levels[j].insert(levels_vec[l]); // save categorical levels
        } else{
          Rcpp::Rcout << "[parse_cat_levels]: j = " << j << " no levels provided" << std::endl;
          Rcpp::stop("cat_levels_list cannot contain null elements");
        }
      }
    } else{
      Rcpp::Rcout << "[parse_cat_levels]: p_cat = " << p_cat;
      Rcpp::Rcout << " cat_levels_list.size() = " << tmp_cat_levels.size() << std::endl;
      Rcpp::stop("cat_levels_list must have size equal to p_cat!");
    } // closes if/else checking that supplied cat_levels_list is the right size
  } // closes if checking that cat_levels are provided
}

/*
 // APRIL 22 DEPRECATE
void parse_cat_levels(std::vector<std::set<int>> &cat_levels, std::vector<int> &K, int &p_cat, Rcpp::List &tmp_cat_levels)
{
  cat_levels.clear();
  cat_levels.resize(p_cat, std::set<int>());
  K.clear();
  K.resize(p_cat);
  if(tmp_cat_levels.size() == p_cat){
    for(int j = 0; j < p_cat; j++){
      Rcpp::IntegerVector levels_vec = Rcpp::as<Rcpp::IntegerVector>(tmp_cat_levels[j]);
      for(int l = 0; l < levels_vec.size(); l++) cat_levels[j].insert(levels_vec[l]);
      K[j] = levels_vec.size();
    }
  } else{
    Rcpp::Rcout << "p_cat = " << p_cat;
    Rcpp::Rcout << "cat_levels_list.size() = " << tmp_cat_levels.size();
    Rcpp::stop("cat_levels_list must have size equal to p_cat!");
  }
}
*/

void parse_nesting(std::vector<hi_lo_map> &nesting,
                   std::set<int> &nest_graph_vertices,
                   std::vector<edge> &nest_graph_edges,
                   int &p_cont,
                   int &p,
                   std::vector<std::set<int>> &cat_levels,
                   Rcpp::Nullable<Rcpp::List> &nest_list)
{
  nesting.clear();
  nest_graph_edges.clear();
  
  // create a vertex for every covariate
  for(int j = 0; j < p; ++j) nest_graph_vertices.insert(j);
  
  if(nest_list.isNotNull()){
    Rcpp::List tmp_nest_list = Rcpp::List(nest_list);
    nesting.resize(tmp_nest_list.size(), hi_lo_map());
    for(int nix = 0; nix < tmp_nest_list.size(); ++nix){
      Rcpp::List tmp_list = tmp_nest_list[nix];
      Rcpp::IntegerVector tmp_pair = tmp_list[0];
      Rcpp::IntegerVector tmp_map = tmp_list[1];
      int hi = tmp_pair[0];
      int lo = tmp_pair[1];
      int K = cat_levels[hi-p_cont].size();
      if( (hi >= p || lo >= p) || (hi < p_cont || lo < p_cont)){
        Rcpp::Rcout << "Element " << nix+1 << " of nest_list:";
        Rcpp::Rcout << " hi-res variable = " << hi << " lo-res variable =" << lo;
        Rcpp::Rcout << " p = " << p << std::endl;
        Rcpp::stop("Variable indices passed in each element of nest_list should be in {p_cont, ..., p-1}!");
      } else if(tmp_map.size() != K){
        Rcpp::Rcout << "Element " << nix+1 << " of nest_list:";
        Rcpp::Rcout << " hi-res variable = " << hi << " lo-res variable =" << lo << std::endl;
        Rcpp::Rcout << " Passed mapping with " << tmp_map.size() << " hi-res values but expected " << K << std::endl;
        Rcpp::stop("Every hi-res value must be mapped to a lo-res value!");
      } else{
        // now we can continue
        nesting[nix].hi = hi;
        nesting[nix].lo = lo;
        std::map<int, std::set<int>>::iterator nest_it = nesting[nix].map.begin();
        for(int hi_val = 0; hi_val < tmp_map.size(); ++hi_val){
          int lo_val = tmp_map[hi_val];
          if(lo_val >= K){
            Rcpp::Rcout << "Element" << nix+1 << " of nest_list:";
            Rcpp::Rcout << " hi-res variable = " << hi << " lo-res variable =" << lo << std::endl;
            Rcpp::Rcout << " hi_res value = " << hi_val << " passed lo-res value = " << lo_val << std::endl;
            Rcpp::Rcout << " There are only " << K << " lo-res values..." << std::endl;
            Rcpp::stop("Passed an invalid lo-res value in this mapping!");
          } else{
            nest_it = nesting[nix].map.find(lo_val);
            if(nest_it == nesting[nix].map.end()){
              // we haven't seen this lo-res value before
              //Rcpp::Rcout << "  adding map element for " << lo_val << std::endl;
              nesting[nix].map.insert(std::pair<int,std::set<int>>(lo_val, std::set<int>()));
            }
            nesting[nix].map.find(lo_val)->second.insert(hi_val);
          } // closes if/else checking if lo-res value is valid
        } // closes loop over values of the hi-res variable
      } // closes if/else checking that current mapping is of right size
      // create an edge for the graph encoding nested edges
      nest_graph_edges.push_back(edge(hi,lo));
    } // closes loop over supplied nest_list
  } // closes if checking that nest_list is not null
}


/*
void parse_training_data(int &n_train, int &p_cont, int &p_cat, Rcpp::NumericMatrix &tX_cont_train, Rcpp::IntegerMatrix &tX_cat_train)
{
  if(tX_cont_train.size() > 1 && tX_cat_train.size() == 1){
    // only continuous predictors are available
    n_train = tX_cont_train.cols();
    p_cont = tX_cont_train.rows();
    p_cat = 0;
  } else if(tX_cont_train.size() == 1 && tX_cat_train.size() > 1){
    n_train = tX_cat_train.cols();
    p_cont = 0;
    p_cat = tX_cat_train.rows();
  } else if(tX_cont_train.size() > 1 && tX_cat_train.size() > 1){
    n_train = tX_cont_train.cols();
    if(tX_cat_train.cols() != n_train){
      Rcpp::Rcout << "X_cat_train has " << tX_cat_train.cols() << " rows but X_cont_train has " << n_train << " rows" << std::endl;
      Rcpp::stop("matrices for continuous and categorical inputs must have same number of rows!");
    }
    p_cont = tX_cont_train.rows();
    p_cat = tX_cat_train.rows();
  } else{
    Rcpp::stop("no covariate data provided!");
  }
}

void parse_testing_data(int &n_test, Rcpp::NumericMatrix &tX_cont_test, Rcpp::IntegerMatrix &tX_cat_test, const int &p_cat, const int &p_cont)
{
  if(tX_cont_test.size() > 1 && p_cont == 0) Rcpp::stop("No continuous preds in training data but continuous preds in testing data!");
  if(tX_cat_test.size() > 1 && p_cat == 0) Rcpp::stop("No categorical preds in training data but categorical preds in testing data!");
  
  if(tX_cont_test.size() > 1 && tX_cat_test.size() == 1) n_test = tX_cont_test.cols();
  else if(tX_cont_test.size() == 1 && tX_cat_test.size() > 1) n_test = tX_cat_test.cols();
  else if(tX_cont_test.size() == 1 && tX_cat_test.size() == 1) n_test = 0;
  else{
    n_test = tX_cont_test.cols();
    if(tX_cat_test.cols() != n_test){
      Rcpp::Rcout << "X_cont_test has " << tX_cont_test.cols() << " rows but X_cat_test has " << tX_cat_test.cols() << " rows" << std::endl;
      Rcpp::stop("matrices for continuous and categorical inputs must have same number of rows!");
    }
  }
}
*/
