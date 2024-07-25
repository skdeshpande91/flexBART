//
//  data_parsing_funs.cpp
//  
//
//  Created by Sameer Deshpande on 11/19/23.
//

#include "data_parsing_funs.h"

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
