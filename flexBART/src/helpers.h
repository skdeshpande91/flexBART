#ifndef GUARD_helper_funs_h
#define GUARD_helper_funs_h

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <cstddef>
#include <vector>


typedef std::vector<int>::iterator int_it; // iterator type for vectors of ints
typedef std::vector<double>::iterator dbl_it; // iterator type vectors of doubles
typedef std::set<int>::iterator set_it; // iterator type for sets of integers
typedef std::map<int, double>::const_iterator rc_it; // iterator type for random combination
typedef std::map<int, std::vector<int>> suff_stat; // sufficient statistic map: key is node id, value is vector of observations that land in that node
typedef suff_stat::iterator suff_stat_it; // iterator for sufficient statistc map

class rule_t{
public:
  bool is_aa; // is it an axis-aligned split
  bool is_cat; // is it a categorical split

  int v_aa;
  std::map<int, double> rc_weight; // weights of random combination
  double c; // cutpoint
  
  int v_cat; // index of variable of the categorical variable on which we split (always between 0 and p_cat)
  std::set<int> l_vals; // holds unique values of levels of v associated w/ left child
  std::set<int> r_vals; // holds unique values of levels of v associated w/ right child
  
  rule_t(){
    is_aa = true;
    is_cat = false;
    v_aa = 0;
    rc_weight = std::map<int,double>();
    c = 0.0;
    v_cat = 0;
    l_vals = std::set<int>();
    r_vals = std::set<int>();
  }
  void clear(){
    is_aa = false;
    is_cat = false;
    v_aa = 0;
    rc_weight.clear();
    c = 0.0;
    v_cat = 0;
    l_vals.clear();
    r_vals.clear();
  }
};

//structures for graph partitioning
struct edge{
  int source;
  int sink;
  double weight;
  edge(){source = 0; sink = 0; weight = 1.0;}
  edge(int source_in, int sink_in){source = source_in; sink = sink_in; weight = 1.0;}
  edge(int source_in, int sink_in, double weight_in){source = source_in; sink = sink_in; weight = weight_in;}
  void copy_edge(edge edge_in){source = edge_in.source; sink = edge_in.sink; weight = edge_in.weight;} // unlikely to ever use this
};

typedef std::map<int, std::vector<edge> > edge_map;
typedef std::vector<edge>::iterator edge_vec_it;
typedef std::map<int, std::vector<edge>>::iterator edge_map_it;


// class holding data dimensions and pointers to the covariate data
class data_info{
public:
  int n; // number of observations
  int p_cont; // number of continuous predictors
  int p_cat; // number of categorical predictors
  int p; // total number of predictors (likely will never every use this)
  double* x_cont; // pointer to the matrix of continuous predictors
  int* x_cat; // pointer to matrix of categorical predictors (levels coded as integers, beginning with 0)
  double* rp; // partial residual;
  data_info(){n = 0; p_cont = 0; p_cat = 0; p = 0; x_cont = 0; x_cat = 0;rp = 0;}
};



// holds hyperparameters for regression tree prior
class tree_prior_info{
public:
  double alpha; // 1st parameter of the branching process prior
  double beta; // 2nd parameter in branching process prior
  
  double prob_bd; // prob of proposing a grow (birth) or prune (death) move. almost always set to 1
  double prob_b; // prob of proposing a grow (birth) move. almost always set to 0.5
  
  std::vector<double>* theta; // prob. that we pick one variable out of p_cont + p_cat
  std::vector<int> *var_count; // counts how many times we split on a single variable
  int* rule_count; // how many total rules are there in the ensemble

  // unif_cuts passed an Rcpp::LogicalVector
  int* unif_cuts; // unif_cuts = 1 to draw cuts uniformly, = 0 to draw from pre-specified cutpoints, = minimum integer for NA
  std::vector<std::set<double> >* cutpoints;
  
  std::vector<std::set<int>> *cat_levels; // holds the levels of the categorical variables
  std::vector<std::vector<edge>> *edges; // vector of edges for the graph-structured categorical levels
  std::vector<int> *K; // number of levels per categorical variable
  double a_cat;
  double b_cat;
  int* graph_split; // do we split categorical variables using the supplied graphs?
  int graph_cut_type; // determines how we generate the partition
  
  bool rc_split;
  double* prob_rc; // prob. of proposing a random combination split. almost always set to 0
  double* theta_rc;
  int* rc_var_count; // the total number of variables used in ALL random combination rules
  int* rc_rule_count; // how many times do we use a random combination rule

  // hyperparameters will go here eventually
  double tau;
  double mu0; // prior mean
  
  // constructor
  tree_prior_info(){
    alpha = 0.95;
    beta = 2.0;
    prob_bd = 1.0;
    prob_b = 0.5;
    theta = 0; // 0 pointer
    var_count = 0; // 0 pointer
    rule_count = 0; // 0 pointer
    
    unif_cuts = 0.0; // 0 pointer
    cutpoints = 0.0; // 0 pointer
    
    cat_levels = 0; // 0 pointer
    edges = 0; // 0 pointer
    K = 0; // 0 pointer
    a_cat = 0; // 0 pointer
    b_cat = 0; // 0 pointer
  
    graph_split = 0; // 0 pointer
    graph_cut_type = 0;
    
    rc_split = false;
    prob_rc = 0; // 0 pointer
    theta_rc = 0; // 0 pointer
    rc_var_count = 0; // 0 pointer
    rc_rule_count = 0; // 0 pointer
    tau = 1.0;
    mu0 = 0.0;
  }
};

// silly class to convert sets of integers into character strings
class set_str_conversion{
public:
  std::map<std::string,char> str_to_hex_lookup;
  std::map<char, std::string> hex_to_str_lookup;
  
  set_str_conversion(){
    str_to_hex_lookup["0000"] = '0';
    str_to_hex_lookup["0001"] = '1';
    str_to_hex_lookup["0010"] = '2';
    str_to_hex_lookup["0011"] = '3';
    str_to_hex_lookup["0100"] = '4';
    str_to_hex_lookup["0101"] = '5';
    str_to_hex_lookup["0110"] = '6';
    str_to_hex_lookup["0111"] = '7';
    str_to_hex_lookup["1000"] = '8';
    str_to_hex_lookup["1001"] = '9';
    str_to_hex_lookup["1010"] = 'a';
    str_to_hex_lookup["1011"] = 'b';
    str_to_hex_lookup["1100"] = 'c';
    str_to_hex_lookup["1101"] = 'd';
    str_to_hex_lookup["1110"] = 'e';
    str_to_hex_lookup["1111"] = 'f';
    
    hex_to_str_lookup['0'] = "0000";
    hex_to_str_lookup['1'] = "0001";
    hex_to_str_lookup['2'] = "0010";
    hex_to_str_lookup['3'] = "0011";
    hex_to_str_lookup['4'] = "0100";
    hex_to_str_lookup['5'] = "0101";
    hex_to_str_lookup['6'] = "0110";
    hex_to_str_lookup['7'] = "0111";
    hex_to_str_lookup['8'] = "1000";
    hex_to_str_lookup['9'] = "1001";
    hex_to_str_lookup['a'] = "1010";
    hex_to_str_lookup['b'] = "1011";
    hex_to_str_lookup['c'] = "1100";
    hex_to_str_lookup['d'] = "1101";
    hex_to_str_lookup['e'] = "1110";
    hex_to_str_lookup['f'] = "1111";
  }
  
  std::string set_to_hex(int &K, std::set<int> &vals){
    // we divide the full set {0, 1, ... , K-1} into blocks of 4
    // block 0 {0,1,2,3}, block 1 {4,5,6,7}, etc.
    // we sweep over each block and see whether or not each element is in the set vals
    // this creates a binary string of length 4, which we then convert into a single character w/ our lookup table
    
    int num_blocks = K/4;
    std::string tmp_str(4,'0'); // temporary string of length 4, overwritten with each block
    std::string hex_str(num_blocks+1,'0'); // our outputted string, initialized for the empty set
    std::map<std::string, char>::iterator str_ch_it; // iterator for looking up in str_to_hex_lookup
    
    for(int blk_id = 0; blk_id <= num_blocks; blk_id++){
      tmp_str.assign(4,'0'); // reset the temporary string to all 0's
      for(int j = 0; j < 4; j++){
        if(vals.count(4*blk_id + j) == 1){
          // if the integer 4*blk_id + j is in the set vals, we make the j-th element of tmp_str = 1
          tmp_str[j] = '1';
        }
      } // closes loop over elements of each block
      str_ch_it = str_to_hex_lookup.find(tmp_str);
      if(str_ch_it == str_to_hex_lookup.end()){
        Rcpp::Rcout << "[set_to_hex]: temporary string " << tmp_str << std::endl;
        Rcpp::stop("string not found in str_to_hex_lookup!");
      } else{
        hex_str[blk_id] = str_ch_it->second;
      }
    } // closes loop over the blocks
    return hex_str;
  }
  
  std::set<int> hex_to_set(int &K, std::string &hex_str){
    
    int num_blocks = K/4;
    if(hex_str.size() != num_blocks+1){
      Rcpp::Rcout << "[hex_to_set]: hex_str = " << hex_str << " is wrong size" << std::endl;
      Rcpp::Rcout << "[hex_to_set]: for K = " << K << " values, hex_str must be of length " << num_blocks+1 << std::endl;
      Rcpp::stop("hex_str is of wrong size!");
    }
    std::map<char, std::string>::iterator ch_str_it; // iterator for looking up in hex_to_str_lookup
    std::string tmp_str;
    std::set<int> vals;
    
    for(int blk_id = 0; blk_id <= num_blocks; blk_id++){
      // std::string's [] lets us look up on a character-by-character basis
      ch_str_it = hex_to_str_lookup.find(hex_str[blk_id]);
      if(ch_str_it == hex_to_str_lookup.end()){
        Rcpp::Rcout << "[hex_to_set]: character " << hex_str[blk_id] << std::endl;
        Rcpp::stop("character not found in hex_to_str_lookup!");
      } else{
        tmp_str = ch_str_it->second;
        for(int j = 0; j < 4; j++){
          if(tmp_str[j] == '1') vals.insert(4*blk_id+j);
        }
      } // closes if/else checking that element of hex_str is a key in hex_to_set_lookup
    } // closes loop over the elements of hex_str
    
    return vals;
  }
  
}
;

inline void parse_cutpoints(std::vector<std::set<double>> &cutpoints, int p_cont, Rcpp::List &tmp_cutpoints, Rcpp::LogicalVector &unif_cuts)
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
      //std::set<double> xi_set;
      //for(int l = 0; l < cutpoints_vec.size(); l++) xi_set.insert(cutpoints_vec[l]);
      //cutpoints.push_back(xi_set);
    }
  } else{
    Rcpp::Rcout << "p_cont = " << p_cont;
    Rcpp::Rcout << "  cutpoints_list.size() = " << tmp_cutpoints.size() << std::endl;
    Rcpp::Rcout << "  unif_cuts.size() = " << unif_cuts.size() << std::endl;
    Rcpp::stop("cutpoints_list & unif_cuts needs to have length p_cont!");
  }
}

inline void parse_cat_levels(std::vector<std::set<int>> &cat_levels, std::vector<int> &K, int &p_cat, Rcpp::List &tmp_cat_levels)
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
// we will pass a list of Rcpp::NumericMatrices, which are computed using igraph::get_data_frame
// this function reads those matrices and builds a vector of edges
inline void parse_edge_mat(std::vector<edge> &edges, Rcpp::NumericMatrix &edge_mat, int &n_vertex)
{
  int n_edges = edge_mat.rows();
  edges.clear();
  if(edge_mat.cols() == 3){
    for(int i = 0; i < n_edges; i++){
      edges.push_back(edge( (int) edge_mat(i,0), (int) edge_mat(i,1), edge_mat(i,2)));
    }
  } else if(edge_mat.cols() == 2){
    for(int i = 0; i < n_edges; i++){
      edges.push_back(edge( (int) edge_mat(i,0), (int) edge_mat(i,1), 1.0));
    }
  } else{
    Rcpp::stop("[parse_edge_mat]: The matrix edge_mat must have 2 columns (unweighted graph) or 3 columns (weighted graph)");
  }
}

// takes in the List of edge_mat's
inline void parse_graphs(std::vector<std::vector<edge>> &edges, int &p_cat, std::vector<int> &K, Rcpp::List &tmp_edge_mats, Rcpp::LogicalVector &graph_split)
{
  edges.clear();
  edges.resize(p_cat, std::vector<edge>());
  if(tmp_edge_mats.size() == p_cat){
    for(int j = 0; j < p_cat; j++){
      if(graph_split(j) == 1){
        Rcpp::NumericMatrix edge_mat = Rcpp::as<Rcpp::NumericMatrix>(tmp_edge_mats[j]);
        parse_edge_mat(edges[j], edge_mat, K[j]);
      } else{
        // do nothing
      }
    }
  } else{
    Rcpp::Rcout << "[parse_graphs]: detected " << p_cat << " categorical variables";
    Rcpp::Rcout << " edge_mat_list has length " << tmp_edge_mats.size() << std::endl;
    Rcpp::stop("edge_mat_list must have length equal to p_cat!");
  }
}



 
inline void parse_training_data(int &n_train, int &p_cont, int &p_cat, Rcpp::NumericMatrix &tX_cont_train, Rcpp::IntegerMatrix &tX_cat_train)
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

inline void parse_testing_data(int &n_test, Rcpp::NumericMatrix &tX_cont_test, Rcpp::IntegerMatrix &tX_cat_test, const int &p_cat, const int &p_cont)
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
#endif
