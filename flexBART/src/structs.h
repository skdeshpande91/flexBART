#ifndef GUARD_structs_h
#define GUARD_structs_h

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




class rule_t{
public:
  bool is_cat; // is it a categorical split
  int v_raw; // for internal use
  int v_aa; // used to access X_cat
  double c; // cutpoint
  int v_cat; // index of variable of the categorical variable on which we split (always between 0 and p_cat)
  std::set<int> l_vals; // holds unique values of levels of v associated w/ left child
  std::set<int> r_vals; // holds unique values of levels of v associated w/ right child
  
  rule_t(){
    is_cat = false;
    v_raw = 0;
    v_aa = 0;
    c = 0.0;
    v_cat = 0;
    l_vals = std::set<int>();
    r_vals = std::set<int>();
  }
  void clear(){
    is_cat = false;
    v_raw = 0;
    v_aa = 0;
    c = 0.0;
    v_cat = 0;
    l_vals.clear();
    r_vals.clear();
  }
};

// structure for diagnostics tracking how often we propose/reject certain types of rules
struct rule_diag_t{
  int grow_prop;
  int prune_prop;
  int grow_rej;
  int prune_rej;
  int aa_prop;
  int aa_rej;
  int cat_prop;
  int cat_rej;
  rule_diag_t(){grow_prop = 0; prune_prop = 0; grow_rej = 0; prune_rej = 0; aa_prop = 0; aa_rej = 0; cat_prop = 0; cat_rej = 0; }
  void reset(){grow_prop = 0; prune_prop = 0; grow_rej = 0; prune_rej = 0; aa_prop = 0; aa_rej = 0; cat_prop = 0; cat_rej = 0;}
};

// class holding data dimensions and pointers to the covariate data
class data_info{
public:
  int n; // number of observations
  int R; // number of ensembles
  int p_cont; // number of continuous predictors
  int p_cat; // number of categorical predictors
  int p; // total number of predictors (likely will never every use this)
  double* z; // for vcbart, these are the linear predictors/weights on ensembles
  double* x_cont; // pointer to the matrix of continuous predictors
  int* x_cat; // pointer to matrix of categorical predictors (levels coded as integers, beginning with 0)
  double* rp; // partial residual;
  data_info(){n = 0; R = 0; p_cont = 0; p_cat = 0; p = 0; z = 0; x_cont = 0; x_cat = 0;rp = 0;}
};


// structure for nested categorical variables:
// lo: low resolution variable index
// hi: high resolution variable index
// map: keys are values of x_lo and value are the x_hi values contained within each value of x_lo
struct hi_lo_map{
  int hi;
  int lo;
  std::map<int, std::set<int>> map;
  hi_lo_map(){hi = 0; lo = 0; map = std::map<int, std::set<int>>();}
  hi_lo_map(int h, int l){hi = h; lo = l; map = std::map<int, std::set<int>>();}
  hi_lo_map(int h, int l, std::map<int, std::set<int>> map_arg){
    hi = h;
    lo = l;
    map = map_arg;
  }

  bool find(int h, int l){
    if( (hi == h && lo == l) || (hi == l && lo == h)) return true;
    else return false;
  }
};


// holds hyperparameters for regression tree prior
class tree_prior_info{
public:
  
  double prob_bd; // prob of proposing a grow (birth) or prune (death) move. almost always set to 1
  double prob_b; // prob of proposing a grow (birth) move. almost always set to 0.5
  std::vector<std::set<double> >* cutpoints; // holds cutpoints if there are any
  std::vector<std::set<int>>* cat_levels; // holds the levels of the categorical variables
  std::vector<std::vector<edge>> *edges; // vector of edges for the graph-structured categorical levels
  std::vector<hi_lo_map>* nesting; // maps between levels of hi and lo resolution variables
  bool nest_v; // should nesting structure influence choice of splitting variable
  bool nest_c; // should nesting structure include choice of cutset
  int nest_v_option; // how should nesting structure influence choice of splitting variable
  int graph_cut_type; // determines how graph is partitioned
  // ensemble specific stuff
  double alpha; // 1st parameter of the branching process prior
  double beta; // 2nd parameter in branching process prior
  // eventually will need stuff about the covariate graph
  double tau;
  double mu0;
  
  std::vector<double>* theta; // prob. that we pick one variable out of p_cont + p_cat
  std::vector<int> *var_count; // counts how many times we split on a single variable
  int* rule_count; // how many total rules are there in the ensemble

  // unif_cuts passed an Rcpp::LogicalVector
  //int* unif_cuts; // unif_cuts = 1 to draw cuts uniformly, = 0 to draw from pre-specified cutpoints, = minimum integer for NA
  
  //std::vector<int> *K; // number of levels per categorical variable

  //int* graph_split; // do we split categorical variables using the supplied graphs?
  //int graph_cut_type; // determines how we generate the partition

  // hyperparameters will go here eventually
  //double tau;
  //double mu0; // prior mean
  
  // constructor
  tree_prior_info(){
    prob_bd = 1.0;
    prob_b = 0.5;
    cutpoints = 0;
    cat_levels = 0;
    edges = 0;
    graph_cut_type = 2;

    nesting = 0;
    nest_v = true;
    nest_v_option = 3;
    nest_c = true;
    
    alpha = 0.95;
    beta = 2.0;
    theta = 0; // 0 pointer
    var_count = 0; // 0 pointer
    rule_count = 0; // 0 pointer
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
#endif
