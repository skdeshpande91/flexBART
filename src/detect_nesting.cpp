#include "graph_funs.h"


// [[Rcpp::export(".detect_nesting")]]
Rcpp::List detect_nesting(arma::mat X_cat,Rcpp::IntegerVector K, int p_cont = 0)
{
  int p_cat = X_cat.n_cols;
  std::vector<std::vector<int>> hi_lo_vec; // records pairs of nested variables
  std::vector<std::pair<int,int>> hi_lo_pair; // vector mapping each lo-res value to its hi-res value
  
  Rcpp::List out; // output container
  
  // if we're missing any categorical level, we need to skip the search for nesting variables entirely
  
  bool missing_level = false;
  for(int j = 0; j < p_cat; ++j){
    arma::vec values = X_cat.col(j);
    for(int k = 0; k < K(j); ++k){
      arma::uvec index = arma::find(values == k);
      if(index.size() == 0){
        missing_level = true;
        break;
      }
    }
    if(missing_level) break;
  }
  
  
  if(missing_level) out["nesting"] = R_NilValue;
  else{
    // every level of every categorical variable is available.
    // loop over all pairs of variables to check nesting
    for(int j = 0; j < p_cat-1; ++j){
      for(int jj = j+1; jj < p_cat; ++jj){
        // hi-resolution variable has more values
        // lo-resolution variable has fewer values
        int hi = 0;
        int lo = 0;
        if(K(j) < K(jj)){
          hi = jj;
          lo = j;
        } else{
          hi = j;
          lo = jj;
        }
        arma::vec hi_values = X_cat.col(hi);
        arma::vec lo_values = X_cat.col(lo);
        bool nested = true;
        std::vector<int> tmp_hi_lo_vec(K(hi));
        std::set<int> used_lo_vals; // stores the lo values that are detected
        
        for(int k = 0; k < K(hi); ++k){
          arma::uvec hi_index = arma::find(hi_values == k);
          if(hi_index.size() == 0){
            // can't find the k-th value of hi
            // shouldn't be hit but we need to be safe
            nested = false;
            break;
          } else{
            arma::vec tmp_lo = arma::unique(lo_values.elem(hi_index));
            if(tmp_lo.size() > 1){
              // same hi-res value occurs with multiple lo-res values
              // so hi is not nested within lo
              nested = false;
              break;
            } else{
              // hi-res value occurs with only one lo-res value
              // we will create the vector recording the lo value to which every hi value belongs
              tmp_hi_lo_vec[k] = tmp_lo(0);
              used_lo_vals.insert(tmp_lo(0));
            } // closes if/else checking if hi-res value occurs with multiple lo-res values
          } // closes if/else checking whether hi-res value appears
        } // closes loop over hi-res variables
        if(nested){
          // we've finished loop. hi appears to be nested within lo
          // need to check whether any lo value appears empty (i.e., this lo value does not appear with a hi value in X_cat)
          if(used_lo_vals.size() != K(lo)){
            // there's a lo value that does not appear with any hi value.
            // conservatively we need to declase this is not nested
            break;
            nested = false;
          } else{
            // every lo value contains at least one hi value (which isn't shared w/ any other lo-value)
            // add the pair to the running list of nested variables
            hi_lo_pair.push_back(std::pair<int,int>(hi,lo));
            hi_lo_vec.push_back(tmp_hi_lo_vec);
          } // closes if/else checking whether all lo values contain a hi-value
        } // closes if checking that j & jj are potentially nested
      } // closes inner loop over jj
    } // closes outer loop over j
    
    if(hi_lo_pair.size() > 0){
      Rcpp::List nest_out(hi_lo_pair.size());
      for(int i = 0; i < hi_lo_pair.size();++i){
        Rcpp::List tmp_nest_info(2);
        Rcpp::IntegerVector tmp_pair(2);
        tmp_pair(0) = hi_lo_pair[i].first + p_cont; // remember to offset by p_cont
        tmp_pair(1) = hi_lo_pair[i].second + p_cont; // remember to offset by p_cont
        tmp_nest_info[0] = tmp_pair;
        Rcpp::IntegerVector tmp_nest_map(K[hi_lo_pair[i].first]);
        for(int ii = 0; ii < hi_lo_vec[i].size(); ++ii) tmp_nest_map[ii] = hi_lo_vec[i][ii];
        tmp_nest_info[1] = tmp_nest_map;
        nest_out[i] = tmp_nest_info;
      }
      out["nesting"] = nest_out;
    } else{
      out["nesting"] = R_NilValue;
    }
  } // closes if/else checking that all levels of all variables present
  
  return(out);
}
