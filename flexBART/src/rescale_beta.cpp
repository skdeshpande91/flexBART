#include<RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat rescale_beta_mean(arma::mat beta_input,
                            double y_mean,
                            double y_sd,
                            Rcpp::NumericVector x_mean,
                            Rcpp::NumericVector x_sd,
                            Rcpp::IntegerVector z_col_id)
{
  int R = z_col_id.size();
  if(beta_input.n_cols != R){
    Rcpp::stop("[rescale_beta_mean]: length(z_col_id) must equal ncol(beta_input)");
  }
  int nd = beta_input.n_rows;
  
  int num_unik = Rcpp::max(z_col_id)+1; // number of unique values
  arma::mat beta_out = arma::zeros<arma::mat>(nd,num_unik);
  
  beta_out.col(0) = beta_input.col(0) * y_sd; // put in the first intercept term
  beta_out.col(0) += y_mean;

  for(int j = 1; j < R; ++j){
    beta_out.col(z_col_id(j)) += y_sd/x_sd[j] * beta_input.col(j);
    beta_out.col(0) -= y_sd/x_sd[j] * x_mean[j] * beta_input.col(j);
  }
  
  return beta_out;
}

// [[Rcpp::export]]
arma::cube rescale_beta(arma::cube beta_input,
                        double y_mean,
                        double y_sd,
                        Rcpp::NumericVector x_mean,
                        Rcpp::NumericVector x_sd,
                        Rcpp::IntegerVector z_col_id)
{
  
  int R = z_col_id.size();
  // beta_input should be nd x n x R

  
  if(beta_input.n_slices != R){
    Rcpp::stop("[rescale_beta]: length(z_col_id) must equal ncol(beta_input)");
  }
  int n = beta_input.n_cols; // number of observations
  int nd = beta_input.n_rows; // number of MCMC draws
  int num_unik = Rcpp::max(z_col_id)+1; // number of unique values
  arma::cube beta_out = arma::zeros<arma::cube>(nd,n,num_unik);
  
  beta_out.slice(0) = beta_input.slice(0) * y_sd;
  beta_out.slice(0) += y_mean;
  
  for(int j = 1; j < R; ++j){
    beta_out.slice(z_col_id(j)) += y_sd/x_sd[j] * beta_input.slice(j);
    beta_out.slice(0) -= y_sd/x_sd[j] * x_mean[j] * beta_input.slice(j);
  }
  return beta_out;
}
