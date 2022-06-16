#include "funs.h"

// [[Rcpp::export]]
Rcpp::NumericMatrix summarize_post_pred(Rcpp::NumericMatrix fit_samples,
                                        Rcpp::NumericVector sigma_samples)
{
  Rcpp::RNGScope scope;
  RNG gen;
  
  int n = fit_samples.ncol();
  int nd = fit_samples.nrow();
  //double l95 = 0.0;
  //double u95 = 0.0;
  Rcpp::NumericMatrix results(n,3);
  arma::vec probs = {0.025, 0.975};
  arma::vec tmp_quantiles(2);
  
  arma::vec tmp_samples = arma::zeros<arma::vec>(nd);
  for(int i = 0; i < n; i++){
    tmp_samples.zeros();
    for(int iter = 0; iter < nd; iter++){
      results(i,0) += fit_samples(iter,i)/( (double) nd); // for the posterior mean
      tmp_samples(iter) = fit_samples(iter,i) + sigma_samples(iter) * gen.normal();
    }
    tmp_quantiles = arma::quantile(tmp_samples, probs);
    results(i,1) = tmp_quantiles(0);
    results(i,2) = tmp_quantiles(1);
  }
  return results;
}
