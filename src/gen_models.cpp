#include "gen_models.h"

#include <cmath>
#include <vector>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

double GenModel::proposal_mu_single(double &m, const int &nid, suff_stat &ss, data_info &di, tree_prior_info &tree_pi)
{
  // NOTE: assuming data_info has pointers to covariate matrix Z, current fit of parameters lambda, and response vector Y_train
  // pass m & v at desired initial values
  // j in the index of the covariate we are updating
  // function updates m & v for proposal distribution
  int iter = 0;
  int max_iter = 100;
  double U = 1;
  double I = 1;

  // Fisher scoring
  while(abs(U) > pow(I, 0.5) / 10 && iter < max_iter){
    // user passes initial value of m
    U = 0;
    I = 0;
    if(ss.count(nid) == 1){
      suff_stat_it ss_it = ss.find(nid);
      if(ss_it->second.size() > 0){
        for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
          // Rcpp::Rcout << "Iteration " << *it << ": lambda = " << di.lambda[*it] << ", rp = " << di.rp[*it] << std::endl;
          double y = inv_link(di.lambda[*it]) + di.rp[*it];
          double theta = inv_link(di.lambda[*it] + m);
          U += score(y, theta);
          I += jacobian(theta);
          // output message for debugging
          // Rcpp::Rcout << "Iteration " << iter << ": y = " << y << ", theta = " << theta << ", U = " << U << ", I = " << I << std::endl;
        }
      }
    }
    U -= m / pow(tree_pi.tau, 2.0);
    I += pow(tree_pi.tau, -2.0);
    m += U/I; 
    iter += 1;
    // output message for debugging
    // Rcpp::Rcout << "Iteration " << iter << ": U = " << U << ", I = " << I << ", m = " << m << std::endl;
  }

  // give a warning if the scoring loop terminates by max iterations
  if (iter >= max_iter) Rcpp::Rcout << "[compute_m]: WARNING! Laplace approximation did not converge after " << iter << " iterations!" << std::endl;

  return m;
}

double GenModel::proposal_mu_multi(double &m, const int &nid, suff_stat &ss, int &r, data_info &di, tree_prior_info &tree_pi)
{
  // NOTE: assuming data_info has pointers to covariate matrix Z, current fit of parameters lambda, and response vector Y_train
  // pass m & v at desired initial values
  // j in the index of the covariate we are updating
  // function updates m & v for proposal distribution
  int iter = 0;
  int max_iter = 100;
  double U = 1;
  double I = 1;

  // compute m with Fisher scoring method
  while(abs(U) > pow(I, 0.5) / 10 && iter < max_iter){
    // user passes initial value of m
    U = 0;
    I = 0;
    if(ss.count(nid) == 1){
      suff_stat_it ss_it = ss.find(nid);
      if(ss_it->second.size() > 0){
        for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
          double y = inv_link(di.lambda[*it]) + di.rp[*it];
          double z = di.z[r + (*it) * di.R];
          double theta = inv_link(di.lambda[*it] + z * m);
          U += z * score(y, theta);
          I += pow(z, 2.0) * jacobian(theta);
        }
      }
    }
    U -= m / pow(tree_pi.tau, 2.0);
    I += pow(tree_pi.tau, -2.0);
    m += U/I; 
    iter += 1;
    // output message for debugging
    // Rcpp::Rcout << "Iteration " << it << ": u = " << u << ", I = " << I << std::endl;
  }

  // give a warning if the scoring loop terminates by max iterations
  if (iter >= max_iter) Rcpp::Rcout << "[compute_m]: WARNING! Laplace approximation did not converge after " << iter << " iterations!" << std::endl;

  return m;
}

double GenModel::proposal_var_single(double m, const int &nid, suff_stat &ss, data_info &di, tree_prior_info &tree_pi)
{
  // compute v from m
  double I = pow(tree_pi.tau, -2);
  if(ss.count(nid) == 1){
    suff_stat_it ss_it = ss.find(nid);
    if(ss_it->second.size() > 0){
      for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
        double theta = inv_link(di.lambda[*it] + m);
        I += jacobian(theta);
      }
    }
  }
  return pow(I, -0.5);
}

double GenModel::proposal_var_multi(double m, const int &nid, suff_stat &ss, int &r, data_info &di, tree_prior_info &tree_pi)
{
  // compute v from m
  double I = pow(tree_pi.tau, -2);
  if(ss.count(nid) == 1){
    suff_stat_it ss_it = ss.find(nid);
    if(ss_it->second.size() > 0){
      for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
        double z = di.z[r + (*it) * di.R];
        double theta = inv_link(di.lambda[*it] + z * m);
        I += pow(z, 2.0) * jacobian(theta);
      }
    }
  }
  return pow(I, -0.5);
}

double GenModel::compute_node_lik_single(double &mu, const int &nid, suff_stat &ss, data_info &di, tree_prior_info &tree_pi)
{
  double theta = 0.0;
  if(ss.count(nid) == 1){
    suff_stat_it ss_it = ss.find(nid);
    if(ss_it->second.size() > 0){
      for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
        double y = inv_link(di.lambda[*it]) + di.rp[*it];
        theta += log_lik(y, di.lambda[*it] + mu);
      }
    }
  }
  return theta - 0.5 * log(2 * M_PI) - log(tree_pi.tau) - 0.5 * pow(((mu - tree_pi.mu0) / tree_pi.tau), 2.0);
}

double GenModel::compute_node_lik_multi(double &mu, const int &nid, suff_stat &ss, int &r, data_info &di, tree_prior_info &tree_pi)
{
  double theta = 0.0;
  if(ss.count(nid) == 1){
    suff_stat_it ss_it = ss.find(nid);
    if(ss_it->second.size() > 0){
      for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
        double y = inv_link(di.lambda[*it]) + di.rp[*it];
        double z = di.z[r + (*it) * di.R];
        theta += log_lik(y, di.lambda[*it] + z * mu);
      }
    }
  }
  return theta - 0.5 * log(2 * M_PI) - log(tree_pi.tau) - 0.5 * pow(((mu - tree_pi.mu0) / tree_pi.tau), 2.0);
}

double Logit::link(double theta) {
    // logit function
    return log(theta) - log(1 - theta);
}

double Logit::inv_link(double lambda) {
    // sigmoid function
    return 1. / (1. + exp(-lambda));
}

double Logit::log_lik(double y, double lambda) {
    return y * log(inv_link(lambda)) + (1 - y) * log(1 - inv_link(lambda));
}

double Logit::score(double y, double theta) {
    return y - theta;
}

double Logit::jacobian(double theta) {
    return theta * (1 - theta);
}

double Poisson::link(double theta) {
    return log(theta);
}

double Poisson::inv_link(double lambda) {
    return exp(lambda);
}

double Poisson::log_lik(double y, double lambda) {
    // use tgamma(y + 1) instead of factorial(y) to avoid overflow
    return y * log(inv_link(lambda)) - inv_link(lambda) - log(tgamma(y + 1));
}

double Poisson::score(double y, double theta) {
    return y - theta;
}

double Poisson::jacobian(double theta) {
    return theta;
}