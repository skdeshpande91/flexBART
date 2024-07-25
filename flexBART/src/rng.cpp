#include "rng.h"

double RNG::uniform(double x, double y ){
  return R::runif(x,y);
}

double RNG::exponential(double lambda){
  return R::rexp(lambda);
}

double RNG::log_uniform(){
  return -1.0 * exponential(1.0);
  // if U ~ Uniform(0,1), then (-1 * log(U)) ~ Exponential(1)
}
// if U ~ uniform(0,1), then -log(u) ~ Exponential(1)
// Furthermore, -log(-log(U)) ~ Gumbel(0,1)
double RNG::gumbel(){
  return -1.0 * log(exponential(1.0));
}

double RNG::normal(double mu, double sd ){
  return R::rnorm(mu,sd);
}

double RNG::gamma(double shape, double scale){
  return R::rgamma(shape, 1)*scale;
}

double RNG::beta(double a1, double a2){
  const double x1 = gamma(a1, 1); return (x1 / (x1 + gamma(a2, 1)));
}

double RNG::chi_square(double df){
  return R::rchisq(df);
}

int RNG::categorical(std::vector<double> &probs){
  int n = probs.size();
  int output = 0;
  double tmp = 0.0;
  if(n > 1){
    double max_quantity = log(probs[0]) + gumbel();
    for(int i = 1; i < n; i++){
      tmp = log(probs[i]) + gumbel();
      if(tmp > max_quantity){
        max_quantity = tmp;
        output = i;
      }
    }
  } else{
    // this is totally redundant
    output = 0;
  }
  return output;
}

int RNG::categorical(std::vector<double>* probs){
  int n = probs->size();
  int output = 0;
  double tmp = 0.0;
  if(n > 1){
    double max_quantity = log(probs->at(0)) + gumbel();
    for(int i = 1; i < n; i++){
      tmp = log(probs->at(i)) + gumbel();
      if(tmp > max_quantity){
        max_quantity = tmp;
        output = i;
      }
    }
  } else{
    output = 0;
  }
  return output;
}

void RNG::dirichlet(std::vector<double> &theta, std::vector<double> &concentration)
{
  if(theta.size() != concentration.size()) Rcpp::stop("[dirichlet]{ Concentration & theta must have same length!");
  int p = theta.size();
  std::vector<double> tmp_gamma(p);
  double tmp_sum = 0.0;
  for(int j = 0; j < p; j++){
    tmp_gamma[j] = gamma(concentration[j],1.0);
    tmp_sum += tmp_gamma[j];
  }
  for(int j = 0; j < p; j++) theta[j] = tmp_gamma[j]/tmp_sum;
}

// we can probably deprecate this
int RNG::multinomial(const int &R, const std::vector<double> &probs){
  int x = 0;
  double cumsum = 0.0;
  double unif = uniform(0,1);
  for(int r = 0; r < R; r++){
    cumsum += probs[r];
    if(unif < cumsum){
      x = r;
      break;
    }
  }
  return(x);
}

// we can likely deprecate this
int RNG::multinomial(const int &R, std::vector<double>* probs){
  int x = 0;
  double cumsum = 0.0;
  double unif = uniform(0.0,1.0);
  for(int r = 0; r < R; r++){
    cumsum += probs->at(r);
    if(unif < cumsum){
      x = r;
      break;
    }
  }
  return(x);
}


arma::vec RNG::std_norm_vec(int d){
  arma::vec results(d);
  for(int i = 0; i < d; i++) results(i) = normal(0.0,1.0);
  return(results);
}

arma::mat RNG::std_norm_mat(int nrow, int ncol){
  arma::mat results(nrow, ncol);
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++) results(i,j) = normal(0.0, 1.0);
  }
  return(results);
}


arma::vec RNG::mvnormal(arma::vec m, arma::mat P){
  int d = m.size();
  if( (P.n_rows != d) | (P.n_cols != d)){
    Rcpp::Rcout << "m.size() = " << d << " P.nrow = " << P.n_rows << " P.n_cols = " << P.n_cols << std::endl;
    Rcpp::stop("[RNG::mvnormal]: m & P have incompatible dimensions!");
  }
  
  arma::mat L = arma::chol(P, "lower"); // P = L L.t()
  arma::vec nu = arma::solve(arma::trimatl(L), m); // nu = L^-1 m.
  arma::vec mu = arma::solve(arma::trimatu(L.t()), nu); // mu = (L')^-1 L^-1 m = P^-1 m
  arma::vec z = arma::solve(arma::trimatu(L.t()), std_norm_vec(d)); // Cov(z) = (L')^-1 ((L')^-1)' = (L')^-1 L^-1 = P^-1
  
  arma::vec results = mu + z;
  return(results);
}


// N(0,1) truncated to be >lo
double RNG::lo_trunc_std_norm(double lo){
  double x;
  if(lo<0) {
    x = R::rnorm(0.0, 1.0);
    while(x<lo) x = R::rnorm(0.0, 1.0);
  } else {
    double a = 0.5*(lo + sqrt(lo*lo + 4.0));
    x = R::rexp(1.0/a) + lo;
    double u = R::runif(0.0, 1.0);
    double diff = (x-a);
    double r = exp(-0.5*diff*diff);
    while(u > r) {
      x = R::rexp(1.0/a) + lo;
      u = R::runif(0.0, 1.0);
      diff = (x-a);
      r = exp(-0.5*diff*diff);
    }
  }
  return x;
}

// N(mu,1) truncated to be > lo
double RNG::lo_trunc_norm(double mean, double lo){
  return mean + lo_trunc_std_norm(lo - mean);
}

double RNG::hi_trunc_norm(double mean, double lo){
  return -1.0 * lo_trunc_norm(-1.0 * mean, -1.0 * lo);
}
