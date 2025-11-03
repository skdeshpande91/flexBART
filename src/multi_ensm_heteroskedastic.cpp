#include "update_tree_heteroskedastic.h"
#include "update_tree_gen.h"
#include "gen_models.h"
#include "data_parsing_funs.h"
#include "funs.h"
// [[Rcpp::export("._multi_fit_heteroskedastic")]]
Rcpp::List multi_fit_heteroskedastic(Rcpp::NumericVector Y_train,
                                    Rcpp::IntegerMatrix cov_ensm,
                                    Rcpp::IntegerMatrix cov_var,
                                    Rcpp::NumericMatrix tZ_train,
                                    Rcpp::NumericMatrix tX_cont_train,
                                    Rcpp::IntegerMatrix tX_cat_train,
                                    Rcpp::NumericMatrix tZ_test,
                                    Rcpp::NumericMatrix tX_cont_test,
                                    Rcpp::IntegerMatrix tX_cat_test,
                                    Rcpp::Nullable<Rcpp::List> cutpoints_list,
                                    Rcpp::Nullable<Rcpp::List> cat_levels_list,
                                    Rcpp::Nullable<Rcpp::List> edge_mat_list,
                                    Rcpp::Nullable<Rcpp::List> nest_list,
                                    int graph_cut_type,
                                    bool sparse, double a_u, double b_u,
                                    bool nest_v, int nest_v_option, bool nest_c,
                                    Rcpp::IntegerVector M_vec,
                                    Rcpp::NumericVector alpha_vec, Rcpp::NumericVector beta_vec,
                                    Rcpp::NumericVector mu0_vec, Rcpp::NumericVector tau_vec,
                                    int nd, int burn, int thin,
                                    int max_iter,
                                    bool save_samples,
                                    bool save_trees,
                                    bool verbose,
                                    int print_every)
{
  Rcpp::RNGScope scope;
  RNG gen;

  GenModel* gmp = new Sigma(); // generalized model pointer
  
  set_str_conversion set_str; // for converting sets of integers into strings
  
  // BEGIN: get dimensions of training data
  int n_train = tZ_train.cols(); // how many training observations
  int R = tZ_train.rows() + 1; // number of ensembles + sigma ensemble
  int p_cont = 0;
  int p_cat = 0;
  if(tX_cont_train.size() > 0) p_cont = tX_cont_train.rows();
  if(tX_cat_train.size() > 0) p_cat = tX_cat_train.rows();
  int p = p_cont + p_cat;
  int n_test = 0;
  if(tZ_test.size() > 0) n_test = tZ_test.cols(); // how many test set observations
  // END: get dimensions of testing data
  
  // BEGIN: set cutpoints & categorical levels + parse network structure
  std::vector<std::set<double>> cutpoints;
  if(p_cont > 0) parse_cutpoints(cutpoints, p_cont, cutpoints_list);
  
  std::vector<std::set<int>> cat_levels;
  std::vector<std::vector<edge>> edges;
  if(p_cat > 0){
    parse_cat_levels(cat_levels, p_cat, cat_levels_list);
    parse_graphs(edges, p_cat, edge_mat_list);
  }
  // END: set cutpoints & categorical levels + parse network structure
  
  // BEGIN: build graph encoding nesting relationships b/w categorical predictors
  std::vector<hi_lo_map> nesting;
  std::vector<edge_map> nest_graph_in;
  std::vector<edge_map> nest_graph_out;
  std::vector<std::map<int, std::set<int>>> nest_graph_components;
  parse_nesting(nesting, nest_graph_in, nest_graph_out, nest_graph_components, p_cont, cov_ensm, cat_levels, nest_list);
  // END: build graph encoding nesting relationships b/w categorical predictors

  // BEGIN: create splitting probabilities
  std::vector<std::vector<int>> var_count(R, std::vector<int>(p, 0));
  std::vector<int> rule_count(R, 0);
  std::vector<std::vector<double>> theta(R, std::vector<double>(p, 0.0));
  std::vector<double> u(R);
  
  // BART ensemble splitting probabilities
  for(int r = 0; r < R; ++r){
    int n_avail_vars = 0;
    for(int j = 0; j < p; ++j){
      if (r < (R-1)){
        if(cov_ensm(j,r) == 1){
          theta[r][j] = 1.0;
          ++n_avail_vars;
        }
      } else{ // sigma ensemble
        if(cov_var(j,0) == 1){
          theta[r][j] = 1.0;
          ++n_avail_vars;
        }
      }
    } // closes loop over variables
    if(n_avail_vars == 0){
      Rcpp::Rcout << "Ensemble r = " << r << " no covariates detected!" << std::endl;
      Rcpp::stop("At least one covariate needed for each ensemble!");
    } else{
      for(int j = 0; j < p; ++j) theta[r][j] /= (double) n_avail_vars;
      u[r] = 1.0/(1.0 + (double) n_avail_vars);
    }
  } // closes loop over ensembles
  // END: create splitting probabilities
  
  // BEGIN: initialize containers for residuals and fit
  double* residual = new double[n_train]; // holds the residual of the mean ensemble
  double* lambda = new double[n_train]; // holds the current estimate of log(sigma(x)^2)
  int tmp_n_test = 1;
  if(n_test > 0) tmp_n_test = n_test;
  double* tmp_fit_test = new double[tmp_n_test]; // for holding test set fits temporarily
  double* tmp_sigma_test = new double[tmp_n_test]; // for holding test set sigmas temporarily
  // END: initialize containers for residuals and fit
  
  // BEGIN: create data_info objects for training & testing
  data_info di_train;
  di_train.n = n_train;
  di_train.p_cont = p_cont;
  di_train.p_cat = p_cat;
  di_train.p = p;
  di_train.R = R - 1; // include all ensembles for z matrix indexing
  di_train.z = tZ_train.begin();
  if(p_cont > 0) di_train.x_cont = tX_cont_train.begin();
  if(p_cat > 0) di_train.x_cat = tX_cat_train.begin();
  di_train.rp = residual;
  di_train.lambda = lambda;
  
  // set up the data info object for testing data
  data_info di_test;
  if(n_test > 0){
    di_test.n = n_test;
    di_test.p_cont = p_cont;
    di_test.p_cat = p_cat;
    di_test.p = p;
    di_test.R = R - 1; // include all ensembles for z matrix indexing
    di_test.z = tZ_test.begin();
    if(p_cont > 0) di_test.x_cont = tX_cont_test.begin();
    if(p_cat > 0)  di_test.x_cat = tX_cat_test.begin();
  }
  // END: create data_info objects for training & testing
 
  // BEGIN: creating tree prior info object
  std::vector<tree_prior_info> tree_pi_vec(R);
  for(int r = 0; r < R; ++r){
    // start by assigning stuff relevant to all ensembles
    if(p_cont > 0) tree_pi_vec[r].cutpoints = &cutpoints;
    if(p_cat > 0){
      tree_pi_vec[r].cat_levels = &cat_levels;
      tree_pi_vec[r].edges = &edges;
      tree_pi_vec[r].graph_cut_type = graph_cut_type;
      tree_pi_vec[r].nesting = &nesting;
      tree_pi_vec[r].nest_in = &(nest_graph_in[r]);
      tree_pi_vec[r].nest_out = &(nest_graph_out[r]);
      tree_pi_vec[r].nest_components = &(nest_graph_components[r]);
    }
    tree_pi_vec[r].nest_v = nest_v;
    tree_pi_vec[r].nest_v_option = nest_v_option;
    tree_pi_vec[r].nest_c = nest_c;
    tree_pi_vec[r].theta = &(theta[r]); // only use theta if nest_v_option = false
    tree_pi_vec[r].var_count = &(var_count[r]);
    tree_pi_vec[r].rule_count = &(rule_count[r]);
    tree_pi_vec[r].alpha = alpha_vec[r];
    tree_pi_vec[r].beta = beta_vec[r];
    tree_pi_vec[r].mu0 = mu0_vec[r];
    tree_pi_vec[r].tau = tau_vec[r];
    tree_pi_vec[r].max_iter = max_iter;
  }
  // END: creating tree prior info object
  
  // BEGIN: initialize stuff for main MCMC loop
  int total_draws = 1 + burn + (nd-1)*thin;
  int sample_index = 0;
  int accept = 0;
  int total_accept = 0; // counts how many trees we change in each iteration
  double tmp_mu; // for holding the value of mu when we're doing the backfitting
  // END: initialize stuff for main MCMC loop
  
  // BEGIN: initialize vector of tree ensembles & maps of observations to leafs
  std::vector<std::vector<tree>> t_vec;
  std::vector<std::vector<suff_stat>> ss_train_vec;
  std::vector<std::vector<suff_stat>> ss_test_vec;
  
  for(int r = 0; r < R; ++r){
    t_vec.push_back(std::vector<tree>(M_vec[r]));
    ss_train_vec.push_back(std::vector<suff_stat>(M_vec[r]));
    ss_test_vec.push_back(std::vector<suff_stat>(M_vec[r]));
  }
  // END: initialize vector of tree ensembles & maps of observations to leafs
  
  for(int i = 0; i < n_train; ++i) residual[i] = Y_train[i]; // start with all trees initialized at zero
  for(int r = 0; r < R-1; ++r){
    for(int m = 0; m < M_vec[r]; ++m){
      tree_traversal(ss_train_vec[r][m], t_vec[r][m], di_train);
      // get the fit of each tree
      for(suff_stat_it l_it = ss_train_vec[r][m].begin(); l_it != ss_train_vec[r][m].end(); ++l_it){
        tmp_mu = t_vec[r][m].get_ptr(l_it->first)->get_mu(); // get the value of mu in the leaf
        if(l_it->second.size() > 0){
          for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
            residual[*it] -= di_train.z[r + (*it) * di_train.R] * tmp_mu; // remove initial fit of each tree from residual
          }
        }
      }
      if(n_test > 0) tree_traversal(ss_test_vec[r][m], t_vec[r][m], di_test);
    }
  }
  
  for(int i = 0; i < n_train; ++i) lambda[i] = 0.0; // start with all trees initialized at zero (log(sigma(x)^2) = 0 implies sigma(x)^2 = 1)
  for(int m = 0; m < M_vec[R-1]; ++m){
    tree_traversal(ss_train_vec[R-1][m], t_vec[R-1][m], di_train);
    // get the fit of each tree
    for(suff_stat_it l_it = ss_train_vec[R-1][m].begin(); l_it != ss_train_vec[R-1][m].end(); ++l_it){
      tmp_mu = t_vec[R-1][m].get_ptr(l_it->first)->get_mu(); // get the value of mu in the leaf
      if(l_it->second.size() > 0){
        for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it) lambda[*it] += tmp_mu;
      }
    }
    if(n_test > 0) tree_traversal(ss_test_vec[R-1][m], t_vec[R-1][m], di_test);
  }
  // END: initialize sigma tree vector & sufficient statistics maps

  // BEGIN: initialize laplace approximation map for sigma ensemble
  std::vector<std::map<int, laplace_approx>> lap_map_vec(M_vec[R-1]);
  for(int m = 0; m < M_vec[R-1]; ++m) compute_laplace_approx_single(lap_map_vec[m], ss_train_vec[R-1][m], di_train, tree_pi_vec[R-1], *gmp);
  // END: initialize laplace approximation map for sigma ensemble
  
  // BEGIN: create output containers
  arma::vec fit_train_mean = arma::zeros<arma::vec>(n_train); // posterior mean for training data
  arma::vec fit_test_mean = arma::zeros<arma::vec>(1); // posterior mean for testing data (if any)
  arma::mat beta_train_mean = arma::zeros<arma::mat>(n_train, R-1); // posterior mean fit for each ensemble for training data
  arma::mat beta_test_mean = arma::zeros<arma::mat>(1,1); // posterior mean for each ensemble for testing data (if any)
  if(n_test > 0){
    fit_test_mean.zeros(n_test);
    beta_test_mean.zeros(n_test, R-1);
  }

  arma::vec sigma_train_mean = arma::zeros<arma::vec>(n_train); // posterior mean for training data
  arma::vec sigma_test_mean = arma::zeros<arma::vec>(1); // posterior mean for testing data (if any)
  if(n_test > 0) sigma_test_mean.zeros(n_test); // arma::set.size can initialize with garbage values
  
  arma::mat fit_train = arma::zeros<arma::mat>(1,1); // posterior samples for training data
  arma::mat fit_test = arma::zeros<arma::mat>(1,1); // posterior samples for testing data (if any)
  arma::cube beta_train = arma::zeros<arma::cube>(1,1,1);
  arma::cube beta_test = arma::zeros<arma::cube>(1,1,1);
  arma::mat sigma_train = arma::zeros<arma::mat>(1,1); 
  arma::mat sigma_test = arma::zeros<arma::mat>(1,1);

  if(save_samples){
    // if we are saving all samples, then we resize the containers accordingly
    fit_train.zeros(nd, n_train);
    beta_train.zeros(nd, n_train, R-1);
    sigma_train.zeros(total_draws, n_train); // total draws so we can assess sampler convergence
    if(n_test > 0){
      fit_test.zeros(nd, n_test);
      sigma_test.zeros(nd, n_test);
      beta_test.zeros(nd, n_test, R-1);
    }
  }
  
  arma::cube var_count_samples(nd, p, R); // how many times was a variable used in each iteration
  arma::mat total_accept_samples(total_draws, R);
  Rcpp::List tree_draws(nd);
  // END: create output containers
  
  //BEGIN: burn-in
  for(int iter = 0; iter < burn; ++iter){
    if(iter % print_every == 0){
      Rcpp::checkUserInterrupt();
      if(verbose){
        if(iter  == 0) Rcpp::Rcout << "  MCMC Iteration: " << iter+1 << " of " << total_draws << "; Warmup" << std::endl;
        else Rcpp::Rcout << "  MCMC Iteration: " << iter << " of " << total_draws << "; Warmup" << std::endl;
      }
    }
    
    for(int r = 0; r < R-1; ++r){
      total_accept = 0;
      for(int m = 0; m < M_vec[r]; ++m){
        //BEGIN: remove fit of m-th tree
        for(suff_stat_it l_it = ss_train_vec[r][m].begin(); l_it != ss_train_vec[r][m].end(); ++l_it){
          // loop over the bottom nodes in m-th tree
          tmp_mu = t_vec[r][m].get_ptr(l_it->first)->get_mu(); // get the value of mu in the leaf
          if(l_it->second.size() > 0){
            for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
              residual[*it] += di_train.z[r + (*it) * di_train.R] * tmp_mu;
            }
          } // closes if checking that leaf is non-empty and computes partial residual
        } // closes loop over leafs
        // END: remove fit of m-th tree
        
        update_tree_heteroskedastic_multi(t_vec[r][m], ss_train_vec[r][m], ss_test_vec[r][m], accept, r, di_train, di_test, tree_pi_vec[r], gen); // update the tree
        total_accept += accept;
      
        // BEGIN: restore fit of m-th tree
        for(suff_stat_it l_it = ss_train_vec[r][m].begin(); l_it != ss_train_vec[r][m].end(); ++l_it){
          tmp_mu = t_vec[r][m].get_ptr(l_it->first)->get_mu();
          if(l_it->second.size() > 0){
            for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
              residual[*it] -= di_train.z[r + (*it) * di_train.R] * tmp_mu;
            }
          } // closes if checking that leaf is non-empty and computing full residual
        } // closes loop over leafs
        // END: restore fit of m-th tree
        if(sparse) update_theta_u(theta[r], u[r], var_count[r], p, a_u, b_u, gen);
        total_accept_samples(iter, r) = total_accept; // how many trees changed in this iteration
      } // closes loop over all of the trees
    } // closes loop over all of the ensembles
      
    // BEGIN: update sigma
    total_accept = 0;
    for(int m = 0; m < M_vec[R-1]; ++m){
      //BEGIN: remove fit of m-th tree
      for(suff_stat_it l_it = ss_train_vec[R-1][m].begin(); l_it != ss_train_vec[R-1][m].end(); ++l_it){
        // loop over the bottom nodes in m-th tree
        tmp_mu = t_vec[R-1][m].get_ptr(l_it->first)->get_mu(); // get the value of mu in the leaf
        if(l_it->second.size() > 0){
          for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it) lambda[*it] -= tmp_mu;
        } // closes if checking that leaf is non-empty and computes partial residual
      } // closes loop over leafs
      // END: remove fit of m-th tree
      
      update_tree_gen_single(t_vec[R-1][m], ss_train_vec[R-1][m], ss_test_vec[R-1][m], lap_map_vec[m], accept, di_train, di_test, tree_pi_vec[R-1], *gmp, gen); // update the tree
      total_accept += accept;

      // BEGIN: restore fit of m-th tree
      for(suff_stat_it l_it = ss_train_vec[R-1][m].begin(); l_it != ss_train_vec[R-1][m].end(); ++l_it){
        tmp_mu = t_vec[R-1][m].get_ptr(l_it->first)->get_mu();
        if(l_it->second.size() > 0){
          for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it) lambda[*it] += tmp_mu;
        } // closes if checking that leaf is non-empty and computing full residual
      } // closes loop over leafs
      // END: restore fit of m-th tree
      if(sparse) update_theta_u(theta[R-1], u[R-1], var_count[R-1], p, a_u, b_u, gen);
      total_accept_samples(iter, R-1) = total_accept; // how many trees changed in this iteration
    } // closes loop over all of the trees
    // save sigma samples
    for(int i = 0; i < n_train; ++i) sigma_train(iter, i) = sqrt(gmp->inv_link(lambda[i]));
    // END: update sigma
  } // closes burn-in
  // END: burn-in
  
  //BEGIN: post-burn-in
  for(int iter = burn; iter < total_draws; ++iter){
    if(iter==total_draws-1){
      if(verbose) Rcpp::Rcout << "  MCMC Iteration: " << iter+1 << " of " << total_draws << "; Sampling" << std::endl;
    } else if(iter%print_every == 0 || (iter==burn)){
      Rcpp::checkUserInterrupt();
      if(verbose) Rcpp::Rcout << "  MCMC Iteration: " << iter << " of " << total_draws << "; Sampling" << std::endl;
    }

    for(int r = 0; r < R-1; ++r){
      total_accept = 0;
      for(int m = 0; m < M_vec[r]; ++m){
        //BEGIN: remove fit of m-th tree
        for(suff_stat_it l_it = ss_train_vec[r][m].begin(); l_it != ss_train_vec[r][m].end(); ++l_it){
          // loop over the bottom nodes in m-th tree
          tmp_mu = t_vec[r][m].get_ptr(l_it->first)->get_mu(); // get the value of mu in the leaf
          if(l_it->second.size() > 0){
            for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
              residual[*it] += di_train.z[r + (*it) * di_train.R] * tmp_mu;
            }
          } // closes if checking that leaf is non-empty and computes partial residual
        } // closes loop over leafs
        // END: remove fit of m-th tree
        
        update_tree_heteroskedastic_multi(t_vec[r][m], ss_train_vec[r][m], ss_test_vec[r][m], accept, r, di_train, di_test, tree_pi_vec[r], gen); // update the tree
        total_accept += accept;
      
        // BEGIN: restore fit of m-th tree
        for(suff_stat_it l_it = ss_train_vec[r][m].begin(); l_it != ss_train_vec[r][m].end(); ++l_it){
          tmp_mu = t_vec[r][m].get_ptr(l_it->first)->get_mu();
          if(l_it->second.size() > 0){
            for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
              residual[*it] -= di_train.z[r + (*it) * di_train.R] * tmp_mu;
            }
          } // closes if checking that leaf is non-empty and computing full residual
        } // closes loop over leafs
        // END: restore fit of m-th tree
        if(sparse) update_theta_u(theta[r], u[r], var_count[r], p, a_u, b_u, gen);
        total_accept_samples(iter, r) = total_accept; // how many trees changed in this iteration
      } // closes loop over all of the trees
    } // closes loop over all of the ensembles
    
    // BEGIN: update sigma
    total_accept = 0;
    for(int m = 0; m < M_vec[R-1]; ++m){
      //BEGIN: remove fit of m-th tree
      for(suff_stat_it l_it = ss_train_vec[R-1][m].begin(); l_it != ss_train_vec[R-1][m].end(); ++l_it){
        // loop over the bottom nodes in m-th tree
        tmp_mu = t_vec[R-1][m].get_ptr(l_it->first)->get_mu(); // get the value of mu in the leaf
        if(l_it->second.size() > 0){
          for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it) lambda[*it] -= tmp_mu;
        } // closes if checking that leaf is non-empty and computes partial residual
      } // closes loop over leafs
      // END: remove fit of m-th tree
      
      update_tree_gen_single(t_vec[R-1][m], ss_train_vec[R-1][m], ss_test_vec[R-1][m], lap_map_vec[m], accept, di_train, di_test, tree_pi_vec[R-1], *gmp, gen); // update the tree
      total_accept += accept;
    
      // BEGIN: restore fit of m-th tree
      for(suff_stat_it l_it = ss_train_vec[R-1][m].begin(); l_it != ss_train_vec[R-1][m].end(); ++l_it){
        tmp_mu = t_vec[R-1][m].get_ptr(l_it->first)->get_mu();
        if(l_it->second.size() > 0){
          for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it) lambda[*it] += tmp_mu;
        } // closes if checking that leaf is non-empty and computing full residual
      } // closes loop over leafs
      // END: restore fit of m-th tree
      if(sparse) update_theta_u(theta[R-1], u[R-1], var_count[R-1], p, a_u, b_u, gen);
      total_accept_samples(iter, R-1) = total_accept; // how many trees changed in this iteration
    } // closes loop over all of the trees

    // save sigma samples
    if(save_samples) for(int i = 0; i < n_train; ++i) sigma_train(iter, i) = sqrt(gmp->inv_link(lambda[i]));
    // END: update sigma

    if( (iter - burn)%thin == 0 ){
      sample_index = (int) ( (iter-burn)/thin);
      for(int r = 0; r < R; ++r){
        for(int j = 0; j < p; ++j) var_count_samples(sample_index,j,r) = var_count[r][j];
      }

      // BEGIN: write tree to string
      if(save_trees){
        Rcpp::List tmp_tree_draws(R);
        for(int r = 0; r < R; ++r){
          Rcpp::CharacterVector tree_string_vec(M_vec[r]);
          for(int m = 0; m < M_vec[r]; ++m) tree_string_vec[m] = write_tree(t_vec[r][m], tree_pi_vec[r], set_str);
          tmp_tree_draws[r] = tree_string_vec;
        }
        tree_draws[sample_index] = tmp_tree_draws;
      }
      // END: write tree to string

      // BEGIN: save samples
      if(save_samples){
        for(int i = 0; i < n_train; ++i){
          fit_train(sample_index,i) = Y_train[i] - residual[i];
          fit_train_mean(i) += Y_train[i] - residual[i];
          sigma_train_mean(i) += gmp->inv_link(lambda[i]);
        }
        // now compute beta_train
        for(int r = 0; r < R-1; ++r){
          for(int m = 0; m < M_vec[r]; ++m){
            for(suff_stat_it l_it = ss_train_vec[r][m].begin(); l_it != ss_train_vec[r][m].end(); ++l_it){
              tmp_mu = t_vec[r][m].get_ptr(l_it->first)->get_mu();
              if(l_it->second.size() > 0){
                for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
                  beta_train(sample_index, *it, r) += tmp_mu;
                  beta_train_mean(*it,r) += tmp_mu;
                }
              }
            }
          }
        }
      } else{
        for(int i = 0; i < n_train; ++i){
          fit_train_mean(i) += Y_train[i] - residual[i];
          sigma_train_mean(i) += gmp->inv_link(lambda[i]);
        }
      }

      if(n_test > 0){
        for(int r = 0; r < R-1; ++r){
          for(int m = 0; m < M_vec[r]; ++m){
            for(suff_stat_it l_it = ss_test_vec[r][m].begin(); l_it != ss_test_vec[r][m].end(); ++l_it){
              tmp_mu = t_vec[r][m].get_ptr(l_it->first)->get_mu();
              if(l_it->second.size() > 0){
                for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
                  int i = *it;
                  fit_test_mean(i) += di_test.z[r + i*di_test.R] * tmp_mu;
                  beta_test_mean(i,r) += tmp_mu;
                  if(save_samples){
                    fit_test(sample_index, i) += di_test.z[r + i*di_test.R] * tmp_mu;
                    beta_test(sample_index,i,r) += tmp_mu;
                  }
                } // closes loop over observations in leaf
              } // closes if checking that leaf is non-empty
            } // closes over leaves in the tree
          } // closes loop over trees
        } // closes loop over ensembles

        for(int i = 0; i < n_test; ++i) tmp_sigma_test[i] = 0.0;
        for(int m = 0; m < M_vec[R-1]; ++m){
          for(suff_stat_it l_it = ss_test_vec[R-1][m].begin(); l_it != ss_test_vec[R-1][m].end(); ++l_it){
            tmp_mu = t_vec[R-1][m].get_ptr(l_it->first)->get_mu();
            if(l_it->second.size() > 0){
              for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
                tmp_sigma_test[*it] += tmp_mu;
              } // closes loop over observations in leaf
            } // closes if checking that leaf is non-empty
          } // closes over leaves in the tree
        } // closes loop over trees

        if(save_samples){
          for(int i = 0; i < n_test; ++i){
            sigma_test(sample_index, i) = gmp->inv_link(tmp_sigma_test[i]);
            sigma_test_mean(i) += gmp->inv_link(tmp_sigma_test[i]);
          }
        } else{
          for(int i = 0; i < n_test; ++i){
            sigma_test_mean(i) += gmp->inv_link(tmp_sigma_test[i]);
          }
        } // closes if/else checking whether we're saving samples or just posterior mean
      } // closes if checking that there are test set observations
    } // closes if that checks whether we should save anything in this iteration
  } // closes post-burn-in loop
  // END: post-burn-in

  if(tree_pi_vec[R-1].convergance_warning) Rcpp::Rcout << "WARNING! At least one Laplace approximation did not converge. Consider increasing 'max_iter'." << std::endl;
  
  fit_train_mean /= ( (double) nd);
  beta_train_mean /= ( (double) nd);
  sigma_train_mean = sqrt(sigma_train_mean) / ( (double) nd);
  if(n_test > 0){
    fit_test_mean /= ( (double) nd);
    beta_test_mean /= ( (double) nd);
    sigma_test_mean = sqrt(sigma_test_mean) / ( (double) nd);
  }
  
  Rcpp::List results;
  
  results["fit_train_mean"] = fit_train_mean;
  results["beta_train_mean"] = beta_train_mean;
  results["sigma_train_mean"] = sigma_train_mean;
  if(save_samples){
    results["fit_train"] = fit_train;
    results["beta_train"] = beta_train;
    results["sigma_train"] = sigma_train;
  }
  
  if(n_test > 0){
    results["fit_test_mean"] = fit_test_mean;
    results["beta_test_mean"] = beta_test_mean;
    results["sigma_test_mean"] = sigma_test_mean;
    if(save_samples){
      results["fit_test"] = fit_test;
      results["beta_test"] = beta_test;
      results["sigma_test"] = sigma_test;
    }
  }
  results["total_accept"] = total_accept_samples;
  results["var_count"] = var_count_samples;
  if(save_trees) results["trees"] = tree_draws;

  return results;
  
}
