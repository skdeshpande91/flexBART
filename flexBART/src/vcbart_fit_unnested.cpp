#include "update_tree.h"
#include "data_parsing_funs.h"
#include "funs.h"
// [[Rcpp::export("._vcbart_fit_unnested")]]
Rcpp::List vcbart_fit(Rcpp::NumericVector Y_train,
                      Rcpp::NumericMatrix tZ_train,
                      Rcpp::NumericMatrix tX_cont_train,
                      Rcpp::IntegerMatrix tX_cat_train,
                      double sigest,
                      Rcpp::IntegerMatrix cov_ensm,
                      Rcpp::Nullable<Rcpp::List> cutpoints_list,
                      Rcpp::Nullable<Rcpp::List> cat_levels_list,
                      Rcpp::Nullable<Rcpp::List> edge_mat_list,
                      Rcpp::NumericMatrix tZ_test,
                      Rcpp::NumericMatrix tX_cont_test,
                      Rcpp::IntegerMatrix tX_cat_test,
                      Rcpp::IntegerVector M_vec,
                      Rcpp::NumericVector alpha_vec, Rcpp::NumericVector beta_vec,
                      Rcpp::NumericVector mu0_vec, Rcpp::NumericVector tau_vec,
                      int graph_cut_type,
                      bool sparse, double a_u, double b_u,
                      double nu, double lambda,
                      int nd, int burn, int thin,
                      bool save_samples,
                      bool save_trees,
                      bool verbose,
                      int print_every)

{
  Rcpp::RNGScope scope;
  RNG gen;
  
  set_str_conversion set_str; // for converting sets of integers into strings
  
  
  // BEGIN PREPROCESSING
  
  // BEGIN: get dimensions of training data
  int n_train = tZ_train.cols(); // how many training observations
  int R = tZ_train.rows(); // how many ensembles
  int p_cont = 0;
  int p_cat = 0;
  if(tX_cont_train.size() > 1){
    p_cont = tX_cont_train.rows(); // how many continuous covariates
  }
  if(tX_cat_train.size() > 1){
    p_cat = tX_cat_train.rows(); // how many categorical covariates
  }
  int p = p_cont + p_cat;
  // END: get dimensions of training data;
  // BEGIN: get dimensions of testing data
  int n_test = 0;
  if(tZ_test.size() > 1) n_test = tZ_test.cols(); // how many test set observations
  // END: get dimensions of testing data


  // BEGIN: set cutpoints & categorical levels + parse network structure
  std::vector<std::set<double>> cutpoints;
  if(p_cont > 0) parse_cutpoints(cutpoints, p_cont, cutpoints_list);
  
  std::vector<std::set<int>> cat_levels;
  std:vector<std::vector<edge>> edges;
  if(p_cat > 0){
    parse_cat_levels(cat_levels, p_cat, cat_levels_list);
    parse_graphs(edges, p_cat, edge_mat_list);
  }
  // END: set cutpoints & categorical levels + parse network structure

  // BEGIN: create splitting probabilities
  std::vector<std::vector<int>> var_count(R, std::vector<int>(p, 0));
  std::vector<int> rule_count(R, 0);
  std::vector<std::vector<double>> theta(R, std::vector<double>(p, 0.0));
  std::vector<double> u(R);
  
  for(int r = 0; r < R; ++r){
    int n_avail_vars = 0;
    for(int j = 0; j < p; ++j){
      if(cov_ensm(j,r) == 1){
        theta[r][j] = 1.0;
        ++n_avail_vars;
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


  // BEGIN: initializing containers for residuals and fit
  double* allfit_train = new double[n_train];
  double* beta_fit_train = new double[n_train * R];
  double* residual = new double[n_train];
  
  int tmp_n_test = 1;
  if(n_test > 0) tmp_n_test = n_test;
  double* allfit_test = new double[tmp_n_test];
  double* beta_fit_test = new double[tmp_n_test * R];
  
  //BEGIN: initialize containers
  for(int i = 0; i < n_train; ++i){
    allfit_train[i] = 0.0;
    for(int r = 0; r < R; ++r) beta_fit_train[r + i * R] = 0.0;
  }
  if(n_test > 0){
    for(int i = 0; i < n_test; ++i){
      allfit_test[i] = 0.0;
      for(int r = 0; r < R; ++r) beta_fit_test[r + i * R] = 0.0;
    }
  }
  // END: initializeing containers for residuals and fit
  
  // BEGIN: creating data info object for training data
  data_info di_train;
  di_train.n = n_train;
  di_train.R = R;
  di_train.p_cont = p_cont;
  di_train.p_cat = p_cat;
  di_train.p = p;
  di_train.z = tZ_train.begin();
  if(p_cont > 0) di_train.x_cont = tX_cont_train.begin();
  if(p_cat > 0) di_train.x_cat = tX_cat_train.begin();
  di_train.rp = residual;
  // END: creating data info object for training data
  
  // BEGIN: creating data info object for testing data (if present)
  data_info di_test;
  if(n_test > 0){
    di_test.n = n_test;
    di_test.R = R;
    di_test.p_cont = p_cont;
    di_test.p_cat = p_cat;
    di_test.p = p;
    di_test.z = tZ_test.begin();
    if(p_cont > 0) di_test.x_cont = tX_cont_test.begin();
    if(p_cat > 0)  di_test.x_cat = tX_cat_test.begin();
  }
  // END: creating data info object for testing data (if present)

    
  // BEGIN: creating tree prior info object
  std::vector<tree_prior_info> tree_pi_vec(R);
  for(int r = 0; r < R; ++r){
    // start by assigning stuff relevant to all ensembles
    if(p_cont > 0) tree_pi_vec[r].cutpoints = &cutpoints;
    if(p_cat > 0){
      tree_pi_vec[r].cat_levels = &cat_levels;
      tree_pi_vec[r].edges = &edges;
      tree_pi_vec[r].graph_cut_type = graph_cut_type;
    }
    
    tree_pi_vec[r].theta = &(theta[r]); // only use theta if nest_v_option = false
    tree_pi_vec[r].var_count = &(var_count[r]);
    tree_pi_vec[r].rule_count = &(rule_count[r]);
    
    tree_pi_vec[r].alpha = alpha_vec[r];
    tree_pi_vec[r].beta = beta_vec[r];
    tree_pi_vec[r].mu0 = mu0_vec[r];
    tree_pi_vec[r].tau = tau_vec[r];
  }
  // END: creating tree prior info object
  
  // BEGIN: initialize sigma
  double sigma = sigest;
  double total_sq_resid = 0.0; // sum of squared residuals
  double scale_post = 0.0;
  double nu_post = 0.0;
  // END : initialize sigma
  
  // BEGIN: initialize stuff for main MCMC loop
  int total_draws = 1 + burn + (nd-1)*thin;
  int sample_index = 0;
  int accept = 0;
  int total_accept = 0; // counts how many trees we change in each iteration
  // END: initialize stuff for main MCMC loop
  
  double tmp_mu; // for holding the value of mu when we're doing the backfitting
  
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
  

  // BEGIN: initialize suff_stat_maps and populate residuals
  for(int r = 0; r < R; ++r){
    for(int m = 0; m < M_vec[r]; ++m){
      // do an initial tree traversal
      // this is kind of silly when t is just a stump
      // but it may help if we were to allow start from an arbitrary ensemble
      tree_traversal(ss_train_vec[r][m], t_vec[r][m], di_train);
      
      // get the fit of each tree
      for(suff_stat_it l_it = ss_train_vec[r][m].begin(); l_it != ss_train_vec[r][m].end(); ++l_it){
        tmp_mu = t_vec[r][m].get_ptr(l_it->first)->get_mu(); // get the value of mu in the leaf
        if(l_it->second.size() > 0){
          for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
            allfit_train[*it] += di_train.z[r + (*it)*R] * tmp_mu;
            beta_fit_train[r + (*it) * R] += tmp_mu;
          } // closes loop over training obs in leaf
        } // closes if checking leaf has training obs
      } // closes loop over leafs
      
      if(n_test > 0){
        tree_traversal(ss_test_vec[r][m], t_vec[r][m], di_test);
        for(suff_stat_it l_it = ss_test_vec[r][m].begin(); l_it != ss_test_vec[r][m].end(); ++l_it){
          if(l_it->second.size() > 0){
            for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
              allfit_test[*it] += di_test.z[r + (*it)*R] * tmp_mu;
              beta_fit_test[r + (*it) * R] += tmp_mu;
            } // closes loop over testing obs in leaf
          } // closes if checking leaf has testing obs
        } // closes loop over leafs
      } // closes if checking for testing data
    } // closes loop over trees in ensemble
  } // closes loop over ensembles
  for(int i = 0; i < n_train; ++i) residual[i] = Y_train[i] - allfit_train[i];
  // BEGIN: initialize suff_stat_maps and populate residuals

  // BEGIN: create output containers
  arma::vec fit_train_mean = arma::zeros<arma::vec>(n_train); // posterior mean for training data
  arma::vec fit_test_mean = arma::zeros<arma::vec>(1); // posterior mean for testing data (if any)
  arma::mat beta_train_mean = arma::zeros<arma::mat>(n_train, R); // posterior mean fit for each ensemble for training data
  arma::mat beta_test_mean = arma::zeros<arma::mat>(1,1); // posterior mean for each ensemble for testing data (if any)
  if(n_test > 0){
    fit_test_mean.zeros(n_test);
    beta_test_mean.zeros(n_test, R);
  }
  
  arma::mat fit_train = arma::zeros<arma::mat>(1,1); // posterior samples for training data
  arma::mat fit_test = arma::zeros<arma::mat>(1,1); // posterior samples for testing data (if any)
  arma::cube beta_train = arma::zeros<arma::cube>(1,1,1);
  arma::cube beta_test = arma::zeros<arma::cube>(1,1,1);
  
  if(save_samples){
    // if we are saving all samples, then we resize the containers accordingly
    fit_train.zeros(nd, n_train);
    beta_train.zeros(nd, n_train, R);
    if(n_test > 0){
      fit_test.zeros(nd, n_test);
      beta_test.zeros(nd, n_test, R);
    }
  }
  
  arma::vec sigma_samples(total_draws);
  arma::cube var_count_samples(total_draws, p, R); // how many times was a variable used in each iteration

  arma::mat total_accept_samples(total_draws, R);
  Rcpp::List tree_draws(nd);
  // END: create output containers
  
  // BEGIN: main MCMC loop

  for(int iter = 0; iter < total_draws; ++iter){
    if(verbose){
      if(iter == 0) Rcpp::Rcout << "  MCMC Iteration: " << iter+1 << " of " << total_draws << "; Warmup" << std::endl;
      else if(iter % print_every == 0){
        Rcpp::checkUserInterrupt();
        Rcpp::Rcout << "  MCMC Iteration: " << iter << " of " << total_draws;
        if(iter < burn) Rcpp::Rcout << "; Warmup";
        else Rcpp::Rcout << "; Sampling";
        Rcpp::Rcout << std::endl;
      } else if(iter==total_draws-1){
        Rcpp::Rcout << "  MCMC Iteration: " << iter+1 << " of " << total_draws << "; Sampling" << std::endl;
      }
    }
    
    if(iter == burn && (n_test > 0)){
      // first sampling iteration, populate allfit_test and beta_fit_test
      for(int r = 0; r < R; ++r){
        for(int m = 0; m < M_vec[r]; ++m){
          for(suff_stat_it l_it = ss_test_vec[r][m].begin(); l_it != ss_test_vec[r][m].end(); ++l_it){
            tmp_mu = t_vec[r][m].get_ptr(l_it->first)->get_mu(); // get value of mu in leaf
            if(l_it->second.size() > 0){
              for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
                allfit_test[*it] += di_test.z[r + (*it) * R] * tmp_mu;
                beta_fit_test[r + (*it) * R] += tmp_mu;
              } // closes loop over observations in leaf
            } // closes if checking that leaf is non-empty
          } // closes loop over leafs in the tree
        } // closes loop over trees in the ensemble
      } // closes loop over ensembles
    } // finish populating allfit_test and beta_fit_test in first sampling iteration after warmup
    
    // BEGIN: update all regression trees
    for(int r = 0; r < R; ++r){
      total_accept = 0;
      
      // BEGIN: loop over all trees in a single ensmble
      for(int m = 0; m < M_vec[r]; ++m){
        // BEGIN: remove fit of a single tree (training)
        for(suff_stat_it l_it = ss_train_vec[r][m].begin(); l_it != ss_train_vec[r][m].end(); ++l_it){
          tmp_mu = t_vec[r][m].get_ptr(l_it->first)->get_mu(); // get the value of mu in the leaf
          if(l_it->second.size() > 0){
            for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
              allfit_train[*it] -= di_train.z[r + (*it) * R] * tmp_mu; // adjust the value of allfit
              beta_fit_train[r + (*it) * R] -= tmp_mu;
              residual[*it] += di_train.z[r + (*it) * R] * tmp_mu;
            }
          }
        } // closes loop over all leafs removing fit
        // END: remove fit of a single tree (training)
        /*
        // SKD: this step is a bit redundant: we can just update the tree
        // BEGIN: remove fit of a single tree (testing). Only occurs post-warmup
        if(iter >= burn && n_test > 0){
          for(suff_stat_it l_it = ss_test_vec[r][m].begin(); l_it != ss_test_vec[r][m].end(); ++l_it){
            tmp_mu = t_vec[r][m].get_ptr(l_it->first)->get_mu();
            if(l_it->second.size() > 0){
              for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
                allfit_test[*it] -= di_test.z[r + (*it) * R] * tmp_mu;
                beta_fit_test[r + (*it) * R] -= tmp_mu;
              }
            }
          }
        } // closes loop removing fit from test
        // END: remove fit of a single tree (testing). Only occurs post-warmup
*/
        // BEGIN: update the tree
        update_tree_unnested(t_vec[r][m], ss_train_vec[r][m], ss_test_vec[r][m], accept, r, sigma, di_train, di_test, tree_pi_vec[r], gen); // update the tree
        total_accept += accept;
        // END: update the tree
        // BEGIN: restore fit of updated tree (training)
        for(suff_stat_it l_it = ss_train_vec[r][m].begin(); l_it != ss_train_vec[r][m].end(); ++l_it){
          tmp_mu = t_vec[r][m].get_ptr(l_it->first)->get_mu();
          if(l_it->second.size() > 0){
            for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
              // add fit of m-th tree back to allfit and subtract it from the value of the residual
              allfit_train[*it] += di_train.z[r + (*it) * R] * tmp_mu;
              beta_fit_train[r + (*it) * R] += tmp_mu;
              residual[*it] -= di_train.z[r + (*it) * R] * tmp_mu;
            }
          }
        } // closes loop restoring fit for training
        // END: restore fit of updated tree (training)
        
        /*
        // BEGIN: restore fit of updated tree (testing). Only occurs post-warmup
        if(iter >= burn && n_test > 0){
          for(suff_stat_it l_it = ss_test_vec[r][m].begin(); l_it != ss_test_vec[r][m].end(); ++l_it){
            tmp_mu = t_vec[r][m].get_ptr(l_it->first)->get_mu();
            if(l_it->second.size() > 0){
              for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
                // add fit of m-th tree back to allfit and subtract it from the value of the residual
                allfit_test[*it] += di_test.z[r + (*it) * R] * tmp_mu;
                beta_fit_test[r + (*it) * R] += tmp_mu;
              }
            }
          } // closes loop restoring fit for training
        }// closes check to see if we need to update testing data fits
        // END: restore fit of updated tree (testing). Only occurs post-warmup
        */

      } // closes loop over all of the trees in this ensemble
      // END: loop over all trees in a single ensmble
      
      // BEGIN: update selection probabilities & save var_count & accept
      if(sparse) update_theta_u_subset(theta[r], u[r], var_count[r], a_u, b_u, gen);
      for(int j = 0; j < p; ++j) var_count_samples(iter, j, r) = var_count[r][j];
      total_accept_samples(iter, r) = total_accept; // how many trees changed in this iteration

    } // closes loop over the ensembles
    // END: update all regression trees
    
    // BEGIN: update sigma
    total_sq_resid = 0.0;
    for(int i = 0; i < n_train; ++i) total_sq_resid += pow(residual[i], 2.0); // sum of squared residuals
    scale_post = lambda * nu + total_sq_resid;
    nu_post = nu + ( (double) n_train);
    sigma = sqrt(scale_post/gen.chi_square(nu_post));
    sigma_samples(iter) = sigma;
    // END: update sigma
    
    // BEGIN: save post-warmup samples
    if( (iter >= burn) && ( (iter - burn)%thin == 0)){
      sample_index = (int) ( (iter-burn)/thin);
      // BEGIN: write tree to string
      if(save_trees){
        Rcpp::List tmp_tree_draws(R);
        for(int r = 0; r < R; ++r){
          Rcpp::CharacterVector tree_string_vec(M_vec[r]);
          for(int m = 0; m < M_vec[r]; ++m){
            tree_string_vec[m] = write_tree(t_vec[r][m], tree_pi_vec[r], set_str);
          }
          tmp_tree_draws[r] = tree_string_vec;
        }
        tree_draws[sample_index] = tmp_tree_draws;
      }
      // END: write tree to string
      
      // BEGIN: save samples and/or posterior means (training)
      if(save_samples){
        for(int i = 0; i < n_train; ++i){
          fit_train(sample_index,i) = allfit_train[i];
          fit_train_mean(i) += allfit_train[i];
          for(int r = 0; r < R; ++r){
            beta_train(sample_index, i, r) = beta_fit_train[r + i*R];
            beta_train_mean(i, r) += beta_fit_train[r + i * R];
          }
        }
      } else{
        for(int i = 0; i < n_train; ++i){
          fit_train_mean(i) += allfit_train[i];
          for(int r = 0; r < R; ++r){
            beta_train_mean(i, r) += beta_fit_train[r + i * R];
          }
        }
      }
      // END: save samples and/or posterior means (training)
      
      // BEGIN: save samples and/or posterior means (testing)
      if(n_test > 0){
        // reset the values
        for(int i = 0; i < n_test; ++i){
          allfit_test[i] = 0.0;
          for(int r = 0; r < R; ++r) beta_fit_test[r + i*R] = 0.0;
        }
        for(int r = 0; r < R; ++r){
          for(int m = 0; m < M_vec[r]; ++m){
            for(suff_stat_it ss_it = ss_test_vec[r][m].begin(); ss_it != ss_test_vec[r][m].end(); ++ss_it){
              tmp_mu = t_vec[r][m].get_ptr(ss_it->first)->get_mu();
              if(ss_it->second.size() > 0){
                for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
                  allfit_test[*it] += di_test.z[r + (*it) * R] * tmp_mu;
                  beta_fit_test[r + (*it) * R] += tmp_mu;
                }
              }
            }
          }
        }
        
        
        if(save_samples){
          for(int i = 0; i < n_test; ++i){
            fit_test(sample_index,i) = allfit_test[i];
            fit_test_mean(i) += allfit_test[i];
            for(int r = 0; r < R; ++r){
              beta_test(sample_index, i, r) = beta_fit_test[r + i*R];
              beta_test_mean(i, r) += beta_fit_test[r + i * R];
            }
          }
        } else{
          for(int i = 0; i < n_test; ++i){
            fit_test_mean(i) += allfit_test[i];
            for(int r = 0; r < R; ++r){
              beta_test_mean(i, r) += beta_fit_test[r + i * R];
            }
          }
        }
      }
      // END: save samples and/or posterior means (testing)
    } // closes if that checks whether we should save anything in this iteration
    // END: save post-warmup samples
  } // closes the main MCMC for loop
  // END: main MCMC loop
  
  // BEGIN: rescale posterior means
  fit_train_mean /= ( (double) nd);
  beta_train_mean /= ( (double) nd);
  if(n_test > 0){
    fit_test_mean /= ( (double) nd);
    beta_test_mean /= ( (double) nd);
  }
  // END: rescale posterior means
  
  
  // BEGIN: save output
  Rcpp::List results;
  results["fit_train_mean"] = fit_train_mean;
  results["beta_train_mean"] = beta_train_mean;
  if(save_samples){
    results["fit_train"] = fit_train;
    results["beta_train"] = beta_train;
  }
  if(n_test > 0){
    results["fit_test_mean"] = fit_test_mean;
    results["beta_test_mean"] = beta_test_mean;
    if(save_samples){
      results["fit_test"] = fit_test;
      results["beta_test"] = beta_test;
    }
  }
  results["sigma"] = sigma_samples;
  results["varcount"] = var_count_samples;
  
  results["total_accept"] = total_accept_samples;
  if(save_trees) results["trees"] = tree_draws;
  return results;
  
}


