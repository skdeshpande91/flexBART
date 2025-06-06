#include "update_tree.h"
#include "data_parsing_funs.h"
#include "funs.h"
// [[Rcpp::export("._single_fit_probit")]]
Rcpp::List single_probit_fit(Rcpp::IntegerVector Y_train,
                             Rcpp::IntegerMatrix cov_ensm,
                             Rcpp::NumericMatrix tX_cont_train,
                             Rcpp::IntegerMatrix tX_cat_train,
                             Rcpp::NumericMatrix tX_cont_test,
                             Rcpp::IntegerMatrix tX_cat_test,
                             Rcpp::Nullable<Rcpp::List> cutpoints_list,
                             Rcpp::Nullable<Rcpp::List> cat_levels_list,
                             Rcpp::Nullable<Rcpp::List> edge_mat_list,
                             Rcpp::Nullable<Rcpp::List> nest_list,
                             int graph_cut_type,
                             bool sparse, double a_u, double b_u,
                             bool nest_v, int nest_v_option, bool nest_c,
                             int M,
                             double alpha, double beta,
                             double mu0, double tau,
                             int nd, int burn, int thin,
                             bool save_samples,
                             bool save_trees,
                             bool verbose, int print_every)
{
  Rcpp::RNGScope scope;
  RNG gen;
  
  set_str_conversion set_str; // for converting sets of integers into strings
  
  // BEGIN: get dimensions of training data
  int n_train = Y_train.size(); // how many training observations
  int R = cov_ensm.cols();// or could hardcode to 1
  int p_cont = 0;
  int p_cat = 0;
  if(tX_cont_train.size() > 1) p_cont = tX_cont_train.rows();
  if(tX_cat_train.size() > 1) p_cat = tX_cat_train.rows();
  int p = p_cont + p_cat;
  int n_test = 0;
  if(p_cont > 0 && tX_cont_test.size() > 0) n_test = tX_cont_test.cols();
  else if(p_cat > 0 && tX_cat_test.size() > 0) n_test = tX_cat_test.cols();
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
  // declare stuff for variable selection
  std::vector<double> theta(p, 1.0/ (double) p);
  double u = 1.0/(1.0 + (double) p);
  std::vector<int> var_count(p, 0); // count how many times a variable has been used in a splitting rule
  int rule_count = 0; // how many total decision rules are there in the ensemble
  // END: create splitting probabilities
  
  // BEGIN: initialize containers for residuals and fit
  double* latent = new double[n_train];
  double* residual = new double[n_train];
  int tmp_n_test = 1;
  if(n_test > 0) tmp_n_test = n_test;
  double* tmp_fit_test = new double[tmp_n_test]; // for holding test set fits temporarily

  // END: initialize containers for residuals and fit
  
  
  // BEGIN: create data_info objects for training & testing
  data_info di_train;
  di_train.n = n_train;
  di_train.p_cont = p_cont;
  di_train.p_cat = p_cat;
  di_train.p = p;
  di_train.R = R;
  if(p_cont > 0) di_train.x_cont = tX_cont_train.begin();
  if(p_cat > 0) di_train.x_cat = tX_cat_train.begin();
  di_train.rp = residual;
  
  // set up the data info object for testing data
  data_info di_test;
  if(n_test > 0){
    di_test.n = n_test;
    di_test.p_cont = p_cont;
    di_test.p_cat = p_cat;
    di_test.p = p;
    di_test.R = R;
    if(p_cont > 0) di_test.x_cont = tX_cont_test.begin();
    if(p_cat > 0)  di_test.x_cat = tX_cat_test.begin();
  }
  // END: create data_info objects for training & testing
 
  // BEGIN: create tree prior info object
  tree_prior_info tree_pi;
  tree_pi.theta = &theta;
  tree_pi.var_count = &var_count;
  tree_pi.rule_count = &rule_count;
  tree_pi.nest_v = nest_v;
  tree_pi.nest_v_option = nest_v_option;
  tree_pi.nest_c = nest_c;
  
  if(p_cont > 0) tree_pi.cutpoints = &cutpoints;
  
  if(p_cat > 0){
    tree_pi.cat_levels = &cat_levels;
    tree_pi.edges = &edges;
    tree_pi.graph_cut_type = graph_cut_type;
    tree_pi.nesting = &nesting;
    tree_pi.nest_in = &(nest_graph_in[0]);
    tree_pi.nest_out = &(nest_graph_out[0]);
    tree_pi.nest_components = &(nest_graph_components[0]);
  }
  tree_pi.alpha = alpha;
  tree_pi.beta = beta;
  tree_pi.mu0 = mu0;
  tree_pi.tau = tau;
  // END: create tree prior info object
  
  // BEGIN: intialize sigma
  double sigma = 1.0; // for probit, sigma is fixed
  // END: initialize sigma
  
  // BEGIN: initialize stuff for main MCMC loop
  int total_draws = 1 + burn + (nd-1)*thin;
  int sample_index = 0;
  int accept = 0;
  int total_accept = 0; // counts how many trees we change in each iteration
  double tmp_fit; // for holding difference b/w latents and residuals
  double tmp_mu; // for holding the value of mu when we're doing the backfitting
  // END: initialize stuff for main MCMC loop
  
  // BEGIN: initialize tree vector & sufficient statistics maps
  std::vector<tree> t_vec(M);
  std::vector<suff_stat> ss_train_vec(M);
  std::vector<suff_stat> ss_test_vec(M);
  // END: initialize tree vector
  
  
  // BEGIN: initialize latents and residual
  double offset = R::qnorm(Rcpp::mean(Y_train), 0.0,1.0, true, false);
  for(int i = 0; i < n_train; ++i){
    if(Y_train[i] == 1) latent[i] = gen.lo_trunc_norm(offset, 0.0);
    else if(Y_train[i] == 0) latent[i] = gen.hi_trunc_norm(offset, 0.0);
    else{
      Rcpp::Rcout << " Outcome for observation i = " << i+1 << " is " << Y_train[i] << std::endl;
      Rcpp::stop("For probit regression, all outcomes must be 1 or 0.");
    }
    residual[i] = latent[i];
  }
  // END: initialize latents and residual

  // BEGIN: initializze tree vector and sufficient statistics maps
  for(int m = 0; m < M; ++m){
    tree_traversal(ss_train_vec[m], t_vec[m], di_train); // populates ss_train_vec[m]
    for(suff_stat_it l_it = ss_train_vec[m].begin(); l_it != ss_train_vec[m].end(); ++l_it){
      tmp_mu = t_vec[m].get_ptr(l_it->first)->get_mu(); // get the value of mu in the leaf
      if(l_it->second.size() > 0){
        for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it) residual[*it] -= tmp_mu; // in this initial sweep, we have to remove the fit of each tree.
      } // closes if checking leaf is non-empty and updating initial value of residual
    }
    if(n_test > 0) tree_traversal(ss_test_vec[m], t_vec[m], di_test);
  }
  // END: initialize tree vector & sufficient statistics maps

  // BEGIN: create output containers
  arma::vec fit_train_mean = arma::zeros<arma::vec>(n_train); // posterior mean for training data
  arma::vec fit_test_mean = arma::zeros<arma::vec>(1); // posterior mean for testing data (if any)
  if(n_test > 0) fit_test_mean.zeros(n_test); // arma::set.size can initialize with garbage values
  
  arma::mat fit_train = arma::zeros<arma::mat>(1,1); // posterior samples for training data
  arma::mat fit_test = arma::zeros<arma::mat>(1,1); // posterior samples for testing data (if any)
  if(save_samples){
    // if we are saving all samples, then we resize the containers accordingly
    fit_train.zeros(nd, n_train);
    if(n_test > 0) fit_test.zeros(nd, n_test);
  }
  
  
  arma::vec total_accept_samples(total_draws);
  arma::mat var_count_samples(nd, p);
  
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
    
    // BEGIN: update latents and residuals
    for(int i = 0; i < n_train; ++i){
      tmp_fit = latent[i] - residual[i]; // current fit for i-th observation
      if(Y_train[i] == 1) latent[i] = gen.lo_trunc_norm(tmp_fit, 0.0);
      else latent[i] = gen.hi_trunc_norm(tmp_fit,0.0);
      residual[i] = latent[i] - tmp_fit; // updates the residual
    }
    // END: update latents and residuals
    
    total_accept = 0;
    for(int m = 0; m < M; ++m){
      //BEGIN: remove fit of m-th tree
      for(suff_stat_it l_it = ss_train_vec[m].begin(); l_it != ss_train_vec[m].end(); ++l_it){
        // loop over the bottom nodes in m-th tree
        tmp_mu = t_vec[m].get_ptr(l_it->first)->get_mu(); // get the value of mu in the leaf
        if(l_it->second.size() > 0){
          for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it) residual[*it] += tmp_mu;
        } // closes if checking leaf is non-empty and computes partial residual
      } // closes loop over leafs
      // END: remove fit of m-th tree
      update_tree_single(t_vec[m], ss_train_vec[m], ss_test_vec[m], accept, sigma, di_train, di_test, tree_pi, gen); // update the tree
      total_accept += accept;
    
      // BEGIN: restore fit of m-th tree
      for(suff_stat_it l_it = ss_train_vec[m].begin(); l_it != ss_train_vec[m].end(); ++l_it){
        tmp_mu = t_vec[m].get_ptr(l_it->first)->get_mu();
        if(l_it->second.size() > 0){
          for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it) residual[*it] -= tmp_mu;
        } // closes if checking leaf is non-empty and computes full residual
      } // closes loop over leafs
      //END: restore fit of m-th tree
    } // closes loop over all of the trees
    
    // save information for the diagnostics
    total_accept_samples(iter) = total_accept; // how many trees changed in this iteration
    if(sparse) update_theta_u(theta, u, var_count, p, a_u, b_u, gen);
  } // closes burn-in
  // END: burn-in
  
  // BEGIN: post-burn-in
  for(int iter = burn; iter < total_draws; ++iter){
    if(iter==total_draws-1){
      if(verbose) Rcpp::Rcout << "  MCMC Iteration: " << iter+1 << " of " << total_draws << "; Sampling" << std::endl;
    } else if(iter%print_every == 0 || (iter==burn)){
      Rcpp::checkUserInterrupt();
      if(verbose) Rcpp::Rcout << "  MCMC Iteration: " << iter << " of " << total_draws << "; Sampling" << std::endl;
    }
    
    // BEGIN: update latents and residuals
    for(int i = 0; i < n_train; ++i){
      tmp_fit = latent[i] - residual[i]; // current fit for i-th observation
      if(Y_train[i] == 1) latent[i] = gen.lo_trunc_norm(tmp_fit, 0.0);
      else latent[i] = gen.hi_trunc_norm(tmp_fit,0.0);
      residual[i] = latent[i] - tmp_fit; // updates the residual
    }
    // END: update latents and residuals
    total_accept = 0;
    for(int m = 0; m < M; ++m){
      //BEGIN: remove fit of m-th tree
      for(suff_stat_it l_it = ss_train_vec[m].begin(); l_it != ss_train_vec[m].end(); ++l_it){
        // loop over the bottom nodes in m-th tree
        tmp_mu = t_vec[m].get_ptr(l_it->first)->get_mu(); // get the value of mu in the leaf
        if(l_it->second.size() > 0){
          for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it) residual[*it] += tmp_mu;
        } // closes if checking leaf is non-empty and computes partial residual
      } // closes loop over leafs
      // END: remove fit of m-th tree
      update_tree_single(t_vec[m], ss_train_vec[m], ss_test_vec[m], accept, sigma, di_train, di_test, tree_pi, gen); // update the tree
      total_accept += accept;
    
      // BEGIN: restore fit of m-th tree
      for(suff_stat_it l_it = ss_train_vec[m].begin(); l_it != ss_train_vec[m].end(); ++l_it){
        tmp_mu = t_vec[m].get_ptr(l_it->first)->get_mu();
        if(l_it->second.size() > 0){
          for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it) residual[*it] -= tmp_mu;
        } // closes if checking leaf is non-empty and computes full residual
      } // closes loop over leafs
      //END: restore fit of m-th tree
    } // closes loop over all of the trees
    
    // BEGIN: update theta (if sparse)
    if(sparse) update_theta_u(theta, u, var_count, p, a_u, b_u, gen);
    // END: update theta (if sparse)
    
    // save information for the diagnostics
    total_accept_samples(iter) = total_accept; // how many trees changed in this iteration
    
    if( (iter - burn)%thin == 0 ){
      sample_index = (int) ( (iter-burn)/thin);
      for(int j = 0; j < p; ++j) var_count_samples(sample_index,j) = var_count[j];

      // time to write each tree as a string
      if(save_trees){
        Rcpp::CharacterVector tree_string_vec(M);
        for(int m = 0; m < M; ++m) tree_string_vec[m] = write_tree(t_vec[m], tree_pi, set_str);
        tree_draws[sample_index] = tree_string_vec; // dump a character vector holding each tree's draws into an element of an Rcpp::List
      }
      if(save_samples){
        for(int i = 0; i < n_train; ++i){
          tmp_fit = latent[i] - residual[i];
          fit_train(sample_index,i) = R::pnorm(tmp_fit, 0.0, 1.0, true, false);
          fit_train_mean(i) += R::pnorm(tmp_fit, 0.0, 1.0, true, false);
        }
      } else{
        for(int i = 0; i < n_train; ++i){
          tmp_fit = latent[i] - residual[i];
          fit_train_mean(i) += R::pnorm(tmp_fit, 0.0, 1.0, true, false);
        }
      }
      if(n_test > 0){
        for(int i = 0; i < n_test; ++i) tmp_fit_test[i] = 0.0;
        
        for(int m = 0; m < M; ++m){
          for(suff_stat_it l_it = ss_test_vec[m].begin(); l_it != ss_test_vec[m].end(); ++l_it){
            tmp_mu = t_vec[m].get_ptr(l_it->first)->get_mu();
            if(l_it->second.size() > 0){
              for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it) tmp_fit_test[*it] += tmp_mu;
            } // closes if checking that leaf is non-empty and increment tmp_fit_test
          } // closes if/else checking whether we're saving samples or just posterior mean
        } // closes loop over trees
        
        if(save_samples){
          for(int i = 0; i < n_test; ++i){
            fit_test(sample_index, i) = R::pnorm(tmp_fit_test[i], 0.0, 1.0, true, false);
            fit_test_mean(i) += R::pnorm(tmp_fit_test[i], 0.0, 1.0, true, false);
          }
        } else{
          for(int i = 0; i < n_test; ++i) fit_test_mean(i) += R::pnorm(tmp_fit_test[i], 0.0, 1.0, true, false);
        } // closes if/else checking whether we're saving samples or just posterior mean
      } // close if checking that there are test set observations
    } // closes if that checks whether we should save anything in this iteration
  } // closes post-burn-in loop
  // END: post-burn-in
  
  fit_train_mean /= ( (double) nd);
  if(n_test > 0) fit_test_mean /= ( (double) nd);
  
  Rcpp::List results;
  
  results["fit_train_mean"] = fit_train_mean;
  if(save_samples) results["fit_train"] = fit_train;
  if(n_test > 0){
    results["fit_test_mean"] = fit_test_mean;
    if(save_samples) results["fit_test"] = fit_test;
  }
  results["total_accept"] = total_accept_samples;
  results["var_count"] = var_count_samples;
  if(save_trees) results["trees"] = tree_draws;

  
  return results;
  
}
