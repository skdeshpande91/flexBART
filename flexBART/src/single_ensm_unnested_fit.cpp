#include "update_tree_fast.h"
#include "data_parsing_funs.h"
#include "funs.h"
// [[Rcpp::export(".single_unnested_fit")]]
Rcpp::List flexBART_fit(Rcpp::NumericVector Y_train,
                        Rcpp::NumericMatrix tZ_train,
                        Rcpp::NumericMatrix tX_cont_train,
                        Rcpp::IntegerMatrix tX_cat_train,
                        Rcpp::NumericMatrix tZ_test,
                        Rcpp::NumericMatrix tX_cont_test,
                        Rcpp::IntegerMatrix tX_cat_test,
                        Rcpp::Nullable<Rcpp::List> cutpoints_list,
                        Rcpp::Nullable<Rcpp::List> cat_levels_list,
                        Rcpp::Nullable<Rcpp::List> edge_mat_list,
                        int graph_cut_type,
                        bool sparse, double a_u, double b_u,
                        double mu0, double tau,
                        double sigest, double lambda, double nu,
                        int M,
                        int nd, int burn, int thin,
                        bool save_samples,
                        bool save_trees,
                        bool verbose, int print_every)
{
  Rcpp::RNGScope scope;
  RNG gen;
  
  set_str_conversion set_str; // for converting sets of integers into strings
  
  // BEGIN: get dimensions of training data
  int n_train = tZ_train.cols(); // how many training observations
  int R = tZ_train.rows(); // how many ensembles
  if(R != 1) Rcpp::stop("[flexBART_fit]: this function assumes R = 1");
  int p_cont = 0;
  int p_cat = 0;
  if(tX_cont_train.size() > 1) p_cont = tX_cont_train.rows();
  if(tX_cat_train.size() > 1) p_cat = tX_cat_train.rows();
  int p = p_cont + p_cat;
  int n_test = 0;
  if(tZ_test.size() > 1) n_test = tZ_test.cols(); // how many test set observations
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

  // BEGIN: create splitting probabilities
  // declare stuff for variable selection
  std::vector<double> theta(p, 1.0/ (double) p);
  double u = 1.0/(1.0 + (double) p);
  std::vector<int> var_count(p, 0); // count how many times a variable has been used in a splitting rule
  int rule_count = 0; // how many total decision rules are there in the ensemble
  // END: create splitting probabilities
  
  // BEGIN: initialize containers for residuals and fit
  //double* allfit_train = new double[n_train];
  double* residual = new double[n_train];
  //int tmp_n_test = 1;
  //if(n_test > 0) tmp_n_test = n_test;
  //double* allfit_test = new double[tmp_n_test];
  //std::vector<double> allfit_test;
  //if(n_test > 0) allfit_test.resize(n_test);
  //for(int i = 0; i < n_train; ++i) allfit_train[i] = 0.0;
  //if(n_test > 0){
  //  for(int i = 0; i < n_test; ++i) allfit_test[i] = 0.0;
  //}
  // END: initialize containers for residuals and fit
  
  
  // BEGIN: create data_info objects for training & testing
  data_info di_train;
  di_train.n = n_train;
  di_train.p_cont = p_cont;
  di_train.p_cat = p_cat;
  di_train.p = p;
  di_train.R = R;
  di_train.z = tZ_train.begin();
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
    di_test.z = tZ_test.begin();
  }
  // END: create data_info objects for training & testing
 
  // BEGIN: create tree prior info object
  tree_prior_info tree_pi;
  tree_pi.theta = &theta;
  tree_pi.var_count = &var_count;
  tree_pi.rule_count = &rule_count;
  
  if(p_cont > 0) tree_pi.cutpoints = &cutpoints;
  
  if(p_cat > 0){
    tree_pi.cat_levels = &cat_levels;
    tree_pi.edges = &edges;
    tree_pi.graph_cut_type = graph_cut_type;
  }
  tree_pi.mu0 = mu0;
  tree_pi.tau = tau;
  // END: create tree prior info object
  
  // BEGIN: intialize sigma
  double sigma = sigest;
  double total_sq_resid = 0.0; // sum of squared residuals
  double scale_post = 0.0;
  double nu_post = 0.0;
  // END: initialize sigma
  
  // BEGIN: initialize stuff for main MCMC loop
  int total_draws = 1 + burn + (nd-1)*thin;
  int sample_index = 0;
  int accept = 0;
  int total_accept = 0; // counts how many trees we change in each iteration
  double tmp_mu; // for holding the value of mu when we're doing the backfitting
  // END: initialize stuff for main MCMC loop
  
  
  //tree::npv bnv; // for checking that our ss map and our trees are not totally and utterly out of sync
  
  // BEGIN: initialize tree vector & sufficient statistics maps
  std::vector<tree> t_vec(M);
  std::vector<suff_stat> ss_train_vec(M);
  std::vector<suff_stat> ss_test_vec(M);
  // END: initialize tree vector
  
  for(int i = 0; i < n_train; ++i) residual[i] = Y_train[i]; // start with all trees initialized at zero
  for(int m = 0; m < M; ++m){
    // do an initial tree traversal
    // this is kind of silly when t is just a stump
    // but it may help if we were to allow start from an arbitrary ensemble
    tree_traversal(ss_train_vec[m], t_vec[m], di_train);
    
    // get the fit of each tree
    for(suff_stat_it ss_it = ss_train_vec[m].begin(); ss_it != ss_train_vec[m].end(); ++ss_it){
      tmp_mu = t_vec[m].get_ptr(ss_it->first)->get_mu(); // get the value of mu in the leaf
      if(ss_it->second.size() > 0){
        for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
          //allfit_train[*it] += tmp_mu;
          residual[*it] -= tmp_mu; // in this initial sweep, we have to remove the fit of each tree.
        }
      }
    }
    if(n_test > 0) tree_traversal(ss_test_vec[m], t_vec[m], di_test);
  }
  
  //for(int i = 0; i < n_train; i++) residual[i] = Y_train[i] - allfit_train[i];
  
  
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
  
  arma::vec sigma_samples(total_draws);
  arma::vec total_accept_samples(total_draws);
  arma::mat var_count_samples(total_draws, p); // always useful to see how often we're splitting on variables in the ensemble
  
  Rcpp::List tree_draws(nd);
  // END: create output containers
  
  // main MCMC loop starts here
  int r = 0; // for this function only, assume there is only one
  for(int iter = 0; iter < total_draws; iter++){
    if(verbose){
    // remember that R is 1-indexed
      if( (iter < burn) && (iter % print_every == 0)){
        Rcpp::Rcout << "  MCMC Iteration: " << iter << " of " << total_draws << "; Warmup" << std::endl;
        Rcpp::checkUserInterrupt();
      } else if(((iter> burn) && (iter%print_every == 0)) || (iter == burn) ){
        Rcpp::Rcout << "  MCMC Iteration: " << iter << " of " << total_draws << "; Sampling" << std::endl;
        Rcpp::checkUserInterrupt();
      } else if( iter == total_draws-1){
        Rcpp::Rcout << "  MCMC Iteration: " << iter+1 << " of " << total_draws << "; Sampling" << std::endl;
        Rcpp::checkUserInterrupt();
      }
    }
    
    // loop over trees
    total_accept = 0;
    for(int m = 0; m < M; ++m){
      for(suff_stat_it ss_it = ss_train_vec[m].begin(); ss_it != ss_train_vec[m].end(); ++ss_it){
        // loop over the bottom nodes in m-th tree
        tmp_mu = t_vec[m].get_ptr(ss_it->first)->get_mu(); // get the value of mu in the leaf
        for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
          // remove fit of m-th tree from allfit: allfit[i] -= tmp_mu
          // for partial residual: we could compute Y - allfit (now that allfit has fit of m-th tree removed)
          // numerically this is exactly equal to adding tmp_mu to the value of residual
          //allfit_train[*it] -= tmp_mu; // adjust the value of allfit
          residual[*it] += tmp_mu;
        }
      } // this whole loop is O(n)
      
      update_tree_unnested(t_vec[m], ss_train_vec[m], ss_test_vec[m], accept, r, sigma, di_train, di_test, tree_pi, gen); // update the tree
      total_accept += accept;
    
      // now we need to update the value of allfit
      for(suff_stat_it ss_it = ss_train_vec[m].begin(); ss_it != ss_train_vec[m].end(); ++ss_it){
        tmp_mu = t_vec[m].get_ptr(ss_it->first)->get_mu();
        if(ss_it->second.size() > 0){
          for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
            // add fit of m-th tree back to allfit and subtract it from the value of the residual
            //allfit_train[*it] += tmp_mu;
            residual[*it] -= tmp_mu;
          }
        }
      } // this loop is also O(n)
    } // closes loop over all of the trees
    // ready to update sigma
    total_sq_resid = 0.0;
    for(int i = 0; i < n_train; i++) total_sq_resid += pow(residual[i], 2.0); // sum of squared residuals
    
    scale_post = lambda * nu + total_sq_resid;
    nu_post = nu + ( (double) n_train);
    sigma = sqrt(scale_post/gen.chi_square(nu_post));
    sigma_samples(iter) = sigma;
    
    // save information for the diagnostics
    total_accept_samples(iter) = total_accept; // how many trees changed in this iteration
    if(sparse) update_theta_u(theta, u, var_count, p, a_u, b_u, gen);
    for(int j = 0; j < p; ++j) var_count_samples(iter,j) = var_count[j];
        
    if( (iter >= burn) && ( (iter - burn)%thin == 0)){
      sample_index = (int) ( (iter-burn)/thin);
      // time to write each tree as a string
      if(save_trees){
        // Option 1 in coatless' answer:
        // https://stackoverflow.com/questions/37502121/assigning-rcpp-objects-into-an-rcpp-list-yields-duplicates-of-the-last-element
        Rcpp::CharacterVector tree_string_vec(M);
        for(int m = 0; m < M; ++m) tree_string_vec[m] = write_tree(t_vec[m], tree_pi, set_str);
        tree_draws[sample_index] = tree_string_vec; // dump a character vector holding each tree's draws into an element of an Rcpp::List
      }
      double tmp_allfit = 0.0;
      if(save_samples){
        for(int i = 0; i < n_train; ++i){
          tmp_allfit = Y_train[i] - residual[i];
          fit_train(sample_index,i) = tmp_allfit;
          fit_train_mean(i) += tmp_allfit;
        }
      } else{
        for(int i = 0; i < n_train; ++i){
          fit_train_mean(i) += Y_train[i] - residual[i];
        }
      }

      if(n_test > 0){
        if(save_samples){
          for(int m = 0; m < M; ++m){
            for(suff_stat_it ss_it = ss_test_vec[m].begin(); ss_it != ss_test_vec[m].end(); ++ss_it){
              tmp_mu = t_vec[m].get_ptr(ss_it->first)->get_mu();
              if(ss_it->second.size() > 0){
                for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
                  fit_test(sample_index, *it) += tmp_mu;
                  fit_test_mean(*it) += tmp_mu;
                }
              }
            }
          }
        } else{
          for(int m = 0; m < M; ++m){
            for(suff_stat_it ss_it = ss_test_vec[m].begin(); ss_it != ss_test_vec[m].end(); ++ss_it){
              tmp_mu = t_vec[m].get_ptr(ss_it->first)->get_mu();
              if(ss_it->second.size() > 0){
                for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
                  fit_test_mean(*it) += tmp_mu;
                }
              }
            }
          }
        }
      }
      

    } // closes if that checks whether we should save anything in this iteration
  } // closes the main MCMC for loop

  fit_train_mean /= ( (double) nd);
  if(n_test > 0) fit_test_mean /= ( (double) nd);
  
  Rcpp::List results;
  
  results["fit_train_mean"] = fit_train_mean;
  if(save_samples){
    results["fit_train"] = fit_train;
  }
  if(n_test > 0){
    results["fit_test_mean"] = fit_test_mean;
    if(save_samples){
      results["fit_test"] = fit_test;
    }
  }
  results["sigma"] = sigma_samples;
  results["total_accept"] = total_accept_samples;
  results["var_count"] = var_count_samples;
  if(save_trees){
    results["trees"] = tree_draws;
  }

  
  return results;
  
}
