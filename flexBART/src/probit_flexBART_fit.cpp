#include "update_tree.h"
#include "data_parsing_funs.h"
#include "funs.h"
// [[Rcpp::export(".probit_flexBART_fit")]]
Rcpp::List probit_flexBART_fit(Rcpp::IntegerVector Y_train,
                               Rcpp::NumericMatrix tX_cont_train,
                               Rcpp::IntegerMatrix tX_cat_train,
                               Rcpp::NumericMatrix tX_cont_test,
                               Rcpp::IntegerMatrix tX_cat_test,
                               Rcpp::LogicalVector unif_cuts,
                               Rcpp::Nullable<Rcpp::List> cutpoints_list,
                               Rcpp::Nullable<Rcpp::List> cat_levels_list,
                               Rcpp::Nullable<Rcpp::List> edge_mat_list,
                               Rcpp::LogicalVector graph_split, int graph_cut_type,
                               bool sparse, double a_u, double b_u,
                               double mu0, double tau,
                               int M,
                               int nd, int burn, int thin,
                               bool save_samples,
                               bool save_trees,
                               bool verbose, int print_every)
{
  Rcpp::RNGScope scope;
  RNG gen;
  
  set_str_conversion set_str; // for converting sets of integers into strings
  
  int n_train = 0;
  int n_test = 0;
  int p_cont = 0;
  int p_cat = 0;
  
  parse_training_data(n_train, p_cont, p_cat, tX_cont_train, tX_cat_train);
  if(Y_train.size() != n_train) Rcpp::stop("Number of observations in Y_train does not match number of rows in training design matrices");
  parse_testing_data(n_test, tX_cont_test, tX_cat_test, p_cat, p_cont);
  
  int p = p_cont + p_cat;
  
  if(verbose){
    Rcpp::Rcout << "n_train = " << n_train << " n_test = " << n_test;
    Rcpp::Rcout << " p_cont = " << p_cont << "  p_cat = " << p_cat << std::endl;
  }
  
  std::vector<std::set<double>> cutpoints;
  if(p_cont > 0){
    if(cutpoints_list.isNotNull()){
      Rcpp::List tmp_cutpoints  = Rcpp::List(cutpoints_list);
      parse_cutpoints(cutpoints, p_cont, tmp_cutpoints, unif_cuts);
    }
  }
  
  std::vector<std::set<int>> cat_levels;
  std::vector<int> K; // number of levels for the different categorical variables
  std::vector<std::vector<edge>> edges;
  
  if(p_cat > 0){
    if(cat_levels_list.isNotNull()){
      Rcpp::List tmp_cat_levels = Rcpp::List(cat_levels_list);
      parse_cat_levels(cat_levels, K, p_cat, tmp_cat_levels);
    } else{
      Rcpp::stop("Must provide list of categorical levels!");
    }
    if(edge_mat_list.isNotNull()){
      Rcpp::List tmp_edge_mat = Rcpp::List(edge_mat_list);
      parse_graphs(edges, p_cat, K, tmp_edge_mat, graph_split);
    }
  }
  
  double* allfit_train = new double[n_train];
  double* residual = new double[n_train];
  double* latents = new double[n_train]; // holds the latent variables
  
  std::vector<double> allfit_test;
  if(n_test > 0) allfit_test.resize(n_test);
  
  // set up the data info object for training data
  data_info di_train;
  di_train.n = n_train;
  di_train.p_cont = p_cont;
  di_train.p_cat = p_cat;
  di_train.p = p;
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
    if(p_cont > 0) di_test.x_cont = tX_cont_test.begin();
    if(p_cat > 0)  di_test.x_cat = tX_cat_test.begin();
  }
  
  // declare stuff for variable selection
  std::vector<double> theta(p, 1.0/ (double) p);
  double u = 1.0/(1.0 + (double) p);
  std::vector<int> var_count(p, 0); // count how many times a variable has been used in a splitting rule
  int rule_count = 0; // how many total decision rules are there in the ensemble

  tree_prior_info tree_pi;
  tree_pi.theta = &theta;
  tree_pi.var_count = &var_count;
  tree_pi.rule_count = &rule_count;
  
  if(p_cont > 0){
    tree_pi.unif_cuts = unif_cuts.begin(); // do we use uniform cutpoints?
    tree_pi.cutpoints = &cutpoints;
  }
  
  if(p_cat > 0){
    tree_pi.cat_levels = &cat_levels;
    tree_pi.edges = &edges;
    tree_pi.K = &K;
    tree_pi.graph_split = graph_split.begin();
    tree_pi.graph_cut_type = graph_cut_type;
  }
  tree_pi.mu0 = mu0;
  tree_pi.tau = tau;
  
  // stuff for MCMC loop
  int total_draws = 1 + burn + (nd-1)*thin;
  int sample_index = 0;
  int accept = 0;
  int total_accept = 0; // counts how many trees we change in each iteration
  rule_diag_t rule_diag;
  

  tree::npv bnv; // for checking that our ss map and our trees are not totally and utterly out of sync
  double tmp_mu; // for holding the value of mu when we're doing the backfitting
  
  // initialize the trees
  std::vector<tree> t_vec(M);
  std::vector<suff_stat> ss_train_vec(M);
  std::vector<suff_stat> ss_test_vec(M);
  
  for(int i = 0; i < n_train; i++) allfit_train[i] = 0.0;
  
  for(int m = 0; m < M; m++){
    // do an initial tree traversal
    // this is kind of silly when t is just a stump
    // but it may help if we were to allow start from an arbitrary ensemble
    tree_traversal(ss_train_vec[m], t_vec[m], di_train);
    
    // get the fit of each tree
    for(suff_stat_it ss_it = ss_train_vec[m].begin(); ss_it != ss_train_vec[m].end(); ++ss_it){
      tmp_mu = t_vec[m].get_ptr(ss_it->first)->get_mu(); // get the value of mu in the leaf
      if(ss_it->second.size() > 0){
        for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
          allfit_train[*it] += tmp_mu;
        }
      }
    }
    if(n_test > 0){
      tree_traversal(ss_test_vec[m], t_vec[m], di_test);
    }
  }
  
  // draw the initial set of latents and compute the residual
  for(int i = 0; i < n_train; i++){
    if(Y_train[i] == 1) latents[i] = gen.lo_trunc_norm(allfit_train[i], 0.0);
    else if(Y_train[i] == 0) latents[i] = gen.hi_trunc_norm(allfit_train[i],0.0);
    else{
      Rcpp::Rcout << " Outcome for observation i = " << i+1 << " is " << Y_train[i] << std::endl;
      Rcpp::stop("For probit regression, all outcomes must be 1 or 0.");
    }
    residual[i] = latents[i] - allfit_train[i];
  }
  
  double sigma = 1.0; // remember for probit, sigma is always 1!
  
  // output containers
  arma::vec fit_train_mean = arma::zeros<arma::vec>(n_train); // posterior mean for training data
  arma::vec fit_test_mean = arma::zeros<arma::vec>(1); // posterior mean for testing data (if any)
  if(n_test > 0) fit_test_mean.zeros(n_test); // set.size() can contain garbage values!!

  arma::mat fit_train = arma::zeros<arma::mat>(1,1); // posterior samples for training data
  arma::mat fit_test = arma::zeros<arma::mat>(1,1); // posterior samples for testing data (if any)
  if(save_samples){
    // if we are saving all samples, then we resize the containers accordingly
    fit_train.zeros(nd, n_train);
    if(n_test > 0) fit_test.zeros(nd, n_test);
  }
  
  arma::vec total_accept_samples(total_draws);
  arma::vec aa_proposed_samples(total_draws); // how many axis-aligned rules proposed in this iteration
  arma::vec aa_rejected_samples(total_draws); // how many axis-aligned rules rejected in this iteration
  arma::vec cat_proposed_samples(total_draws); // how many categorical rules proposed in this iteration
  arma::vec cat_rejected_samples(total_draws); // how many categorical rules rejected in this iteration

  arma::mat theta_samples(1,1); // unless we're doing DART, no need to waste space
  if(sparse) theta_samples.set_size(total_draws, p);
  arma::mat var_count_samples(total_draws, p); // always useful to see how often we're splitting on variables in the ensemble
  
  Rcpp::List tree_draws(nd);
  
  // main MCMC loop starts here
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
    // at the start of the loop, we need to update the latents
    for(int i = 0; i < n_train; i++){
      // residual is latent - allfit
      if(Y_train[i] == 1) latents[i] = gen.lo_trunc_norm(allfit_train[i], 0.0);
      else if(Y_train[i] == 0) latents[i] = gen.hi_trunc_norm(allfit_train[i], 0.0);
      else{
        Rcpp::Rcout << " Outcome for observation i = " << i+1 << " is " << Y_train[i] << std::endl;
        Rcpp::stop("For probit regression, all outcomes must be 1 or 0.");
      }
      residual[i] = latents[i] - allfit_train[i];
    }
    
    // loop over trees
    total_accept = 0;
    rule_diag.reset();
    for(int m = 0; m < M; m++){
      for(suff_stat_it ss_it = ss_train_vec[m].begin(); ss_it != ss_train_vec[m].end(); ++ss_it){
        // loop over the bottom nodes in m-th tree
        tmp_mu = t_vec[m].get_ptr(ss_it->first)->get_mu(); // get the value of mu in the leaf
        if(ss_it->second.size() > 0){
          for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
            // remove fit of m-th tree from allfit: allfit[i] -= tmp_mu
            // for partial residual: we could compute Y - allfit (now that allfit has fit of m-th tree removed)
            // numerically this is exactly equal to adding tmp_mu to the value of residual
            allfit_train[*it] -= tmp_mu; // adjust the value of allfit
            residual[*it] += tmp_mu;
          }
        }
      } // this whole loop is O(n)
      
      update_tree(t_vec[m], ss_train_vec[m], ss_test_vec[m], accept, rule_diag, sigma, di_train, di_test, tree_pi, gen); // update the tree
      total_accept += accept;
    
      // now we need to update the value of allfit
      for(suff_stat_it ss_it = ss_train_vec[m].begin(); ss_it != ss_train_vec[m].end(); ++ss_it){
        tmp_mu = t_vec[m].get_ptr(ss_it->first)->get_mu();
        if(ss_it->second.size() > 0){
          for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
            // add fit of m-th tree back to allfit and subtract it from the value of the residual
            allfit_train[*it] += tmp_mu;
            residual[*it] -= tmp_mu;
          }
        }
      } // this loop is also O(n)
    } // closes loop over all of the trees
   
    
    if(sparse){
      update_theta_u(theta, u, var_count, p, a_u, b_u, gen);
      for(int j = 0; j < p; j++){
        theta_samples(iter, j) = theta[j];
        var_count_samples(iter,j) = var_count[j];
      }
    } else{
      for(int j = 0; j < p; j++) var_count_samples(iter, j) = var_count[j];
    }

    // save information for the diagnostics
    total_accept_samples(iter) = total_accept; // how many trees changed in this iteration
    aa_proposed_samples(iter) = rule_diag.aa_prop;
    aa_rejected_samples(iter) = rule_diag.aa_rej;
    cat_proposed_samples(iter) = rule_diag.cat_prop;
    cat_rejected_samples(iter) = rule_diag.cat_rej;

    if( (iter >= burn) && ( (iter - burn)%thin == 0)){
      sample_index = (int) ( (iter-burn)/thin);
      // time to write each tree as a string
      if(save_trees){
        // Option 1 in coatless' answer:
        // https://stackoverflow.com/questions/37502121/assigning-rcpp-objects-into-an-rcpp-list-yields-duplicates-of-the-last-element
        Rcpp::CharacterVector tree_string_vec(M);
        for(int m = 0; m < M; m++){
          tree_string_vec[m] = write_tree(t_vec[m], tree_pi, set_str);
        }
        tree_draws[sample_index] = tree_string_vec; // dump a character vector holding each tree's draws into an element of an Rcpp::List
      }
      
      if(save_samples){
        for(int i = 0; i < n_train; i++){
          fit_train(sample_index,i) = R::pnorm(allfit_train[i], 0.0, 1.0, true, false);
          fit_train_mean(i) += R::pnorm(allfit_train[i], 0.0, 1.0, true, false);
        }
      } else{
        for(int i = 0; i < n_train; i++) fit_train_mean(i) += R::pnorm(allfit_train[i], 0.0, 1.0, true, false);
      }
      
      if(n_test > 0){
        for(int i = 0; i < n_test; i++) allfit_test[i] = 0.0; // reset the value of allfit_test
        
        for(int m = 0; m < M; m++){
          for(suff_stat_it ss_it = ss_test_vec[m].begin(); ss_it != ss_test_vec[m].end(); ++ss_it){
            tmp_mu = t_vec[m].get_ptr(ss_it->first)->get_mu(); // get the value of mu in the corresponding leaf
            if(ss_it->second.size() > 0){
              for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it) allfit_test[*it] += tmp_mu;
            }
          } // loop over the keys in the m-th sufficient stat map
        } // closes loop over trees
        
        if(save_samples){
          for(int i = 0; i < n_test; i++){
            fit_test(sample_index,i) = R::pnorm(allfit_test[i], 0.0, 1.0, true, false);
            fit_test_mean(i) += R::pnorm(allfit_test[i], 0.0, 1.0, true, false);
          }
        } else{
          for(int i = 0; i < n_test; i++) fit_test_mean(i) += R::pnorm(allfit_test[i], 0.0, 1.0, true, false);
        }
      } // closes if checking whether we have any test set observations
    } // closes if that checks whether we should save anything in this iteration
  } // closes the main MCMC for loop
  
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
  results["aa_proposed"] = aa_proposed_samples;
  results["aa_rejected"] = aa_rejected_samples;
  results["cat_proposed"] = cat_proposed_samples;
  results["cat_rejected"] = cat_rejected_samples;
  results["var_count"] = var_count_samples;
  if(save_trees) results["trees"] = tree_draws;
  if(sparse) results["theta"] = theta_samples;
  
  return results;
  
}
