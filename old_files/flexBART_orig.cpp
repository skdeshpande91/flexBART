//  flex_bart_fast.cpp

#include "update_tree.h"

Rcpp::List flexBART_fit(Rcpp::NumericVector Y_train,
                        Rcpp::NumericMatrix tX_cont_train,
                        Rcpp::IntegerMatrix tX_cat_train,
                        Rcpp::NumericMatrix tX_cont_test,
                        Rcpp::IntegerMatrix tX_cat_test,
                        Rcpp::LogicalVector unif_cuts,
                        Rcpp::Nullable<Rcpp::List> cutpoints_list,
                        Rcpp::Nullable<Rcpp::List> cat_levels_list,
                        Rcpp::LogicalVector graph_split, int graph_cut_type,
                        Rcpp::Nullable<Rcpp::List> adj_support_list,
                        bool rc_split, double prob_rc, double a_rc, double b_rc,
                        bool sparse, double a_u, double b_u,
                        double mu0, double tau,
                        double lambda, double nu,
                        int M,
                        int nd, int burn, int thin,
                        bool save_trees,
                        bool verbose, int print_every,
                        bool check_ss_map = false) // check that our sufficient stat map has the right keys
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
  std::vector<std::set<int>> cat_levels;
  std::vector<int> K; // number of levels for the different categorical variables
  std::vector<std::vector<unsigned int>> adj_support;
  
  if(p_cont > 0){
    if(cutpoints_list.isNotNull()){
      Rcpp::List tmp_cutpoints = Rcpp::List(cutpoints_list);
      parse_cutpoints(cutpoints, p_cont, tmp_cutpoints, unif_cuts);
    }
  }
  
  if(p_cat > 0){
    if(cat_levels_list.isNotNull()){
      Rcpp::List tmp_cat_levels = Rcpp::List(cat_levels_list);
      parse_cat_levels(cat_levels, K, p_cat, tmp_cat_levels);
    } else{
      Rcpp::stop("Must provide categorical levels.");
    }
    if(adj_support_list.isNotNull()){
      Rcpp::List tmp_adj_support = Rcpp::List(adj_support_list);
      parse_cat_adj(adj_support, p_cat, tmp_adj_support, graph_split);
    }
  }
  
  double* allfit_train = new double[n_train];
  double* residual = new double[n_train];
  
  std::vector<double> allfit_test;
  if(n_test > 0) allfit_test.resize(n_test); // can wrap this as an Rcpp::NumericVector later on
  
  // set up our data info object
  data_info di_train;
  di_train.n = n_train;
  di_train.p_cont = p_cont;
  di_train.p_cat = p_cat;
  di_train.p = p;
  if(p_cont > 0) di_train.x_cont = tX_cont_train.begin();
  if(p_cat > 0){
    di_train.x_cat = tX_cat_train.begin();
    di_train.cat_levels = &cat_levels;
    di_train.K = &K;
    di_train.adj_support = &adj_support;
  }
  di_train.rp = residual;
  
  data_info di_test;
  if(n_test > 0){
    di_test.n = n_test;
    di_test.p_cont = p_cont;
    di_test.p_cat = p_cat;
    di_test.p = p;
    if(p_cont > 0) di_test.x_cont = tX_cont_test.begin();
    if(p_cat > 0){
      di_test.x_cat = tX_cat_test.begin();
      di_test.cat_levels = &cat_levels;
      di_test.K = &K;
      di_test.adj_support = &adj_support;
    }
  }
  
  // stuff for variable selection
  std::vector<double> theta(p, 1.0/ (double) p);
  double u = 1.0/(1.0 + (double) p);
  std::vector<int> var_count(p, 0); // count how many times a variable has been used in a splitting rule
  int rule_count = 0; // how many total decision rules are there in the ensemble
  int rc_rule_count = 0; // how many random combination rules are there in the ensemble
  int rc_var_count = 0; // only when we are using random combination rules
  double theta_rc = 0.0; // prob of including a variable in a random combination rule
  if(p_cont >= 2 && rc_split){
    theta_rc = 2.0/( (double) p_cont);
  }
  
  tree_prior_info tree_pi;
  tree_pi.theta = &theta;
  tree_pi.var_count = &var_count;
  tree_pi.rule_count = &rule_count;
  
  tree_pi.unif_cuts = unif_cuts.begin(); // do we use uniform cutpoints?
  tree_pi.cutpoints = &cutpoints;

  
  tree_pi.graph_split = graph_split.begin();
  tree_pi.graph_cut_type = graph_cut_type;
  tree_pi.rc_split = rc_split;
  tree_pi.prob_rc = &prob_rc;
  tree_pi.theta_rc = &theta_rc;
  tree_pi.rc_var_count = &rc_var_count;
  tree_pi.rc_rule_count = &rc_rule_count;
  tree_pi.mu0 = mu0;
  tree_pi.tau = tau;
  
  // stuff for sigma
  double sigma = 1.0;
  double total_sq_resid = 0.0; // sum of squared residuals
  double scale_post = 0.0;
  double nu_post = 0.0;
  
  // stuff for MCMC loop
  int total_draws = 1 + burn + (nd-1)*thin;
  int sample_index = 0;
  int accept = 0;
  int total_accept = 0; // counts how many trees we change in each iteration
  tree::npv bnv; // for checking that our ss map and our trees are not totally and utterly out of sync
  double tmp_mu; // for holding the value of mu when we're doing the backfitting
  
  // initialize the trees
  std::vector<tree> t_vec(M);
  std::vector<suff_stat> ss_vec(M);
  
  for(int i = 0; i < n_train; i++) allfit_train[i] = 0.0;
  
  for(int m = 0; m < M; m++){
    // do an initial tree traversal
    // this is kind of silly when t is just a stump
    // but it may help if we were to allow start from an arbitrary ensemble
    tree_traversal(ss_vec[m], t_vec[m], di_train);
    
    // get the fit of each tree
    for(suff_stat_it ss_it = ss_vec[m].begin(); ss_it != ss_vec[m].end(); ++ss_it){
      tmp_mu = t_vec[m].get_ptr(ss_it->first)->get_mu(); // get the value of mu in the leaf
      for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
        allfit_train[*it] += tmp_mu;
      }
    }
  }
  for(int i = 0; i < n_train; i++) residual[i] = Y_train[i] - allfit_train[i];
  
  // output containers
  
  //arma::mat fit_train = arma::zeros<arma::mat>(n_train, nd);
  arma::mat fit_train = arma::zeros<arma::mat>(nd, n_train); // similar to regular BART
  arma::mat fit_test = arma::zeros<arma::mat>(1,1);
  //if(n_test > 0) fit_test.set_size(n_test, nd);
  if(n_test > 0) fit_test.set_size(nd, n_test);
  arma::vec sigma_samples(total_draws);
  arma::vec total_accept_samples(nd);
  arma::mat theta_samples(1,1); // unless we're doing DART, no need to waste space
  if(sparse) theta_samples.set_size(total_draws, p);
  arma::mat var_count_samples(total_draws, p); // always useful to see how often we're splitting on variables in the ensemble
  arma::vec theta_rc_samples(1);
  arma::vec rc_rule_count_samples(1);
  arma::vec rc_var_count_samples(1);
  if(rc_split){
    theta_rc_samples.set_size(total_draws);
    rc_rule_count_samples.set_size(total_draws);
    rc_var_count_samples.set_size(total_draws);
  }
  
  Rcpp::List tree_draws(nd);
  
  // main MCMC loop goes here
  for(int iter = 0; iter < total_draws; iter++){
    if(verbose){
      if( (iter < burn) && (iter % print_every == 0)){
        Rcpp::Rcout << "  MCMC Iteration: " << iter << " of " << total_draws << "; Warmup" << std::endl;
        Rcpp::checkUserInterrupt();
      } else if(((iter> burn) && (iter%print_every == 0)) || (iter == burn)){
        Rcpp::Rcout << "  MCMC Iteration: " << iter << " of " << total_draws << "; Sampling" << std::endl;
        Rcpp::checkUserInterrupt();
      }
    }
    // loop over trees
    total_accept = 0;
    for(int m = 0; m < M; m++){
      //Rcpp::Rcout << m << " "; // DEBUG
      for(suff_stat_it ss_it = ss_vec[m].begin(); ss_it != ss_vec[m].end(); ++ss_it){
        // loop over the bottom nodes in m-th tree
        tmp_mu = t_vec[m].get_ptr(ss_it->first)->get_mu(); // get the value of mu in the leaf
        for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
          // remove fit of m-th tree from allfit: allfit[i] -= tmp_mu
          // for partial residual: we could compute Y - allfit (now that allfit has fit of m-th tree removed)
          // numerically this is exactly equal to adding tmp_mu to the value of residual
          allfit_train[*it] -= tmp_mu; // adjust the value of allfit
          residual[*it] += tmp_mu;
        }
      } // this whole loop is O(n)
      
      update_tree(t_vec[m], ss_vec[m], accept, sigma, di_train, tree_pi, gen); // update the tree!
      total_accept += accept;
    
      if(check_ss_map){
        // this checks whether our sufficient stat map is correct.
        bnv.clear();
        t_vec[m].get_bots(bnv); // get the bottom nodes of the tree
        if(bnv.size() != ss_vec[m].size()){
          Rcpp::Rcout << "# leafs in tree = " << bnv.size() << " # elements in suff stat map = " << ss_vec[m].size() << std::endl;
          t_vec[m].print();
          Rcpp::Rcout << "suff stat map. nid(# observations):";
          for(suff_stat_it it = ss_vec[m].begin(); it != ss_vec[m].end(); ++it){
            Rcpp::Rcout << " " << it->first << "(" << it->second.size() << ")";
          }
          Rcpp::Rcout << std::endl;
          Rcpp::stop("number of leaves and number of elements in suff stat map must be equal!");
        } else{
          // ss map has the right number of elements
          for(tree::npv_it bn_it = bnv.begin(); bn_it != bnv.end(); ++bn_it){
            if(ss_vec[m].count( (*bn_it)->get_nid() ) != 1){
              // could not find an element in ss with key equal to the nid of a bottom node
              Rcpp::Rcout << "bottom node nid " << (*bn_it)->get_nid() << " is not a key in suff stat map!" << std::endl;
              Rcpp::stop("there is a mismatch between bottom node id's and the keys of the suff stat map!");
            }
          }
        }
      } // closes if for checking suff_stat_map and tree are in sync.
      
      // now we need to update the value of allfit
      for(suff_stat_it ss_it = ss_vec[m].begin(); ss_it != ss_vec[m].end(); ++ss_it){
        tmp_mu = t_vec[m].get_ptr(ss_it->first)->get_mu();
        for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
          // add fit of m-th tree back to allfit and subtract it from the value of the residual
          allfit_train[*it] += tmp_mu;
          residual[*it] -= tmp_mu;
        }
      } // this loop is also O(n)
    } // closes loop over all of the trees
    //Rcpp::Rcout << std::endl; // DEBUG
    // ready to update sigma
    total_sq_resid = 0.0;
    for(int i = 0; i < n_train; i++) total_sq_resid += pow(residual[i], 2.0); // sum of squared residuals
    
    scale_post = lambda * nu + total_sq_resid;
    nu_post = nu + ( (double) n_train);
    sigma = sqrt(scale_post/gen.chi_square(nu_post));
    sigma_samples(iter) = sigma;
    
    if(sparse){
      update_theta_u(theta, u, var_count, p, a_u, b_u, gen);
      for(int j = 0; j < p; j++){
        theta_samples(iter, j) = theta[j];
        var_count_samples(iter,j) = var_count[j];
      }
    } else{
      for(int j = 0; j < p; j++) var_count_samples(iter, j) = var_count[j];
    }
    if(rc_split){
      update_theta_rc(theta_rc, rc_var_count, rc_rule_count, a_rc, b_rc, p_cont, gen);
      theta_rc_samples(iter) = theta_rc;
      rc_var_count_samples(iter) = rc_var_count;
      rc_rule_count_samples(iter) = rc_rule_count;
    }
    
    
    if( (iter >= burn) && ( (iter - burn)%thin == 0)){
      sample_index = (int) ( (iter-burn)/thin);
      total_accept_samples(sample_index) = total_accept; // how many trees changed in this iteration
      // time to write each tree as a string
      if(save_trees){
        // Option 1 in coatless' answer:
        // https://stackoverflow.com/questions/37502121/assigning-rcpp-objects-into-an-rcpp-list-yields-duplicates-of-the-last-element
        Rcpp::CharacterVector tree_string_vec(M);
        for(int m = 0; m < M; m++){
          tree_string_vec[m] = write_tree(t_vec[m], di_train, set_str);
        }
        tree_draws[sample_index] = tree_string_vec; // dump a character vector holding each tree's draws into an element of an Rcpp::List
      }
      
      for(int i = 0; i < n_train; i++) fit_train(sample_index, i) = allfit_train[i];
      if(n_test > 0){
        fit_ensemble(allfit_test, t_vec, di_test);
        for(int i = 0; i < n_test; i++) fit_test(sample_index, i) = allfit_test[i];
      }
    } // closes if that checks whether we should save anything in this iteration
  } // closes the main MCMC for loop

  Rcpp::List results;
  results["fit_train"] = fit_train;
  if(n_test > 0) results["fit_test"] = fit_test;
  results["sigma"] = sigma_samples;
  results["total_accept"] = total_accept_samples;
  results["var_count"] = var_count_samples;
  if(save_trees) results["trees"] = tree_draws;
  if(sparse) results["theta"] = theta_samples;
  if(rc_split){
    results["theta_rc_samples"] = theta_rc_samples;
    results["rc_var_count_samples"] = rc_var_count_samples;
    results["rc_rule_count_samples"] = rc_rule_count_samples;
  }
  return results;
}
