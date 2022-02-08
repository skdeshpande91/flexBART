//  flex_bart_fast.cpp

#include "update_tree.h"

// [[Rcpp::export(".flexBARTfast_fit")]]
Rcpp::List flexBARTfast_fit(Rcpp::NumericVector Y_train,
                            Rcpp::NumericMatrix tX_cont_train,
                            Rcpp::IntegerMatrix tX_cat_train,
                            Rcpp::NumericMatrix tX_cont_test,
                            Rcpp::IntegerMatrix tX_cat_test,
                            Rcpp::Nullable<Rcpp::List> cutpoints_list,
                            Rcpp::Nullable<Rcpp::List> cat_levels_list,
                            Rcpp::Nullable<Rcpp::List> adj_support_list,
                            bool unif_cuts,
                            bool mst_split, bool mst_reweight,
                            double mu0, double tau,
                            double prob_aa, double prob_rc, // probs of axis-aligned and random-combination splits
                            double lambda, double nu,
                            int M = 50,
                            int nd = 1000, int burn = 1000, int thin = 1,
                            bool save_trees = true,
                            bool verbose = true, int print_every = 50,
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
  std::vector<int> K; // number of levels for each categorical variable
  std::vector<std::vector<unsigned int>> adj_support;
  
  if(p_cont > 0){
    if(cutpoints_list.isNotNull() && !unif_cuts){
      Rcpp::List tmp_cutpoints = Rcpp::List(cutpoints_list);
      parse_cutpoints(cutpoints, p_cont, tmp_cutpoints);
    }
  }
  
  if(p_cat > 0){
    if(cat_levels_list.isNotNull() && adj_support_list.isNotNull()){
      Rcpp::List tmp_cat_levels = Rcpp::List(cat_levels_list);
      Rcpp::List tmp_adj_support = Rcpp::List(adj_support_list);
      parse_categorical(cat_levels, adj_support, K, p_cat, tmp_cat_levels, tmp_adj_support);
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
  di_train.unif_cuts = unif_cuts; // do we use uniform cutpoints?
  if(p_cont > 0){
    di_train.x_cont = tX_cont_train.begin();
    di_train.cutpoints = &cutpoints;
  }
  if(p_cat > 0){
    di_train.x_cat = tX_cat_train.begin();
    di_train.cat_levels = &cat_levels;
    di_train.K = &K;
    di_train.adj_support = &adj_support;
  }
  di_train.rp = residual;
  
  data_info di_test;
  di_test.n = n_test;
  if(n_test > 0){
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
  
  tree_prior_info tree_pi;
  if(p_cont == 0){
    // no continuous variables so no axis-aligned or random combination rules
    tree_pi.prob_aa = 0.0;
    tree_pi.prob_rc = 0.0;
  } else{
    if(p_cat == 0){
      // prob_rc + prob_aa need to sum to 1
      tree_pi.prob_aa = prob_aa/(prob_aa + prob_rc);
      tree_pi.prob_rc = prob_rc/(prob_aa + prob_rc);
    } else{
      tree_pi.prob_aa = prob_aa;
      tree_pi.prob_rc = prob_rc;
    }
  }
  
  tree_pi.mst_split = mst_split;
  tree_pi.mst_reweight = mst_reweight;
  
  // stuff for variable selection
  int aa_rule_count = 0;
  int rc_rule_count = 0;
  int cat_rule_count = 0;
  
  tree_pi.aa_rule_count = &aa_rule_count;
  tree_pi.rc_rule_count = &rc_rule_count;
  tree_pi.cat_rule_count = &cat_rule_count;
    
  std::vector<int> aa_var_count;
  int rc_var_count = 0;
  std::vector<int> cat_var_count;
  
  std::vector<double> theta_aa;
  double theta_rc = 0.0;
  std::vector<double> theta_cat;
  
  if(p_cont > 0){
    aa_var_count.resize(p_cont, 0);
    rc_var_count = 0;
    theta_aa.resize(p_cont, 1.0/( (double) p_cont));
    theta_rc = 2.0/( (double) p_cont);
    
    tree_pi.theta_aa = &theta_aa;
    tree_pi.theta_rc = &theta_rc;
    tree_pi.aa_var_count = &aa_var_count;
    tree_pi.rc_var_count = &rc_var_count;
    
  }
  if(p_cat > 0){
    cat_var_count.resize(p_cat, 0);
    theta_cat.resize(p_cat, 1.0/( (double) p_cat));
    
    tree_pi.cat_var_count = &cat_var_count;
    tree_pi.theta_cat = &theta_cat;
  }
 
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
      for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
        allfit_train[*it] += tmp_mu;
      }
    }
    
    if(n_test > 0){
      tree_traversal(ss_test_vec[m], t_vec[m], di_test);
    }
    
    
  }
  for(int i = 0; i < n_train; i++) residual[i] = Y_train[i] - allfit_train[i];
  
  // output containers
  
  arma::mat fit_train = arma::zeros<arma::mat>(n_train, nd);
  arma::mat fit_test = arma::zeros<arma::mat>(1,1);
  if(n_test > 0) fit_test.set_size(n_test, nd);
  arma::vec sigma_samples(total_draws);
  arma::vec total_accept_samples(nd);
  
  Rcpp::List tree_draws(nd);
  
  // main MCMC loop goes here
  for(int iter = 0; iter < total_draws; iter++){
    if( (iter < burn) && (iter % print_every == 0)){
      Rcpp::Rcout << "  MCMC Iteration: " << iter << " of " << total_draws << "; Burn-in" << std::endl;
      Rcpp::checkUserInterrupt();
    } else if(((iter> burn) && (iter%print_every == 0)) || (iter == burn)){
      Rcpp::Rcout << "  MCMC Iteration: " << iter << " of " << total_draws << "; Sampling" << std::endl;
      Rcpp::checkUserInterrupt();
    }
    
    // loop over trees
    total_accept = 0;
    for(int m = 0; m < M; m++){
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
      
      //update_tree(t_vec[m], ss_vec[m], accept, sigma, di_train, tree_pi, gen); // update the tree!
      update_tree(t_vec[m], ss_train_vec[m], ss_test_vec[m], accept, sigma, di_train, di_test, tree_pi, gen);
      total_accept += accept;
    
      if(check_ss_map){
        // this checks whether our sufficient stat map is correct.
        bnv.clear();
        t_vec[m].get_bots(bnv); // get the bottom nodes of the tree
        if(bnv.size() != ss_train_vec[m].size()){
          Rcpp::Rcout << "# leafs in tree = " << bnv.size() << " # elements in training suff stat map = " << ss_train_vec[m].size() << std::endl;
          t_vec[m].print();
          Rcpp::Rcout << "suff stat map. nid(# observations):";
          for(suff_stat_it it = ss_train_vec[m].begin(); it != ss_train_vec[m].end(); ++it){
            Rcpp::Rcout << " " << it->first << "(" << it->second.size() << ")";
          }
          Rcpp::Rcout << std::endl;
          Rcpp::stop("number of leaves and number of elements in suff stat map must be equal!");
        } else{
          // ss map has the right number of elements
          for(tree::npv_it bn_it = bnv.begin(); bn_it != bnv.end(); ++bn_it){
            if(ss_train_vec[m].count( (*bn_it)->get_nid() ) != 1){
              // could not find an element in ss with key equal to the nid of a bottom node
              Rcpp::Rcout << "bottom node nid " << (*bn_it)->get_nid() << " is not a key in training suff stat map!" << std::endl;
              Rcpp::stop("there is a mismatch between bottom node id's and the keys of the training suff stat map!");
            }
          }
        }
        if(n_test > 0){
          if(bnv.size() != ss_test_vec[m].size()){
            Rcpp::Rcout << "# leafs in tree = " << bnv.size() << " # elements in training suff stat map = " << ss_test_vec[m].size() << std::endl;
            t_vec[m].print();
            Rcpp::Rcout << "suff stat map. nid(# observations):";
            for(suff_stat_it it = ss_test_vec[m].begin(); it != ss_test_vec[m].end(); ++it){
              Rcpp::Rcout << " " << it->first << "(" << it->second.size() << ")";
            }
            Rcpp::Rcout << std::endl;
            Rcpp::stop("number of leaves and number of elements in suff stat map must be equal!");
          } else{
            // ss map has the right number of elements
            for(tree::npv_it bn_it = bnv.begin(); bn_it != bnv.end(); ++bn_it){
              if(ss_test_vec[m].count( (*bn_it)->get_nid() ) != 1){
                // could not find an element in ss with key equal to the nid of a bottom node
                Rcpp::Rcout << "bottom node nid " << (*bn_it)->get_nid() << " is not a key in training suff stat map!" << std::endl;
                Rcpp::stop("there is a mismatch between bottom node id's and the keys of the training suff stat map!");
              }
            }
          }
        }
      } // closes if for checking suff_stat_map and tree are in sync.
      
      // now we need to update the value of allfit
      for(suff_stat_it ss_it = ss_train_vec[m].begin(); ss_it != ss_train_vec[m].end(); ++ss_it){
        tmp_mu = t_vec[m].get_ptr(ss_it->first)->get_mu();
        for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
          // add fit of m-th tree back to allfit and subtract it from the value of the residual
          allfit_train[*it] += tmp_mu;
          residual[*it] -= tmp_mu;
        }
      } // this loop is also O(n)
    } // closes loop over all of the trees
    // ready to update sigma
    total_sq_resid = 0.0;
    //for(int i = 0; i < n_train; i++) total_sq_resid += pow(Y_train[i] - allfit_train[i], 2.0);
    for(int i = 0; i < n_train; i++) total_sq_resid += pow(residual[i], 2.0); // sum of squared residuals
    
    scale_post = lambda * nu + total_sq_resid;
    nu_post = nu + ( (double) n_train);
    sigma = sqrt(scale_post/gen.chi_square(nu_post));
    sigma_samples(iter) = sigma;
    
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
      
      for(int i = 0; i < n_train; i++) fit_train(i, sample_index) = allfit_train[i];
      if(n_test > 0){
        for(int i = 0; i < n_train; i++) allfit_test[i] = 0.0; // reset the value of allfit_test
        for(int m = 0; m < M; m++){
          for(suff_stat_it ss_it = ss_test_vec[m].begin(); ss_it != ss_test_vec[m].end(); ++ss_it){
            tmp_mu = t_vec[m].get_ptr(ss_it->first)->get_mu(); // get the value of mu in the corresponding leaf
            for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it) allfit_test[*it] += tmp_mu;
          } // loop over the keys in the m-th sufficient stat map
        } // closes loop over trees
        //fit_ensemble(allfit_test, t_vec, di_test);
        //for(int i = 0; i < n_test; i++) fit_test(i, sample_index) = allfit_test[i];
      } // closes loop checking if we actually have test set observations.
    } // closes if that checks whether we should save anything in this iteration
  } // closes the main MCMC for loop

  Rcpp::List results;
  results["fit_train"] = fit_train;
  if(n_test > 0) results["fit_test"] = fit_test;
  results["sigma"] = sigma_samples;
  results["total_accept"] = total_accept_samples;
  if(save_trees) results["trees"] = tree_draws;
  return results;
}
