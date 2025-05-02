#include "update_tree_multi.h"

void compute_jump_posterior(std::map<int, jump_post> &jp_map, tree &t, suff_stat &ss, int &r, double &sigma, data_info &di, tree_prior_info &tree_pi)
{
  // reminder posterior of jump is N(P^-1 Theta, P^-1)
  
  // first we need to remove the fit from each leaf and as we do that, we update running estimate of P and Theta
  int i = 0;
  double z = 0.0;
  double tmp_mu = 0.0;
  jp_map.clear();
  for(suff_stat_it l_it = ss.begin(); l_it != ss.end(); ++l_it){
    tmp_mu = t.get_ptr(ss_it->first)->get_mu(); // get mu for this leaf
    jp_map.insert(std::pair<int, jump_post>(l_it->first, jump_post())); // create element in jp_map for leaf
    
    std::map<int, jump_post>::iterator jp_it = jp_map.find(l_it->first); // iterator at new element in jp_map for leaf
    jp_it->second.P = 1.0/pow(tree_pi.tau, 2.0); // prior leaf precision
    jp_it->second.Theta = tree_pi.mu0/pow(tree_pi.tau, 2.0); // contribution from prior leaf
    if(l_it->second.size() > 0){
      for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
        i = *it;
        z = *(di.z + i*di.R + r);
        di.rp[i] += z * tmp_mu; // removes fit
        jp_it->second.P += pow(z, 2.0)/pow(sigma, 2.0);
        jp_it->second.Theta += z * di.rp[i]/pow(sigma, 2.0);
      } // closes loop over observations in leaf
    } // closes if checking that there are elements in the leaf
  } // closes loop over leaf/ss_map elements
}

double compute_lil(int &nid, std::map<int, jump_post> &jp_map)
{
  if(jp_map.count(nid) == 0) Rcpp::stop("[compute_lil]: Did not find element for node in jp_map");
  std::map<int, jump_post>::iterator jp_it = jp_map.find(nid);
  return(-0.5 * log(jp_it->second.P) + 0.5 * pow(jp_it->second.Theta, 2.0)/jp_it->second.P);
}

void draw_mu(tree &t, suff_stat &ss, std::map<int, jump_post> &jp_map, data_info &di, RNG &gen)
{
  int i = 0;
  double z = 0.0;
  double post_mean = 0.0;
  double post_sd = 0.0;
  double tmp_mu = 0.0;
  for(suff_stat_it l_it = ss.begin(); l_it != ss.end(); ++l_it){
    std::map<int, jump_post>::iterator jp_it = jp_map.find(l_it->first);
    post_mean = jp_it->second.Theta/jp_it->second.P;
    post_sd = sqrt(1.0/jp_it->second.P);
    tmp_mu = gen.normal(post_mean, post_sd);
    t.get_ptr(l_it->first)->set_mu(tmp_mu);
    if(l_it->second.size() > 0){
      for(std::vector<int>::iterator it = l_it->second.begin(); it != l_it->second.end(); ++it){
        i = *it;
        z = *(di.z + i*di.R + r);
        di.rp[i] -= z * tmp_mu; // restores fit
      } // closes loop over observations in leaf
    } // closes if checking that leaf has observations
  } // loop over leafs
}

// new version that update jump posterior values and ss map for new leaf nodes in the same loop
void compute_ss_grow(suff_stat &ss, std::map<int, jump_post> &jp_map, int &nx_nid, rule_t &rule, int &r, double &sigma, data_info &di, tree_prior_info &tree_pi)
{
  jp_map.clear();
  
  int i = 0;
  double z = 0.0;
  int nxl_nid = 2*nx_nid;
  int nxr_nid = 2*nx_nid+1;
  
  // check that our data structures have element for nx but not for nxl or nxr (the proposed children)
  if(ss.count(nx_nid) == 0 || ss.count(nxl_nid) == 1 || ss.count(nxr_nid) == 1) Rcpp::stop("[compute_ss_grow_train]: something is wrong with ss_train");
  
  // if we get to here, data structure are fine
  
  // add elements to ss and jp_map for nxl
  ss.insert(std::pair<int, std::vector<int>>(nxl_nid, std::vector<int>()));
  suff_stat_it nxl_it = tmp_ss.find(nxl_nid); // points to nxl's element in ss

  jp_map.insert(std::pair<int, jump_post>(nxl_nid, jump_post()));
  std::map<int, jump_post>::iterator jpl_it = jp_map.find(nxl_nid); // points to nxl's element in jp_map
  jpl_it->second.P = 1.0/pow(tree_pi.tau, 2.0); // prior leaf precision
  jpl_it->second.Theta = tree_pi.mu0/pow(tree_pi.tau, 2.0); // contribution from prior leaf
  
  // add elements to ss and jp_map for nxr
  tmp_ss.insert(std::pair<int, std::vector<int>>(nxr_nid, std::vector<int>()));
  suff_stat_it nxr_it = tmp_ss.find(nxr_nid); // points to nxr's element in ss

  jp_map.insert(std::pair<int, jump_post>(nxr_nid, jump_post()));
  std::map<int, jump_post>::iterator jpr_it = jp_map.find(nxr_nid); // points to nxr's element in jp_map
  jpr_it->second.P = 1.0/pow(tree_pi.tau, 2.0); // prior leaf precision
  jpr_it->second.Theta = tree_pi.mu0/pow(tree_pi.tau, 2.0); // contribution from prior leaf
  
  
  for(suff_stat_it l_it = ss.begin(); l_it != ss.end(); ++l_it){
    if(l_it->first == nx_nid){
      if(l_it->second.size() > 0){
        if(!rule.is_cat){
          //double* xx_cont = 0;
          double xx_cont = *(di.x_cont + i*di.p_cont + rule.v_aa);
          for(std::vector<int>::iterator it = l_it->second.begin(); it != l_it->second.end(); ++it){
            i = *it;
            z = *(di.z + i*di.R + r);
            xx_cont = di.x_cont + i*di.p_cont;
            if(xx_cont < rule.c){
              nxl_it->second.push_back(i);
              jpl_it->second.P += pow(z, 2.0)/pow(sigma, 2.0);
              jpl_it->second.Theta += z * di.rp[i]/pow(sigma, 2.0);
            } else if(xx_cont >= rule.c){
              nxr_it->second.push_back(i);
              jpr_it->second.P += pow(z, 2.0)/pow(sigma, 2.0);
              jpr_it->second.Theta += z * di.rp[i]/pow(sigma, 2.0);
            } else{
              Rcpp::Rcout << "  i = " << i << " v = " << rule.v_aa+1 << "  value = " << xx_cont << " cutpoint = " << rule.c << std::endl;
              Rcpp::stop("[compute_ss_grow_train]: could not assign observation to left or right child in axis-aligned split!");
            }
          } // closes loop over observations in nx
        } else{
          int xx_cat = 0;
          for(std::vector<int>::iterator it = l_it->second.begin(); it != l_it->second.end(); ++it){
            i = *it;
            z = *(di.z + i*di.R + r);
            xx_cat = *(di.x_cat + i*di.p_cat + rule.v_cat);
            int l_count = rule.l_vals.count(xx_cat);
            int r_count = rule.r_vals.count(xx_cat);
            if(l_count != r_count){
              if(l_count == 1){
                nxl_it->second.push_back(i);
                jpl_it->second.P += pow(z, 2.0)/pow(sigma, 2.0);
                jpl_it->second.Theta += z * di.rp[i]/pow(sigma, 2.0);
              } else{
                nxr_it->second.push_back(i);
                jpr_it->second.P += pow(z, 2.0)/pow(sigma, 2.0);
                jpr_it->second.Theta += z * di.rp[i]/pow(sigma, 2.0);
              }
            } else{
              Rcpp::Rcout << "i = " << i << "v = " << rule.v_cat+1 << "  value = " << xx_cat << std::endl;
              Rcpp::Rcout << "left values:";
              for(set_it levels_it = rule.l_vals.begin(); levels_it != rule.l_vals.end(); ++levels_it) Rcpp::Rcout << " " << *levels_it;
              Rcpp::Rcout << std::endl;
              
              Rcpp::Rcout << "right values:";
              for(set_it levels_it = rule.r_vals.begin(); levels_it != rule.r_vals.end(); ++levels_it) Rcpp::Rcout << " " << *levels_it;
              Rcpp::Rcout << std::endl;

              Rcpp::stop("[compute_ss_grow]: could not assign observation to left or right child in categorical split!");
            } // closes if/else checking whether observation i goes to nxl or nxr
          } // closes loop over observations in nx
        } // closes if/else checking whether rule is axis-aligned or categorical
      } // closes if checking that leaf is non-empty
    } else{
      // we're not in a leaf being modified, so let's just compute posterior jump parameters
      jp_map.insert(std::pair<int, jump_post>(l_it->first, jump_post()));
      std::map<int, jump_post>::iterator jp_it = jp_map.find(l_it->first);
      jp_it->second.P = 1.0/pow(tree_pi.tau, 2.0); // prior leaf precision
      jp_it->second.Theta = tree_pi.mu0/pow(tree_pi.tau, 2.0); // contribution from prior leaf
      if(l_it->second.size() > 0){
        for(int_it it = l_it->second.begin(); it != l_it->second.end(); ++it){
          i = *it;
          z = *(di.z + i*di.R + r);
          jp_it->second.P += pow(z, 2.0)/pow(sigma, 2.0);
          jp_it->second.Theta += z * di.rp[i]/pow(sigma, 2.0);
        } // closes loop over observations in leaf
      } // closes if checking that there are elements in the leaf
    } // closes if/else checking whether we're at the node being grown
  } // closes loop over the nodes
}

// in-place modification without worrying about jump_post objects (i.e. for testing)
void compute_ss_grow(suff_stat &ss, int &nx_nid, rule_t &rule, data_info &di)
{
  int i = 0;
  int nxl_nid = 2*nx_nid;
  int nxr_nid = 2*nx_nid+1;
  
  // check that our data structures have element for nx but not for nxl or nxr (the proposed children)
  if(ss.count(nx_nid) == 0 || ss.count(nxl_nid) == 1 || ss.count(nxr_nid) == 1) Rcpp::stop("[compute_ss_grow_train]: something is wrong with ss_train");

  
  // if we get to here, data structure are fine
  suff_stat_it nx_it = ss.find(nx_nid); // points to nx's element in ss
  
  // add elements to ss and jp_map for nxl
  ss.insert(std::pair<int, std::vector<int>>(nxl_nid, std::vector<int>()));
  suff_stat_it nxl_it = ss.find(nxl_nid); // points to nxl's element in ss
  
  // add elements to ss and jp_map for nxr
  ss.insert(std::pair<int, std::vector<int>>(nxr_nid, std::vector<int>()));
  suff_stat_it nxr_it = ss.find(nxr_nid); // points to nxr's element in ss

  if(nx_it->second.size() > 0){
    if(!rule.is_cat){
      double xx_cont = 0;
      for(std::vector<int>::iterator it = nx_it->second.begin(); it != nx_it->second.end(); ++it){
        i = *it;
        xx_cont = *(di.x_cont + i*di.p_cont + rule.v_aa);
        if(xx_cont < rule.c){
          nxl_it->second.push_back(i);
        } else if(xx_cont >= rule.c){
          nxr_it->second.push_back(i);
        } else{
          Rcpp::Rcout << "  i = " << i << " v = " << rule.v_aa+1 << "  value = " << xx_cont << " cutpoint = " << rule.c << std::endl;
          Rcpp::stop("[compute_ss_grow]: could not assign observation to left or right child in axis-aligned split!");
        }
      } // closes loop over observations in nx
    } else{
      int* xx_cat = 0;
      for(std::vector<int>::iterator it = nx_it->second.begin(); it != nx_it->second.end(); ++it){
        i = *it;
        xx_cat = *(di.x_cat + i*di.p_cat + rule.v_cat);
        int l_count = rule.l_vals.count(xx_cat);
        int r_count = rule.r_vals.count(xx_cat);
        if(l_count != r_count){
          if(l_count == 1){
            nxl_it->second.push_back(i);
          } else{
            nxr_it->second.push_back(i);
          }
        } else{
          Rcpp::Rcout << "i = " << i << "v = " << rule.v_cat+1 << "  value = " << xx_cat << std::endl;
          Rcpp::Rcout << "left values:";
          for(set_it levels_it = rule.l_vals.begin(); levels_it != rule.l_vals.end(); ++levels_it) Rcpp::Rcout << " " << *levels_it;
          Rcpp::Rcout << std::endl;
          
          Rcpp::Rcout << "right values:";
          for(set_it levels_it = rule.r_vals.begin(); levels_it != rule.r_vals.end(); ++levels_it) Rcpp::Rcout << " " << *levels_it;
          Rcpp::Rcout << std::endl;
          
          Rcpp::stop("[compute_ss_grow]: could not assign observation to left or right child in categorical split!");
        } // closes if/else checking whether observation i goes to nxl or nxr
      } // closes loop over observations in nx
    } // closes if/else checking whether rule is axis-aligned or categorical
  } // closes if checking that leaf is non-empty
}

// in-place modification of ss and jp for prune moves
void compute_ss_prune(suff_stat &ss, std::map<int, jump_post> &jp_map, int &nxl_nid, int &nxr_nid, int &nx_nid, int &r, double &sigma, data_info &di, tree_prior_info &tree_pi)
{
  int i = 0;
  double z = 0.0;
  // check that our data structures have element for nxl and nxr but not for nx
  if(ss.count(nxl_nid) == 0 || ss.count(nxr_nid) == 0 || ss.count(nx_nid) == 1) Rcpp::stop("[compute_ss_prune]: something's wrong with ss_train");
  if(jp_map.count(nxl_nid) == 0 || jp_map.count(nxr_nid) == 0 || jp_map.count(nx_nid) == 1) Rcpp::stop("[compute_ss_prune]: something's wrong with jp_map");
  
  // if we get to here, data structure are fine
  suff_stat_it nxl_it = ss.find(nxl_nid); // points to nxl's element in ss
  suff_stat_it nxr_it = ss.find(nxr_nid); // points to nxr's element in ss
  
  ss.insert(std::pair<int, std::vector<int>>(nx_nid, std::vector<int>())); // add element for nx in ss
  suff_stat_it nx_it = ss.find(nx_nid); // points to nx's element in ss
  
  jp_map.insert(std::pair<int, jump_post>(nx_nid, jump_post())); // add element for nx in jp_map
  std::map<int, jump_post>::iterator jp_it = jp_map.find(nx_nid); // points to nx's element in jp_map
  jp_it->second.P = 1.0/pow(tree_pi.tau, 2.0); // prior leaf precision
  jp_it->second.Theta = tree_pi.mu0/pow(tree_pi.tau, 2.0); // contribution from prior leaf
  
  if(nxl_it->second.size() > 0){
    for(std::vector<int>::iterator it = nxl_it->second.begin(); it != nxl_it->second.end(); ++it){
      i = *it;
      z = *(di.z + i*di.R + r);
      nx_it->second.push_back(i);
      jp_it->second.P += pow(z, 2.0)/pow(sigma, 2.0);
      jp_it->second.Theta += z * di.rp[i]/pow(sigma, 2.0);
    } // closes loop over observations in nxl
  } // closes if checking that nxl is non-empty
  if(nxr_it->second.size() > 0){
    for(std::vector<int>::iterator it = nxr_it->second.begin(); it != nxr_it->second.end(); ++it){
      i = *it;
      z = *(di.z + i*di.R + r);
      nx_it->second.push_back(i);
      jp_it->second.P += pow(z, 2.0)/pow(sigma, 2.0);
      jp_it->second.Theta += z * di.rp[i]/pow(sigma, 2.0);
    } // closes loop over observations in nxr
  } // closes if checking that nxr is non-empty
}

// in-place modification of ss without worrying about jump posterior (i.e., for test data)
void compute_ss_prune(suff_stat &ss, int &nxl_nid, int &nxr_nid, int &nx_nid, data_info &di)
{

  // check that our data structures have element for nxl and nxr but not for nx
  if(ss.count(nxl_nid) == 0 || ss.count(nxr_nid) == 0 || ss.count(nx_nid) == 1) Rcpp::stop("[compute_ss_prune]: something's wrong with ss_train");
  
  // if we get to here, data structure are fine
  suff_stat_it nxl_it = ss.find(nxl_nid); // points to nxl's element in ss
  suff_stat_it nxr_it = ss.find(nxr_nid); // points to nxr's element in ss
  
  ss.insert(std::pair<int, std::vector<int>>(nx_nid, std::vector<int>())); // add element for nx in ss
  suff_stat_it nx_it = ss.find(nx_nid); // points to nx's element in ss
  
  if(nxl_it->second.size() > 0){
    for(std::vector<int>::iterator it = nxl_it->second.begin(); it != nxl_it->second.end(); ++it) nx_it->second.push_back(*it);
  } // closes if checking that nxl is non-empty
  if(nxr_it->second.size() > 0){
    for(std::vector<int>::iterator it = nxr_it->second.begin(); it != nxr_it->second.end(); ++it) nx_it->second.push_back(*it);
  } // closes if checking that nxr is non-empty
}


void grow_tree(tree &t, suff_stat &ss_train, suff_stat &ss_test, std::map<int, jump_post> &jp_map, int &accept, int &r, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
{
  std::vector<int> bn_nid_vec; // vector to hold the id's of all of the bottom nodes in the tree
  for(suff_stat_it ss_it = ss_train.begin(); ss_it != ss_train.end(); ++ss_it) bn_nid_vec.push_back(ss_it->first);
  
  int ni = floor(gen.uniform() * bn_nid_vec.size()); // randomly pick the index of the node from which we will grow
  int nx_nid = bn_nid_vec[ni]; // id of the node from which we are growing.
  tree::tree_p nx = t.get_ptr(nx_nid); // pointer to the node from which we are growing. refer to this node as nx
  tree::tree_cp nxp = nx->get_p(); // pointer to parent of nx in tree
  
  // we are ready to compute the log transition ratio:
  double q_grow_old = tree_pi.prob_b; // transition prob. of growing old tree into new tree
  double q_prune_new = 1.0 - tree_pi.prob_b; // transition prob. of pruning new true into old tree
  
  int nleaf_old = t.get_nbots(); // number of leaves in old tree
  int nnog_old = t.get_nnogs(); // number of nodes in old tree with no grandchildren (nog node)
  int nnog_new = nnog_old; // number of nodes in new tree with no grandchildren
  
  if(nxp == 0){
    // nx is the root node so transition always propose growing it
    q_grow_old = 1.0; //
    nnog_new = 1; // nx has no grandchildren in new tree
  } else if(nxp->is_nog()){
    // parent of nx has no grandchildren in old tree
    // in new tree nxp has grandchildren but nx does not
    // hence nnog_new = nnod_old
    nnog_new = nnog_old;
  } else{
    // parent of nx has grandchildren in old tree and will continue to do so in new tree
    // nx has no grandchildren in the new tree
    nnog_new = 1 + nnog_old;
  }
  
  // numerator of transition ratio: P(uniformly pick a nog node in new tree) * P(decide to prune new tree)
  // denominator of transition rate: P(uniformly pick a leaf node in old tree) * P(decide to grow old tree)
  
  double log_trans_ratio = (log(q_prune_new) - log( (double) nnog_new)) - (log(q_grow_old) - log( (double) nleaf_old));

  // for prior ratio:
  // numerator: p(grow at nx) * (1 - p(grow at nxl)) * (1 - p(grow at nxr))
  // denominator: (1 - p(grow at nx))
  // we need 1 - P(grow at nx in old tree) = 1 - alpha(1 + depth(nx))^(-beta) in denominator
  // we need P(grow at nx in new) (1 - P(grow at nxl in Tnew))(1 - P(grow at nxr in Tnew)) in numerator
  
  double p_grow_nx = tree_pi.alpha/pow(1.0 + (double) nx->get_depth(), tree_pi.beta); // prior prob of growing tree at nx
  double p_grow_nxl = tree_pi.alpha/pow(2.0 + (double) nx->get_depth(), tree_pi.beta); // prior prob of growing tree at nxl. remember depth of nxl is 1 + depth of nx
  double p_grow_nxr = tree_pi.alpha/pow(2.0 + (double) nx->get_depth(), tree_pi.beta); // prior prob of growing tree at nxr. remember depth of nxr is 1 + depth of nx
  
  double log_prior_ratio = log(p_grow_nx) + log(1.0 - p_grow_nxl) + log(1.0 - p_grow_nxr) - log(1.0 - p_grow_nx);
  
  // we now are ready to draw a decision rule
  rule_t rule;
  draw_rule(rule, t, nx_nid, di_train, tree_pi, gen); // draw the actual rule

  compute_ss_grow(ss_train, jp_map, nx_nid, rule, r, sigma, di_train, tree_pi);
  int nxl_nid = 2*nx_nid;
  int nxr_nid = 2*nx_nid+1;
  
  double nxl_lil = compute_lil(nxl_nid, jp_map); // nxl's contribution to log marginal likelihood of new tree
  double nxr_lil = compute_lil(nxr_nid, jp_map); // nxr's contribution to log marginal likelihood of new tree
  double nx_lil = compute_lil(nx_nid, jp_map); // nx's contribution to log marginal likelihood of old tree

  // likelihood ratio also needs to include some constants from prior on jumps condition on tree
  // in GROW move, the new tree has one extra leaf so there's an additional factor of tau^(-1) * exp(-mu0^2/2tau^2) from leaf prior in the numerator
  double log_like_ratio = nxl_lil + nxr_lil - nx_lil - 1.0 * log(tree_pi.tau) - 0.5 * pow(tree_pi.mu0/tree_pi.tau,2.0);
  
  double log_alpha = log_like_ratio + log_prior_ratio + log_trans_ratio; // MH ratio
  if(log_alpha > 0) log_alpha = 0.0; // if MH ratio greater than 1, we set it equal to 1. this is almost never needed
  if(gen.log_uniform() <= log_alpha){
    // accept the transition!
    
    ++(*tree_pi.rule_count); // increment running count of total number of splitting rules
    if(!rule.is_cat){
      ++(tree_pi.var_count->at(rule.v_aa)); // in our bookkeeping, continuous variables come first
    } else {
      int v_raw = rule.v_cat + di_train.p_cont;
      ++(tree_pi.var_count->at(v_raw));
    }
    
    // we are accepting the grow move, so we can eliminate the element for nx in ss_train
    ss_train.erase(nx_nid);
    jp_map.erase(nx_nid);
    
    if(di_test.n > 0){
      compute_ss_grow(ss_test, nx_nid, rule, di_test);
      // at this point, ss_test has elements for nx, nxl, and nxr
      ss_test.erase(nx_nid);
    }
    t.birth(nx_nid, rule); // actually do the birth
    accept = 1;
  } else{
    accept = 0;
    // we did not accept the move and so we need to remove element for nxl and nxr from ss_train
    ss_train.erase(nxl_nid);
    ss_train.erase(nxr_nid);
    
    jp_map.erase(nxl_nid);
    jp_map.erase(nxr_nid);
  }
}

void prune_tree(tree &t, suff_stat &ss_train, suff_stat &ss_test, std::map<int, jump_post> &jp_map, int &accept, int &r, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
{
  // first we randomly select a nog node
  tree::npv nogs_vec; // vector of pointers to nodes w/ no grandchildren
  t.get_nogs(nogs_vec);

  int ni = floor(gen.uniform() * nogs_vec.size());
  tree::tree_p nx = nogs_vec[ni]; // pointer to node whose children we will prune
  tree::tree_p nxl = nx->get_l(); // left child that will be pruned
  tree::tree_p nxr = nx->get_r(); // right child that will be pruned
  
  // transition ratio stuff
  double q_prune_old = 1.0 - tree_pi.prob_b; // transition prob that we prune old tree
  double q_grow_new = tree_pi.prob_b; // transition prob that we grow new tree
  tree::tree_p nxp = nx->get_p(); // pointer to parent node of nx in old tree
  if(nxp == 0) q_grow_new = 1.0; // nx is top node so new tree is just root and we'd always propose a GROW when encountering new tree
  else q_grow_new = tree_pi.prob_b; // nx is not top node so  given T_new, we propose grow with prob tree_pi.pb

  int nleaf_new = t.get_nbots() - 1; // new tree has one less leaf node than old tree
  int nnog_old = t.get_nnogs(); // number of nodes with no grandchildren in old tree
  
  // numerator of transition ratio: P(uniformly pick a leaf node in new tree) * P(decide to grow newtree)
  // denominator of transition ratio: P(uniformly pick a nog node in old tree) * P(decide to prune old tree)

  double log_trans_ratio = (log(q_grow_new) - log(nleaf_new)) - (log(q_prune_old) - log(nnog_old));

  // prior ratio
  // numerator: we need [1 - P(grow at nx in Tnew)] = 1 - tree_pi.alpha/pow(1 + nx->get_depth(), tree_pi.beta)
  // denom: we need [P(grow at nx in Told)] x [1 - P(grow at nxl in Told)] x [1 - P(grow at nxr in Told)]
  double p_grow_nx = tree_pi.alpha/pow(1.0 + (double) nx->get_depth(), tree_pi.beta); // prior prob of growing tree at nx
  double p_grow_nxl = tree_pi.alpha/pow(2.0 + (double) nx->get_depth(), tree_pi.beta); // prior prob of growing tree at nxl, left child of nx
  double p_grow_nxr = tree_pi.alpha/pow(2.0 + (double) nx->get_depth(), tree_pi.beta); // prior prob of growing tree nxr, right child of nx
  double log_prior_ratio = log(1.0 - p_grow_nx) - (log(1.0 - p_grow_nxl) + log(1.0 - p_grow_nxr) + log(p_grow_nx));
  
  int nx_nid = nx->get_nid(); // id for nx
  int nxl_nid = nxl->get_nid(); // id for nxl
  int nxr_nid = nxr->get_nid(); // id for nxr
  compute_ss_prune(ss_train,jp_map, nxl_nid, nxr_nid, nx_nid, r, sigma, di_train, tree_pi);
  
  double nxl_lil = compute_lil(nxl_nid, jp_map);
  double nxr_lil = compute_lil(nxr_nid, jp_map);
  double nx_lil = compute_lil(nx_nid, jp_map);
  
  // old tree has one more leaf node than new tree so there is additional factor of
  // (tau)^(-1) * exp(-mu0^2/(2 * tau^2)) in denominator of likelihood ratio that comes from the prior on leafs
  double log_like_ratio = nx_lil - nxl_lil - nxr_lil + log(tree_pi.tau) + 0.5 * pow(tree_pi.mu0/tree_pi.tau,2.0);
  
  double log_alpha = log_like_ratio + log_prior_ratio + log_trans_ratio;
  if(log_alpha > 0) log_alpha = 0; // if MH greater than we, set it equal to 1
  if(gen.log_uniform() <= log_alpha){
    // accept the proposal!
    accept = 1;
    // need to decrement several counters
    --(*tree_pi.rule_count);
    
    if(!nx->get_is_cat()){
      // we pruned away an axis-aligned or categorical rule
     --(tree_pi.var_count->at(nx->get_v_aa()));
    } else{
      int v_raw = di_train.p_cont + nx->get_v_cat();
      --(tree_pi.var_count->at(v_raw));
    }

    // need to adjust ss
    // ss_train contains element for nxl and nxr and these need to be removed
    ss_train.erase(nxl_nid);
    ss_train.erase(nxr_nid);
    jp_map.erase(nxl_nid);
    jp_map.erase(nxr_nid);
        
    if(di_test.n > 0){
      compute_ss_prune(ss_test, nxl_nid, nxr_nid, nx_nid, di_train);
      ss_test.erase(nxl_nid);
      ss_test.erase(nxr_nid);
    }
    
    t.death(nx_nid); // actually perform the death
  } else{
    accept = 0;
    // ss_train contains nx, nxl, and nxr
    // since prune was rejected, we have to get rid of nx!
    ss_train.erase(nx_nid);
    jp_map.erase(nx_nid);
  }
}

void update_tree(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, int &r, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
{
  accept = 0; // initialize indicator of MH acceptance to 0 (reject)
  double PBx = tree_pi.prob_b; // prob of proposing a birth move (typically 0.5)
  if(t.get_treesize() == 1) PBx = 1.0; // if tree is just the root, we must always GROW
  
  std::map<int, jump_post> jp_map;
  compute_jump_posterior(jp_map, ss_train, r, sigma, di_train, tree_pi);
  
  if(gen.uniform() < PBx) grow_tree_unnested(t, ss_train, ss_test, jp_map, accept, r, sigma, di_train, di_test,tree_pi, gen);
  else prune_tree(t, ss_train, ss_test, jp_map, accept, r, sigma, di_train, di_test, tree_pi, gen);
  
  draw_mu(t, jp_map, gen);

}

