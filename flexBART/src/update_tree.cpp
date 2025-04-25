#include "update_tree.h"

/*
// in-place modification of ss in grow move
void compute_ss_grow(suff_stat &ss, int &nx_nid, rule_t &rule, data_info &di)
{
  int i;
  int nxl_nid = 2*nx_nid;
  int nxr_nid = 2*nx_nid+1;
  
  // may want to add some checks that we don't already have element for left & right children and that there is something in nx_nid
  if(ss.count(nx_nid) == 1 && ss.count(nxl_nid) == 0 && ss.count(nxr_nid) == 0){
    suff_stat_it nx_it = ss.find(nx_nid); // iterator at element for nx in ss
    ss.insert(std::pair<int, std::vector<int>>(nxl_nid, std::vector<int>()));
    ss.insert(std::pair<int, std::vector<int>>(nxr_nid, std::vector<int>()));
    
    suff_stat_it nxl_it = ss.find(nxl_nid);
    suff_stat_it nxr_it = ss.find(nxr_nid);
    
    if(nx_it->second.size() > 0){
      
      if(!rule.is_cat){
        double* xx_cont = 0;
        for(std::vector<int>::iterator it = nx_it->second.begin(); it != nx_it->second.end(); ++it){
          i = *it;
          xx_cont = di.x_cont + i*di.p_cont;
          if(xx_cont[rule.v_aa] < rule.c) nxl_it->second.push_back(i);
          else if(xx_cont[rule.v_aa] >= rule.c) nxr_it->second.push_back(i);
          else{
            Rcpp::Rcout << "  i = " << i << " v = " << rule.v_aa+1 << "  value = " << xx_cont[rule.v_aa] << " cutpoint = " << rule.c << std::endl;
            Rcpp::stop("[compute_ss_grow]: could not assign observation to left or right child in axis-aligned split!");
          } // closes if/else checking whether observation i goes to nxl or nxr
        } // closes loop over observations in nx
      } else{
        int* xx_cat = 0;
        for(std::vector<int>::iterator it = nx_it->second.begin(); it != nx_it->second.end(); ++it){
          i = *it;
          xx_cat = di.x_cat + i*di.p_cat;
          int l_count = rule.l_vals.count(xx_cat[rule.v_cat]);
          int r_count = rule.r_vals.count(xx_cat[rule.v_cat]);
          if(l_count != r_count){
            if(l_count == 1) nxl_it->second.push_back(i);
            else nxr_it->second.push_back(i);
          } else{
            Rcpp::Rcout << "i = " << i << "v = " << rule.v_cat+1 << "  value = " << xx_cat[rule.v_cat] << std::endl;
            Rcpp::Rcout << "left values:";
            for(set_it levels_it = rule.l_vals.begin(); levels_it != rule.l_vals.end(); ++levels_it) Rcpp::Rcout << " " << *levels_it;
            Rcpp::Rcout << std::endl;
            
            Rcpp::Rcout << "right values:";
            for(set_it levels_it = rule.r_vals.begin(); levels_it != rule.r_vals.end(); ++levels_it) Rcpp::Rcout << " " << *levels_it;
            Rcpp::Rcout << std::endl;
            
            Rcpp::stop("[compute_ss_grow]: could not assign observation to left or right child in categorical split!");
          } // closes if/else checking whether observation i goes to nxl or nxr
        } // closes loop over observations in nx
      } // closes if/else checking whether rule is categorical or continuous
    } // closes if checking that there are observations in the leave
  } else{
    for(suff_stat_it ss_it = ss.begin(); ss_it != ss.end(); ++ss_it) Rcpp::Rcout << " " << ss_it->first;
    Rcpp::Rcout << std::endl;
    Rcpp::stop("[compute_ss_grow]: couldn't find parent nx or already had children nxl & nxr in ss");
  } // closes if/else checking that ss is valid.
}
*/


void compute_suff_stat_grow(suff_stat &orig_suff_stat, suff_stat &new_suff_stat, int &nx_nid, rule_t &rule, data_info &di)
{
  double* xx_cont;
  int* xx_cat;
  int i;
  int l_count;
  int r_count;
  
  // we are growing tree from node nx, which has id of nx_nid
  
  int nxl_nid = 2*nx_nid; // id of proposed left child of nx
  int nxr_nid = 2*nx_nid+1; // id of proposed right child of nx
  
  suff_stat_it nx_it = orig_suff_stat.find(nx_nid); // iterator at element for nx in original sufficient statistic map
  new_suff_stat.clear();
  
  // copy orig_suff_stat into new_suff_stat
  for(suff_stat_it it = orig_suff_stat.begin(); it != orig_suff_stat.end(); ++it){
    new_suff_stat.insert(std::pair<int,std::vector<int>>(it->first, it->second));
  }
  
  // now we manipulate new_suff_stat to drop nx and add nxl and nxr
  new_suff_stat.insert(std::pair<int,std::vector<int>>(nxl_nid, std::vector<int>())); // create map element for left child of nx
  new_suff_stat.insert(std::pair<int,std::vector<int>>(nxr_nid, std::vector<int>())); // create map element for right child of nx
  new_suff_stat.erase(nx_nid); // remove map element for nx as it is not a bottom leaf node in new tree
  
  suff_stat_it nxl_it = new_suff_stat.find(nxl_nid); // iterator at element for nxl in new sufficient stat map
  suff_stat_it nxr_it = new_suff_stat.find(nxr_nid); // iterator at element for nxr in new sufficient stat map
  
  // loop over all observation that were assigned to nx in original tree
  // note:
  //   nx_it->first is just the node id for nx (nx_nid)
  //   nx_it->second is a vector of integers containing the indicies of observations that land in nx
  // in helper.h we defined int_it as std::vector<int>::iterator
  if(nx_it->second.size() > 0){
    for(int_it it = nx_it->second.begin(); it != nx_it->second.end(); ++it){
      i = *it;
      if(di.x_cont != 0) xx_cont = di.x_cont + i * di.p_cont;
      if(di.x_cat != 0) xx_cat = di.x_cat + i * di.p_cat;
      
      if(!rule.is_cat){
        // axis-aligned rule
        if(xx_cont[rule.v_aa] < rule.c) nxl_it->second.push_back(i);
        else if(xx_cont[rule.v_aa] >= rule.c) nxr_it->second.push_back(i);
        else{
          Rcpp::Rcout << "  i = " << i << " v = " << rule.v_aa+1 << "  value = " << xx_cont[rule.v_aa] << " cutpoint = " << rule.c << std::endl;
          Rcpp::stop("[compute_ss_grow]: could not assign observation to left or right child in axis-aligned split!");
        }
      } else{
        // categorical rule
        // we need to see whether i-th observation's value of the categorical pred goes to left or right
        // std::set.count returns 1 if the value is in the set and 0 otherwise
        l_count = rule.l_vals.count(xx_cat[rule.v_cat]);
        r_count = rule.r_vals.count(xx_cat[rule.v_cat]);
        if(l_count == 1 && r_count == 0) nxl_it->second.push_back(i);
        else if(l_count == 0 && r_count == 1) nxr_it->second.push_back(i);
        else if(l_count == 1 && r_count == 1) Rcpp::stop("[compute_ss_grow]: observation goes to both left & right child...");
        else{
          Rcpp::Rcout << "i = " << i << "v = " << rule.v_cat+1 << "  value = " << xx_cat[rule.v_aa] << std::endl;
          Rcpp::Rcout << "left values:";
          for(set_it levels_it = rule.l_vals.begin(); levels_it != rule.l_vals.end(); ++levels_it) Rcpp::Rcout << " " << *levels_it;
          Rcpp::Rcout << std::endl;
          
          Rcpp::Rcout << "right values:";
          for(set_it levels_it = rule.r_vals.begin(); levels_it != rule.r_vals.end(); ++levels_it) Rcpp::Rcout << " " << *levels_it;
          Rcpp::Rcout << std::endl;
          
          Rcpp::stop("[compute_ss_grow]: could not assign observation to left or right child in categorical split!");
        }
      }
    } // closes loop over all entries in nx
  }
}

/*
void compute_ss_prune(suff_stat &ss, int &nxl_nid, int &nxr_nid, int &nx_nid, data_info &di)
{
  if(ss.count(nxl_nid) == 1 && ss.count(nxr_nid) == 1 && ss.count(nx_nid) == 0){
    suff_stat_it nl_it = ss.find(nxl_nid);
    suff_stat_it nr_it = ss.find(nxr_nid);
    
    ss.insert(std::pair<int, std::vector<int>>(nx_nid, std::vector<int>()));
    suff_stat_it nx_it = ss.find(nx_nid);
    
    if(nl_it->second.size() > 0){
      for(std::vector<int>::iterator it = nl_it->second.begin(); it != nl_it->second.end(); ++it) nx_it->second.push_back(*it);
    }
    if(nr_it->second.size() > 0){
      for(std::vector<int>::iterator it = nr_it->second.begin(); it != nr_it->second.end(); ++it) nx_it->second.push_back(*it);
    }
  } else{
    Rcpp::Rcout << "nxl_nid = " << nxl_nid << " nxr_nid = " << nxr_nid << " nx_nid = " << nx_nid << std::endl;
    for(suff_stat_it ss_it = ss.begin(); ss_it != ss.end(); ++ss_it) Rcpp::Rcout << " " << ss_it->first;
    Rcpp::Rcout << std::endl;
    Rcpp::stop("[compute_ss_prune]: did not find left or right node or already had parent in ss");
  } // closes if/else checking that ss is valid
}
*/

void compute_suff_stat_prune(suff_stat &orig_suff_stat, suff_stat &new_suff_stat, int &nl_nid, int &nr_nid, int &np_nid, data_info &di)
{
  //int i;
  if(orig_suff_stat.count(nl_nid) != 1) Rcpp::stop("[compute_ss_prune]: did not find left node in suff stat map");
  if(orig_suff_stat.count(nr_nid) != 1) Rcpp::stop("[compute_ss_prune]: did not find right node in suff stat map");
  
  suff_stat_it nl_it = orig_suff_stat.find(nl_nid); // iterator at element for nl in original suff stat map
  suff_stat_it nr_it = orig_suff_stat.find(nr_nid); // iterator at element for nr in original suff stat map
  
  new_suff_stat.clear();
  // this makes a completely new copy of orig_suff_stat
  for(suff_stat_it ss_it = orig_suff_stat.begin(); ss_it != orig_suff_stat.end(); ++ss_it){
    new_suff_stat.insert(std::pair<int,std::vector<int>>(ss_it->first, ss_it->second));
  }
  new_suff_stat.insert(std::pair<int,std::vector<int>>(np_nid, std::vector<int>())); // add element for np in new suff stat map
  new_suff_stat.erase(nl_nid); // delete element for nl in new suff stat map since nl has been pruned
  new_suff_stat.erase(nr_nid); // delete element for nr in new suff stat map since nr has been pruned
  
  if(new_suff_stat.count(np_nid) != 1) Rcpp::stop("[compute_ss_prune]: didn't create element in new suff stat map for np correctly");
  suff_stat_it np_it = new_suff_stat.find(np_nid); // iterator at element for np in new suff stat map
  
  // time to populate np_it
  // first let's add the elements from nl_it
  for(int_it it = nl_it->second.begin(); it != nl_it->second.end(); ++it) np_it->second.push_back( *it );
  for(int_it it = nr_it->second.begin(); it != nr_it->second.end(); ++it) np_it->second.push_back( *it );
}


double compute_lil(suff_stat &ss, int &nid, int &r, double &sigma, data_info &di, tree_prior_info &tree_pi)
{
  // reminder posterior of jump mu is N(P^-1 Theta, P^-1)
  if(ss.count(nid) != 1) Rcpp::stop("[compute_lil]: did not find node in suff stat map!");
  suff_stat_it ss_it = ss.find(nid);
  
  double P = 1.0/pow(tree_pi.tau, 2.0); // prior leaf precision
  double Theta = tree_pi.mu0/pow(tree_pi.tau, 2.0); // contribution from prior leaf mean
  
  if(ss_it->second.size() > 0){
    for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
      int i = *it;
      P += pow(di.z[r + i * di.R], 2.0)/pow(sigma, 2.0);
      Theta += di.z[r + i * di.R] * di.rp[i]/pow(sigma, 2.0);
    }
  }
  
  return(-0.5 * log(P) + 0.5 * pow(Theta,2.0) / P);
}

void draw_mu(tree &t, suff_stat &ss, int &r, double &sigma, data_info &di, tree_prior_info &tree_pi, RNG &gen)
{
  //int i;
  double P;
  double Theta;
  double post_sd;
  double post_mean;
  tree::tree_p bn; // we are modifying bn so we need a pointer not a constant pointer
  
  for(suff_stat_it ss_it = ss.begin(); ss_it != ss.end(); ++ss_it){
    bn = t.get_ptr(ss_it->first);
    if(bn == 0) Rcpp::stop("[draw_mu]: could not find node that is in suff stat map in the tree");
    else{
      P = 1.0/pow(tree_pi.tau, 2.0); // prior leaf precision
      Theta = tree_pi.mu0/pow(tree_pi.tau, 2.0); // contribution from prior leaf mean
      
      if(ss_it->second.size() > 0){
        for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
          int i = *it;
          P += pow(di.z[r + i * di.R], 2.0)/pow(sigma, 2.0);
          Theta += di.z[r + i * di.R] * di.rp[i]/pow(sigma, 2.0);
        }
      }
      post_sd = sqrt(1.0/P);
      post_mean = Theta/P;
      bn->set_mu(gen.normal(post_mean, post_sd));
    }
  }
}



void grow_tree_unnested(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, int &r, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
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
  draw_rule_unnested(rule, t, nx_nid, di_train, tree_pi, gen); // draw the actual rule

  // at this point we have the proposed rule and are ready to update our sufficient statistic map
  suff_stat prop_ss_train;
  compute_suff_stat_grow(ss_train, prop_ss_train, nx_nid, rule, di_train); // figure out which training observations from nx move to nxl and nxr
  //compute_ss_grow(ss_train, nx_nid, rule, di_train);
  // at this point, ss_train contains an element for nx, nxl, and nxr
  
  int nxl_nid = 2*nx_nid; // id for the left child of nx
  int nxr_nid = 2*nx_nid+1; // id for right child of nx
  
  double nxl_lil = compute_lil(prop_ss_train, nxl_nid, r, sigma, di_train, tree_pi); // nxl's contribution to log marginal likelihood of new tree
  double nxr_lil = compute_lil(prop_ss_train, nxr_nid, r, sigma, di_train, tree_pi); // nxr's contribution to log marginal likelihood of new tree
  double nx_lil = compute_lil(ss_train, nx_nid, r, sigma, di_train, tree_pi); // nx's contribution to log marginal likelihood of old tree
  
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
    /*
    // we are accepting the grow move, so we can eliminate the element for nx in ss_train
    ss_train.erase(nx_nid);
    if(di_test.n > 0){
      compute_ss_grow(ss_test, nx_nid, rule, di_test);
      // at this point, ss_test has elements for nx, nxl, and nxr
      ss_test.erase(nx_nid);
    }
   */

    
    suff_stat_it nxl_it = prop_ss_train.find(nxl_nid); // iterator at element for nxl in the proposed suff_stat map
    suff_stat_it nxr_it = prop_ss_train.find(nxr_nid); // iterator at element for nxr in the proposed suff_stat map
    
    if(nxl_it == prop_ss_train.end() || nxr_it == prop_ss_train.end()){
      // couldn't find a key in prop_ss_train equal to nxl_nid or nxr_nid
      Rcpp::Rcout << "[grow_tree]: sufficient stat map for training data not updated correctly in grow move!" << std::endl;
      Rcpp::Rcout << "  left child id = " << nxl_nid << "  right child = " << nxr_nid << std::endl;
      Rcpp::Rcout << "  available ids in map:";
      for(suff_stat_it it = prop_ss_train.begin(); it != prop_ss_train.end(); ++it) Rcpp::Rcout << " " << it->first;
      Rcpp::Rcout << std::endl;
      Rcpp::stop("missing id for either left or right child in proposed suff_stat_map!");
    }
    
    ss_train.insert(std::pair<int,std::vector<int>>(nxl_nid, nxl_it->second)); // add element for nxl in sufficient stat map
    ss_train.insert(std::pair<int,std::vector<int>>(nxr_nid, nxr_it->second)); // add element for nxr in sufficient stat map
    ss_train.erase(nx_nid); // remove element for nx in sufficient stat map
    if(di_test.n > 0){
      suff_stat prop_ss_test;
      compute_suff_stat_grow(ss_test, prop_ss_test, nx_nid, rule, di_test); // figure out which testing observation from nx more to nxl and nxr
      nxl_it = prop_ss_test.find(nxl_nid);
      nxr_it = prop_ss_test.find(nxr_nid);
      if(nxl_it == prop_ss_test.end() || nxr_it == prop_ss_test.end()){
        // couldn't find a key in prop_ss_train equal to nxl_nid or nxr_nid
        Rcpp::Rcout << "[grow_tree]: sufficient stat map for testing data not updated correctly in grow move!" << std::endl;
        Rcpp::Rcout << "  left child id = " << nxl_nid << "  right child = " << nxr_nid << std::endl;
        Rcpp::Rcout << "  available ids in map:";
        for(suff_stat_it it = prop_ss_test.begin(); it != prop_ss_test.end(); ++it) Rcpp::Rcout << " " << it->first;
        Rcpp::Rcout << std::endl;
        Rcpp::stop("missing id for either left or right child in proposed suff_stat_map!");
      }
      ss_test.insert(std::pair<int,std::vector<int>>(nxl_nid, nxl_it->second));
      ss_test.insert(std::pair<int,std::vector<int>>(nxr_nid, nxr_it->second));
      ss_test.erase(nx_nid);
    }
    
    t.birth(nx_nid, rule); // actually do the birth
    accept = 1;
  } else{
    accept = 0;
    // we did not accept the move and so we need to remove element for nxl and nxr from ss_train
    //ss_train.erase(nxl_nid);
    //ss_train.erase(nxr_nid);
  }
}


void grow_tree_nested(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, int &r, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
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
  draw_rule_nested(rule, t, nx_nid, di_train, tree_pi, gen); // draw the actual rule

  // at this point we have the proposed rule and are ready to update our sufficient statistic map
  suff_stat prop_ss_train;
  compute_suff_stat_grow(ss_train, prop_ss_train, nx_nid, rule, di_train); // figure out which training observations from nx move to nxl and nxr
  //compute_ss_grow(ss_train, nx_nid, rule, di_train);
  // at this point, ss_train contains an element for nx, nxl, and nxr
  
  int nxl_nid = 2*nx_nid; // id for the left child of nx
  int nxr_nid = 2*nx_nid+1; // id for right child of nx
  
  double nxl_lil = compute_lil(prop_ss_train, nxl_nid, r, sigma, di_train, tree_pi); // nxl's contribution to log marginal likelihood of new tree
  double nxr_lil = compute_lil(prop_ss_train, nxr_nid, r, sigma, di_train, tree_pi); // nxr's contribution to log marginal likelihood of new tree
  double nx_lil = compute_lil(ss_train, nx_nid, r, sigma, di_train, tree_pi); // nx's contribution to log marginal likelihood of old tree
  
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
   
    // we need to update ss, the sufficient statistic object
    /*
    // we are accepting the grow move, so we can eliminate the element for nx in ss_train
    ss_train.erase(nx_nid);
    
    if(di_test.n > 0){
      compute_ss_grow(ss_test, nx_nid, rule, di_test);
      // at this point, ss_test has elements for nx, nxl, and nxr
      ss_test.erase(nx_nid);
    }
     */

    suff_stat_it nxl_it = prop_ss_train.find(nxl_nid); // iterator at element for nxl in the proposed suff_stat map
    suff_stat_it nxr_it = prop_ss_train.find(nxr_nid); // iterator at element for nxr in the proposed suff_stat map
    
    if(nxl_it == prop_ss_train.end() || nxr_it == prop_ss_train.end()){
      // couldn't find a key in prop_ss_train equal to nxl_nid or nxr_nid
      Rcpp::Rcout << "[grow_tree]: sufficient stat map for training data not updated correctly in grow move!" << std::endl;
      Rcpp::Rcout << "  left child id = " << nxl_nid << "  right child = " << nxr_nid << std::endl;
      Rcpp::Rcout << "  available ids in map:";
      for(suff_stat_it it = prop_ss_train.begin(); it != prop_ss_train.end(); ++it) Rcpp::Rcout << " " << it->first;
      Rcpp::Rcout << std::endl;
      Rcpp::stop("missing id for either left or right child in proposed suff_stat_map!");
    }
    
    ss_train.insert(std::pair<int,std::vector<int>>(nxl_nid, nxl_it->second)); // add element for nxl in sufficient stat map
    ss_train.insert(std::pair<int,std::vector<int>>(nxr_nid, nxr_it->second)); // add element for nxr in sufficient stat map
    ss_train.erase(nx_nid); // remove element for nx in sufficient stat map
    if(di_test.n > 0){
      suff_stat prop_ss_test;
      compute_suff_stat_grow(ss_test, prop_ss_test, nx_nid, rule, di_test); // figure out which testing observation from nx more to nxl and nxr
      nxl_it = prop_ss_test.find(nxl_nid);
      nxr_it = prop_ss_test.find(nxr_nid);
      if(nxl_it == prop_ss_test.end() || nxr_it == prop_ss_test.end()){
        // couldn't find a key in prop_ss_train equal to nxl_nid or nxr_nid
        Rcpp::Rcout << "[grow_tree]: sufficient stat map for testing data not updated correctly in grow move!" << std::endl;
        Rcpp::Rcout << "  left child id = " << nxl_nid << "  right child = " << nxr_nid << std::endl;
        Rcpp::Rcout << "  available ids in map:";
        for(suff_stat_it it = prop_ss_test.begin(); it != prop_ss_test.end(); ++it) Rcpp::Rcout << " " << it->first;
        Rcpp::Rcout << std::endl;
        Rcpp::stop("missing id for either left or right child in proposed suff_stat_map!");
      }
      ss_test.insert(std::pair<int,std::vector<int>>(nxl_nid, nxl_it->second));
      ss_test.insert(std::pair<int,std::vector<int>>(nxr_nid, nxr_it->second));
      ss_test.erase(nx_nid);
    }
    t.birth(nx_nid, rule); // actually do the birth
    accept = 1;
  } else{
    accept = 0;
    //ss_train.erase(nxl_nid);
    //ss_train.erase(nxr_nid);
  }
}


void prune_tree(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, int &r, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
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
  
  // likelihood ratio
  suff_stat prop_ss_train;
  int nx_nid = nx->get_nid(); // id for nx
  int nxl_nid = nxl->get_nid(); // id for nxl
  int nxr_nid = nxr->get_nid(); // id for nxr
  
  compute_suff_stat_prune(ss_train, prop_ss_train, nxl_nid, nxr_nid, nx_nid, di_train); // create a sufficient statistic map for the new tree
  //compute_ss_prune(ss_train, nxl_nid, nxr_nid, nx_nid, di_train);
  
  double nxl_lil = compute_lil(ss_train, nxl_nid, r, sigma, di_train, tree_pi);
  double nxr_lil = compute_lil(ss_train, nxr_nid, r, sigma, di_train, tree_pi);
  double nx_lil = compute_lil(prop_ss_train, nx_nid, r, sigma, di_train, tree_pi);
  
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
/*
    // need to adjust ss
    // ss_train contains element for nxl and nxr and these need to be removed
    ss_train.erase(nxl_nid);
    ss_train.erase(nxr_nid);
        
    if(di_test.n > 0){
      compute_ss_prune(ss_test, nxl_nid, nxr_nid, nx_nid, di_train);
      ss_test.erase(nxl_nid);
      ss_test.erase(nxr_nid);
    }
*/
    
    suff_stat_it nx_it = prop_ss_train.find(nx_nid); // iterator at element for nx in suff stat map for new tree
    if(nx_it == prop_ss_train.end()){
      // did not find nx_nid in the keys of prop_ss_train
      Rcpp::Rcout << "[prune_tree]: did not find id of new leaf node in the keys of training sufficient stat map" << std::endl;
      Rcpp::Rcout << "  id of new leaf: " << nx_nid << std::endl;
      Rcpp::Rcout << "  ids in map:";
      for(suff_stat_it it = prop_ss_train.begin(); it != prop_ss_train.end(); ++it) Rcpp::Rcout << " " << it->first;
      Rcpp::Rcout << std::endl;
      Rcpp::stop("missing id for new leaf node in prune move in training sufficient stat map");
    } else{
      ss_train.erase(nxl_nid); // delete entry for nxl in suff stat map
      ss_train.erase(nxr_nid); // delete entry for nxr in suff stat map
      ss_train.insert(std::pair<int, std::vector<int>>(nx_nid, nx_it->second)); // add an entry for nx in suff stat map
    }
    
    if(di_test.n > 0){
      suff_stat prop_ss_test;
      compute_suff_stat_prune(ss_test, prop_ss_test, nxl_nid, nxr_nid, nx_nid, di_test);
      nx_it = prop_ss_test.find(nx_nid);
      if(nx_it == prop_ss_test.end()){
        // did not find nx_nid in the keys of prop_ss_test
        Rcpp::Rcout << "[prune_tree]: did not find id of new leaf node in the keys of testing sufficient stat map" << std::endl;
        Rcpp::Rcout << "  id of new leaf: " << nx_nid << std::endl;
        Rcpp::Rcout << "  ids in map:";
        for(suff_stat_it it = prop_ss_test.begin(); it != prop_ss_test.end(); ++it) Rcpp::Rcout << " " << it->first;
        Rcpp::Rcout << std::endl;
        Rcpp::stop("missing id for new leaf node in prune move in training sufficient stat map");
      } else{
        ss_test.erase(nxl_nid); // delete entry for nxl in suff stat map
        ss_test.erase(nxr_nid); // delete entry for nxr in suff stat map
        ss_test.insert(std::pair<int, std::vector<int>>(nx_nid, nx_it->second)); // add an entry for nx in suff stat map
      }
    }
    t.death(nx_nid); // actually perform the death
  } else{
    accept = 0;
    // ss_train contains nx, nxl, and nxr
    // since prune was rejected, we have to get rid of nx!
    //ss_train.erase(nx_nid);
  }
}

void update_tree_unnested(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, int &r, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
{
  accept = 0; // initialize indicator of MH acceptance to 0 (reject)
  double PBx = tree_pi.prob_b; // prob of proposing a birth move (typically 0.5)
  if(t.get_treesize() == 1) PBx = 1.0; // if tree is just the root, we must always GROW
  
  if(gen.uniform() < PBx) grow_tree_unnested(t, ss_train, ss_test, accept, r, sigma, di_train, di_test,tree_pi, gen);
  else prune_tree(t, ss_train, ss_test, accept, r, sigma, di_train, di_test, tree_pi, gen);

  // by this point, the decision tree has been updated so we can draw new jumps.
  draw_mu(t, ss_train, r, sigma,di_train, tree_pi, gen);
}


void update_tree_nested(tree &t, suff_stat &ss_train, suff_stat &ss_test, int &accept, int &r, double &sigma, data_info &di_train, data_info &di_test, tree_prior_info &tree_pi, RNG &gen)
{
  accept = 0; // initialize indicator of MH acceptance to 0 (reject)
  double PBx = tree_pi.prob_b; // prob of proposing a birth move (typically 0.5)
  if(t.get_treesize() == 1) PBx = 1.0; // if tree is just the root, we must always GROW
  
  if(gen.uniform() < PBx) grow_tree_nested(t, ss_train, ss_test, accept, r, sigma, di_train, di_test,tree_pi, gen);
  else prune_tree(t, ss_train, ss_test, accept, r, sigma, di_train, di_test, tree_pi, gen);

  // by this point, the decision tree has been updated so we can draw new jumps.
  draw_mu(t, ss_train, r, sigma,di_train, tree_pi, gen);
}
