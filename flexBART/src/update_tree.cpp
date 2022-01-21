#include "update_tree.h"

// when grow_tree is called: ss maps observations to the leaf nodes of the current tree
// if the grow move is accepted, ss gets updated
void grow_tree(tree &t, suff_stat &ss, int &accept, double &sigma, data_info &di, tree_prior_info &tree_pi, RNG &gen)
{
  
  std::vector<int> bn_nid_vec; // vector to hold the id's of all of the bottom nodes in the tree
  for(suff_stat_it ss_it = ss.begin(); ss_it != ss.end(); ++ss_it) bn_nid_vec.push_back(ss_it->first);
  
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
  rule_t rule; // rule is currently empty; we will populate it now.

  int rule_counter = 0; // we are allowed multiple tries to draw a valid random combination or categorical rule
  double c_upper = 1.0; // upper bound for range of cutpoints in axis aligned split
  double c_lower = -1.0; // lower bound for range of cutpoints in axis aligned split
  double tmp_weight = 0.0; // weights of random combination
  double c_max = 1.0; // upper bound for absolute value of cutpoint in random combination split
  
  double unif = gen.uniform();
  if(unif < tree_pi.prob_aa){
    // axis aligned split
    rule.is_cat = false;
    rule.is_rc = false;
    rule.v_aa = gen.multinomial(di.p_cont, tree_pi.theta_aa);
    nx->get_rg_aa(rule.v_aa, c_lower, c_upper); // what is the range of valid cutpoints for the selected variable (based only on tree topology)
    if(c_lower >= c_upper){
      // this should really throw an error... but there's no harm in proposing a trivial split.
      c_lower = -1.0;
      c_upper = 1.0;
    }
    rule.c = gen.uniform(c_lower, c_upper); // uniformly select the cutpoint.
  } else if(unif < tree_pi.prob_aa + tree_pi.prob_rc){
    // random combination split
    rule.is_cat = false;
    rule.is_rc = true;

    while( (rule.rc_weight.size() < 2) && (rule_counter < 1000) ){
      rule.rc_weight.clear();
      c_max = 0.0;
      for(int j = 0; j < di.p_cont; j++){
        if(gen.uniform() < (*tree_pi.theta_rc)){
          tmp_weight = gen.uniform(-1.0,1.0); // Breiman used Uniform(-1,1) weights and so shall we
          rule.rc_weight.insert(std::pair<int,double>(j,tmp_weight));
          c_max += fabs(tmp_weight);
        }
      }
      ++(rule_counter);
    }
    if(rule.rc_weight.size() < 2) Rcpp::stop("[propose_rule]: failed to generate a valid random combination rule in 1000 attempts!");
    else{
      rule.c = gen.uniform(-1.0,1.0) * c_max;
    }
  } else{
    // categorical split
    rule.is_cat = true;
    rule.is_rc = false;
    
    rule.v_cat = gen.multinomial(di.p_cat, tree_pi.theta_cat); // pick the categorical variable on which to split
    std::set<int> avail_levels = di.cat_levels->at(rule.v_cat); // get the full set of levels for this variable
    nx->get_rg_cat(rule.v_cat, avail_levels); // determine the set of levels available at nx.

    // if there is only one level left for this variable at nx, we will just propose a trivial split
    // and will reset the value of avail_levels to be the full set of all levels for the variable
    if(avail_levels.size() <= 1) avail_levels = di.cat_levels->at(rule.v_cat);
    
    rule.l_vals.clear();
    rule.r_vals.clear();
    
    if(tree_pi.mst_split){
      // do a split based on pruning a random MST of the adjacency graph
      graph_partition(avail_levels, rule.l_vals, rule.r_vals, di.adj_support->at(rule.v_cat), di.K->at(rule.v_cat), tree_pi.mst_reweight, gen);
    } else{
      // we can split the levels independently w/ prob 0.5 to go to each child
      rule_counter = 0;
      while( ((rule.l_vals.size() == 0) || (rule.r_vals.size() == 0)) && rule_counter < 1000 ){
        rule.l_vals.clear();
        rule.r_vals.clear();
        for(set_it it = avail_levels.begin(); it != avail_levels.end(); ++it){
          if(gen.uniform() <= 0.5) rule.l_vals.insert(*it);
          else rule.r_vals.insert(*it);
        }
        ++(rule_counter);
      }
      if(rule_counter == 1000){
        Rcpp::stop("failed to generate valid categorical split in 1000 attempts"); // this should almost surely not get triggered.
      }
    }
    if( (rule.l_vals.size() == 0) || (rule.r_vals.size() == 0) ) Rcpp::stop("proposed an invalid categorical rule!");
  }
    
  // at this point we have the proposed rule and are ready to update our sufficient statistic map
  suff_stat prop_ss;
  compute_suff_stat_grow(ss, prop_ss, nx_nid, rule, t, di); // figure out which observations from nx move to nxl and nxr
  int nxl_nid = 2*nx_nid; // id for the left child of nx
  int nxr_nid = 2*nx_nid+1; // id for right child of nx
  
  double nxl_lil = compute_lil(prop_ss, nxl_nid, sigma, di, tree_pi); // nxl's contribution to log marginal likelihood of new tree
  double nxr_lil = compute_lil(prop_ss, nxr_nid, sigma, di, tree_pi); // nxr's contribution to log marginal likelihood of new tree
  double nx_lil = compute_lil(ss, nx_nid, sigma, di, tree_pi); // nx's contribution to log marginal likelihood of old tree
  
  // likelihood ratio also needs to include some constants from prior on jumps condition on tree
  // in GROW move, the new tree has one extra leaf so there's an additional factor of tau^(-1) * exp(-mu0^2/2tau^2) from leaf prior in the numerator
  double log_like_ratio = nxl_lil + nxr_lil - nx_lil - 1.0 * log(tree_pi.tau) - 0.5 * pow(tree_pi.mu0/tree_pi.tau,2.0);
  
  double log_alpha = log_like_ratio + log_prior_ratio + log_trans_ratio; // MH ratio
  if(log_alpha > 0) log_alpha = 0.0; // if MH ratio greater than 1, we set it equal to 1. this is almost never needed
  if(gen.log_uniform() <= log_alpha){
    // accept the transition!
    if(!rule.is_cat && !rule.is_rc){
      // axis aligned rule
      ++(*tree_pi.aa_rule_count); // increment running count of total number of axis aligned splitting rules
      ++(tree_pi.aa_var_count->at(rule.v_aa)); // increment running count of number of times we've split on the v-th continuous predictor
    } else if(!rule.is_cat && rule.is_rc){
      // random combination rule
      ++(*tree_pi.rc_rule_count); // increase running count of the total number of random combination splitting rules
      for(rc_it it = rule.rc_weight.begin(); it != rule.rc_weight.end(); ++it){
        ++(*tree_pi.rc_var_count); // update the *total* number of variables that are involved in a random combination rule
      }
    } else if(rule.is_cat && !rule.is_rc){
      //categorical rule
      ++(*tree_pi.cat_rule_count); // increment running count of the total number of categorical splitting rules
      ++(tree_pi.cat_var_count->at(rule.v_cat)); // increment running count of number of times we've split on the v-th categorical predictor
    } else{
      Rcpp::Rcout << "[grow_tree]: accepted a birth at node " << nx_nid << " but unable to figure out the rule type" << std::endl;
      t.print();
      Rcpp::stop("[grow tree]: cannot resolve rule type");
    }
    
    // we need to update ss, the sufficient statistic object
    // this accounting is checked in test_grow_tree();

    suff_stat_it nxl_it = prop_ss.find(nxl_nid); // iterator at element for nxl in the proposed suff_stat map
    suff_stat_it nxr_it = prop_ss.find(nxr_nid); // iterator at element for nxr in the proposed suff_stat map
    ss.insert(std::pair<int,std::vector<int>>(nxl_nid, nxl_it->second)); // add element for nxl in sufficient stat map
    ss.insert(std::pair<int,std::vector<int>>(nxr_nid, nxr_it->second)); // add element for nxr in sufficient stat map
    ss.erase(nx_nid); // remove element for nx in sufficient stat map
    t.birth(nx_nid, rule); // actually do the birth
    accept = 1;
  } else{
    accept = 0;
    // don't do anything with rule counters or variable splitting counters etc.
  }
}


void prune_tree(tree &t, suff_stat &ss, int &accept, double &sigma, data_info &di, tree_prior_info &tree_pi, RNG &gen)
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
  suff_stat prop_ss;
  int nx_nid = nx->get_nid(); // id for nx
  int nxl_nid = nxl->get_nid(); // id for nxl
  int nxr_nid = nxr->get_nid(); // id for nxr
  
  compute_suff_stat_prune(ss, prop_ss, nxl_nid, nxr_nid, nx_nid, t, di); // create a sufficient statistic map for the new tree
  double nxl_lil = compute_lil(ss, nxl_nid, sigma, di, tree_pi);
  double nxr_lil = compute_lil(ss, nxr_nid, sigma, di, tree_pi);
  double nx_lil = compute_lil(prop_ss, nx_nid, sigma, di, tree_pi);
  
  // old tree has one more leaf node than new tree so there is additional factor of
  // (tau)^(-1) * exp(-mu0^2/(2 * tau^2)) in denominator of likelihood ratio that comes from the prior on leafs
  double log_like_ratio = nx_lil - nxl_lil - nxr_lil + log(tree_pi.tau) + 0.5 * pow(tree_pi.mu0/tree_pi.tau,2.0);
  
  double log_alpha = log_like_ratio + log_prior_ratio + log_trans_ratio;
  if(log_alpha > 0) log_alpha = 0; // if MH greater than we, set it equal to 1
  if(gen.log_uniform() <= log_alpha){
    // accept the proposal!
    accept = 1;
    // need to decrement several counters
    
    if(!nx->get_is_cat() && !nx->get_is_rc()){
      // nx was an axis-aligned rules
      --(*tree_pi.aa_rule_count);
      --(tree_pi.aa_var_count->at(nx->get_v_aa()));
    } else if(!nx->get_is_cat() && nx->get_is_rc()){
      // nx was random combination rule
      std::map<int,double> rc_weight = nx->get_rc_weight();
      --(*tree_pi.rc_rule_count);
      for(rc_it it = rc_weight.begin(); it != rc_weight.end(); ++it){
        --(*tree_pi.rc_var_count);
      }
    } else if(nx->get_is_cat() && !nx->get_is_rc()){
      // nx was a categorical rule
      --(*tree_pi.cat_rule_count);
      --(tree_pi.cat_var_count->at(nx->get_v_cat()));
    } else{
      Rcpp::Rcout << "[prune_tree]: accepted a prune at nog node " << nx_nid << " but unable to figure out its rule type" << std::endl;
      t.print();
      Rcpp::stop("[prune tree]: cannot resolve rule type");
    }

    // need to adjust ss
    // this accounting is checked in test_prune_tree();
    suff_stat_it nx_it = prop_ss.find(nx_nid); // iterator at element for nx in suff stat map for new tree
    ss.erase(nxl_nid); // delete entry for nxl in suff stat map
    ss.erase(nxr_nid); // delete entry for nxr in suff stat map
    ss.insert(std::pair<int, std::vector<int>>(nx_nid, nx_it->second)); // add an entry for nx in suff stat map
    t.death(nx_nid); // actually perform the death
  } else{
    accept = 0;
  }
  
}

void update_tree(tree &t, suff_stat &ss, int &accept, double &sigma, data_info &di, tree_prior_info &tree_pi, RNG &gen)
{
  accept = 0; // initialize indicator of MH acceptance to 0 (reject)
  double PBx = tree_pi.prob_b; // prob of proposing a birth move (typically 0.5)
  if(t.get_treesize() == 1) PBx = 1.0; // if tree is just the root, we must always GROW
  
  if(gen.uniform() < PBx) grow_tree(t, ss, accept, sigma, di,tree_pi, gen);
  else prune_tree(t, ss, accept, sigma, di, tree_pi, gen);

  // by this point, the decision tree has been updated so we can draw new jumps.
  draw_mu(t, ss, sigma, di, tree_pi, gen);
}
