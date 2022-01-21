#include "draw_tree.h"



void draw_tree(tree &t, data_info &di, tree_prior_info &tree_pi, bool mst_splits, RNG &gen)
{
  if(t.get_treesize() > 1){
    t.to_null(); // prune tree back to root
  }
  tree::npv bnv;
  int dnx; // depth of node nx
  int max_depth = 0; // depth of deepest leaf
  int prev_max_depth = 0; // max depth from previous iteration
  
  double PGnx = 0.0; // probability of growing at node nx
  bool grow = true;
  int counter = 0;
  
  // stuff for decision rules
  rule_t rule;
  int rule_counter;
  double c_upper = 1.0; // for axis-aligned rules
  double c_lower = -1.0; // for axis-aligned rules
  double tmp_weight = 0.0;
  double c_max = 1.0; //
  
  
  // when we cut the edge from the MST, we can either pick an edge uniformly (reweight = false)
  // or we can delete an edge with prob. proportional to the size of the smallest cluster that results (reweight = true)
  bool reweight = true; // if false we produce lots of singletons
  
  // stop growing after 100 attempts of growing the bottom nodes (this is a gross upper bound)
  while(grow && counter < 100){
    prev_max_depth = max_depth;
    bnv.clear();
    t.get_bots(bnv); // get the bottom nodes of the tree (could be done more efficiently but not super important)
    
    Rcpp::Rcout << "Starting round " << counter;
    Rcpp::Rcout << "  tree size = " << t.get_treesize();
    
    
    for(tree::npv_it l_it = bnv.begin(); l_it != bnv.end(); ++l_it){
      dnx = (*l_it)->get_depth(); // remember l_it is a pointer to an element in bnv, which is itself a pointer, hence the need for (*)->
      if(dnx > max_depth) max_depth = dnx; // the node we're at is deeper than the maximum depth of the tree in the previous iteration
    }
    
    Rcpp::Rcout << "  max depth = " << max_depth << std::endl;
    
    if( (max_depth < prev_max_depth) || (max_depth > 1 + prev_max_depth) ){
      // each time through the loop we can only grow the deepest leaf nodes
      // we should *never* encounter this condition but it's here to be safe
      Rcpp::Rcout << "max_depth = " << max_depth << " prev_max_depth = " << prev_max_depth << std::endl;
      Rcpp::stop("[draw_tree]: max depth should be prev_max_depth or prev_max_depth + 1!");
    } else if(max_depth == prev_max_depth && max_depth != 0){
      // tree didn't grow in the last iteration so we should break out of the loop
      break;
    } else {
      grow = false;
      //Rcpp::Rcout << "max_depth = " << max_depth << " prev_max_depth = " << prev_max_depth << std::endl;
      for(tree::npv_it l_it = bnv.begin(); l_it != bnv.end(); ++l_it){
        dnx = (*l_it)->get_depth();
        //Rcpp::Rcout << "trying node " << (*l_it)->get_nid() << "at depth " << dnx << std::endl;
        if(dnx == max_depth){
          // current node nx is at the maximum depth, we will try to grow the tree from nx
          PGnx = tree_pi.alpha/pow(1.0 + (double) dnx, tree_pi.beta);
          double tmp_unif = gen.uniform();
          
          //Rcpp::Rcout << " PGnx = " << PGnx << "tmp_unif = " << tmp_unif << std::endl;
          Rcpp::Rcout << "  node " << (*l_it)->get_nid() << " PGnx = " << PGnx << " tmp_unif = " << tmp_unif;

          if(tmp_unif < PGnx){
            Rcpp::Rcout << " can grow...";
            grow = true;
            // we're actually going to grow the tree!
            rule.clear(); // clear out the rule
            double unif = gen.uniform();
            if(unif < tree_pi.prob_aa){
              // axis-aligned split
              c_lower = -1.0;
              c_upper = 1.0;
              
              rule.is_cat = false;
              rule.is_rc = false;
              rule.v_aa = gen.multinomial(di.p_cont, tree_pi.theta_aa); // theta_aa is just a pointer
              (*l_it)->get_rg_aa(rule.v_aa, c_lower, c_upper);
              if(c_lower >= c_upper){
                c_lower = -1.0;
                c_upper = 1.0;
              }
              rule.c = gen.uniform(c_lower, c_upper);
              Rcpp::Rcout << " aa rule...";
            } else if(unif < tree_pi.prob_aa + tree_pi.prob_rc){
              // random combination split
              rule.is_cat = false;
              rule.is_rc = true;
              
              rule_counter = 0;
              while( (rule.rc_weight.size() < 2) && (rule_counter < 1000) ){
                rule.rc_weight.clear();
                c_max = 0.0;
                for(int j = 0; j < di.p_cont; j++){
                  if(gen.uniform() < (*tree_pi.theta_rc) ){
                    tmp_weight = gen.uniform(-1.0, 1.0); // Breiman used Uniform(-1,1) weights and so shall we
                    rule.rc_weight.insert(std::pair<int,double>(j, tmp_weight));
                    c_max += fabs(tmp_weight);
                  }
                }
                ++(rule_counter);
              }
              if(rule.rc_weight.size() < 2) Rcpp::stop("failed to generate a valid random combination rule in 1000 attempts!");
              else{
                rule.c = gen.uniform(-1.0, 1.0) * c_max;
              }
              Rcpp::Rcout << "  rc rule...";
            } else{
              // categorical split
              rule.is_cat = true;
              rule.is_rc = false;
              
              rule.v_cat = gen.multinomial(di.p_cat, tree_pi.theta_cat);
              std::set<int> avail_levels = di.cat_levels->at(rule.v_cat);
              (*l_it)->get_rg_cat(rule.v_cat, avail_levels); // get the available levels at the current node
              Rcpp::Rcout << " # avail levels " << avail_levels.size();
              if(avail_levels.size() >= 2){
                // this is only here in draw trees
                if(mst_splits){
                  graph_partition(avail_levels, rule.l_vals, rule.r_vals, di.adj_support->at(rule.v_cat), di.K->at(rule.v_cat), reweight, gen);
                }
                else{
                  rule_counter = 0;
                  while( ((rule.l_vals.size() == 0) || (rule.r_vals.size() == 0) ) && rule_counter < 1000 ){
                    rule.l_vals.clear();
                    rule.r_vals.clear();
                    for(set_it it = avail_levels.begin(); it != avail_levels.end(); ++it){
                      if(gen.uniform() <= 0.5) rule.l_vals.insert(*it);
                      else rule.r_vals.insert(*it);
                    }
                    ++(rule_counter);
                  }
                }
              } else{
                Rcpp::stop("not enough levels found!");
              }
            } // closes if/else's determining what type of rule to propose
            t.birth( (*l_it)->get_nid(), rule);
          } // closes if checking that we're actually trying to grow the tree
        } else{
          Rcpp::Rcout << "  node " << (*l_it)->get_nid() << " not at max depth. moving on";
        }// closes if/else checking that we're at node at the deepest level of the tree
        Rcpp::Rcout << std::endl;
      } // closes for loop over all bottom nodes in the tree
    } // closes if/else checking that max depth is valid
    ++(counter);
  } // closes main while loop
}

/*
void draw_tree(tree &t, data_info &di, tree_prior_info &tree_pi, RNG &gen){
  
  // di.cat_levels is a pointer to a vector of unordered sets that hold
  // the unique values of each categorical variable
  
  
  
  tree::npv bnv; // all bottom nodes
  size_t dnx; // depth of the node nx
  size_t max_depth = 0; // depth of deepest leaf
  size_t prev_max_depth = 0; // max depth from the previous iteration
  
  double PGnx = 0.0; // probability of growing at node nx
  bool grow = true; // can we continue growing?
  int counter = 0; // we will only attempt growing 100 times
  int cat_counter = 0; // how many times do we try to generate the left & right splits of a categorical variable
  int cont_counter = 0; // how many times do we try to draw a projection vector
  double c_max = 0.0;
  double tmp_phi = 0.0;
  
  
  while(grow && counter < 100){
    prev_max_depth = max_depth;
    bnv.clear();
    t.get_bots(bnv); // get the bottom nodes
    
    for(size_t l = 0; l < bnv.size(); l++){
      dnx = bnv[l]->get_depth();
      //Rcpp::Rcout << "l = " << l << " depth = " << dnx << std::endl;
      if(dnx > max_depth) max_depth = dnx;
    } // figure out the max depth of the tree
    
    if(max_depth < prev_max_depth){
      Rcpp::Rcout << "max_depth = " << max_depth << " prev_max_depth = " << prev_max_depth << std::endl;
      Rcpp::stop("[draw_tree]: max depth should never be less than prev_max_depth");
    } else if(max_depth > 1 + prev_max_depth){
      Rcpp::Rcout << "max_depth = " << max_depth << " prev_max_depth = " << prev_max_depth << std::endl;
      Rcpp::stop("[draw_tree]: max_depth should not be more than 1 + prev_max_depth");
    } else if(max_depth == prev_max_depth && max_depth != 0){
      Rcpp::Rcout << "max_depth = " << max_depth << " prev_max_depth = " << prev_max_depth << std::endl;
      break;
    }// can no longer grow the tree
    
    // if we get to this point, then max_depth = 1 + prev_max_depth and we attempt to grow the tree
    grow = false;
    for(size_t l = 0; l < bnv.size(); l++){
      dnx = bnv[l]->get_depth();
      if(dnx == max_depth){
        PGnx = tree_pi.alpha/pow(1.0 + (double) dnx, tree_pi.beta);
        if(gen.uniform() < PGnx){
          // we are going to propose a grow move!
          grow = true;
          rule_t rule;
          if(gen.uniform() <= tree_pi.cont_rule_prob){
            // propose a continuous decision rule
            rule.is_cat = false;
            rule.proj_vec.clear();
            rule.c = 0.0;
            rule.v_cat=0;
            rule.l_vals.clear();
            rule.r_vals.clear();
            tmp_phi = 0.0;
            
            cont_counter = 0;
            while( (rule.proj_vec.size() < 1) && (cont_counter < 100)){
              rule.proj_vec.clear();
              c_max = 0.0;
              for(int j = 0; j < di.p_cont; j++){
                if(gen.uniform() < tree_pi.v_cont_probs->at(j)){
                  tmp_phi = gen.uniform(-1.0, 1.0);
                  rule.proj_vec.insert(std::pair<int,double>(j,tmp_phi));
                  c_max += fabs(tmp_phi);
                }
              }
              cont_counter++;
            }
            //Rcpp::Rcout << "  proj depends on" << rule.proj_vec.size() << " variables" << std::endl;
            if(rule.proj_vec.size() < 1 || cont_counter == 100) Rcpp::stop("[draw_tree]: Could't get non-null projection vector");
            if(rule.proj_vec.size() == 1){
              // projection contains only one variable
              //std::map<int,double>::iterator it = rule.proj_vec.begin();
              //it->second = 1.0;
              //Rcpp::Rcout << " projection has support 1" << std::endl;
              rule.proj_vec.begin()->second = 1.0;
              c_max = 1.0;
            }
            rule.c = gen.uniform(-1.0 * c_max, 1.0 * c_max); // draw cutpoint uniformly from appropriate range
          } else{
            rule.is_cat = true;
            rule.proj_vec.clear();
            rule.c = 0.0;
            
            rule.v_cat = gen.multinomial(di.p_cat, tree_pi.v_cat_probs);
            std::set<int> avail_levels = di.cat_levels->at(rule.v_cat);
            bnv[l]->get_rg_cat(rule.v_cat, avail_levels);
            // if there are only trivial rules left, we need to propose a split of all possible levels in the original data (this produces one empty leaf node)
            if(avail_levels.size() <= 1) avail_levels = di.cat_levels->at(rule.v_cat);
            
            cat_counter = 0;
            while( ((rule.l_vals.size() == 0) || (rule.r_vals.size() == 0)) && cat_counter < 100){
              rule.l_vals.clear();
              rule.r_vals.clear();
              for(set_it it = avail_levels.begin(); it != avail_levels.end(); ++it){
                if(gen.uniform() < 0.5) rule.l_vals.insert(*it);
                else rule.r_vals.insert(*it);
              }
              cat_counter++;
            }
            if(cat_counter == 100) Rcpp::stop("Unable to assign at least one level to each child in 100 attempts!"); // this should never get triggered...
          }
          // by this time, we should have a rule ready to go
          t.birth(bnv[l]->get_nid(), rule);
        } // closes if/else checking whether we are continuing to grow or not
      } // closes if checking that we are trying to split a node at maximal depth
    } // closes loop over all bottom nodes
    counter++;
  } // closes while
  
  bnv.clear();
  t.get_bots(bnv);
  for(size_t l = 0; l < bnv.size(); l++) bnv[l]->set_mu(tree_pi.mu0 + tree_pi.tau * gen.normal());
}
*/
