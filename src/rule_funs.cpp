#include "rule_funs.h"


// draw cutpoint for axis-aligned rules
void draw_aa_cutpoint(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen)
{
  double c_upper = 1.0;
  double c_lower = -1.0;
  tree::tree_p nx = t.get_ptr(nid); // at what node are we proposing this rule.
  if(tree_pi.cutpoints->at(rule.v_aa).size() > 0){
    // draw the cutpoint from the supplied cutpoints
    std::set<double> cutpoints = tree_pi.cutpoints->at(rule.v_aa);
    c_lower = *(cutpoints.begin()); // returns smallest element in set
    c_upper = *(cutpoints.rbegin()); // reverse iterator, returns largest value in set
    nx->get_rg_aa(rule.v_aa, c_lower, c_upper);
    if(c_lower >= c_upper){
      // this is a weird tree and we'll just propose a trivial split
      c_lower = *(cutpoints.begin()); // returns smallest element in set
      c_upper = *(cutpoints.rbegin()); // reverse iterator, returns largest value in set
    }
    std::vector<double> valid_cutpoints;
    if(cutpoints.count(c_lower) != 1 || cutpoints.count(c_upper) != 1){
      // c_lower and c_upper were not found in the set of available cutpoints
      Rcpp::Rcout << "[draw_rule]: attempting to select a cutpoint from given set" << std::endl;
      Rcpp::Rcout << "  lower bound is: " << c_lower << " count in set is " << cutpoints.count(c_lower) << std::endl;
      Rcpp::Rcout << "  upper bound is: " << c_upper << " count in set is " << cutpoints.count(c_upper) << std::endl;
      
      Rcpp::Rcout << "rule.v_aa = " << rule.v_aa << " nid = " << nid << std::endl;
      t.print();
      
      Rcpp::stop("we should never have a c that is outside the pre-defined set of cutpoints!");
    }
    // we want to draw from the cutpoints exclusive of c_lower & c_upper;
    // i.e. we want to start with the one just after c_lower and just before c_upper
    // std::set::lower_bound: iterator at first element that is not considered to come before
    // std::set::upper_bound: iterator at first element considered to come after
    // if value is not in set, lower_bound and upper_bound give same result
    // if value is in set: lower bound returns the value, upper bound returns the next value
    for(std::set<double>::iterator it = cutpoints.upper_bound(c_lower); it != cutpoints.lower_bound(c_upper); ++it){
      valid_cutpoints.push_back(*it);
    }
    int num_cutpoints = valid_cutpoints.size();
    if(num_cutpoints < 1){
      // no valid splits are available; we will just pick something, all of the observations will go to one child anyway...
      valid_cutpoints.clear();
      for(std::set<double>::iterator it = cutpoints.begin(); it != cutpoints.end(); ++it){
        valid_cutpoints.push_back(*it);
      }
      num_cutpoints = valid_cutpoints.size();
    }
    // at this point, valid cutpoints is a vector containing the available cutpoints at this node. we pick one uniformly.
    rule.c = valid_cutpoints[floor(gen.uniform() * num_cutpoints)];
  } else{
    // draw cutpoints uniformly
    nx->get_rg_aa(rule.v_raw, c_lower, c_upper);
    if(c_lower >= c_upper){
      c_lower = -1.0;
      c_upper = 1.0;
    }
    rule.c = gen.uniform(c_lower, c_upper);
  }
}


// helper function that partitions categorical levels
void partition_levels(rule_t &rule, std::set<int> &avail_levels, tree_prior_info &tree_pi, RNG &gen)
{
  rule.l_vals.clear();
  rule.r_vals.clear();
  
  if(avail_levels.size() == 1){
    // only one level available at the node. we will send it to the left and all other levels to the right
    rule.l_vals.insert(*avail_levels.begin());
    if(tree_pi.cat_levels->at(rule.v_raw).size() == 1){
      std::set<int> orig_levels = tree_pi.cat_levels->at(rule.v_raw);
      for(std::set<int>::iterator it = orig_levels.begin(); it != orig_levels.end(); ++it){
        if(rule.l_vals.count(*it) == 0) rule.r_vals.insert(*it);
      }
    } else{
      Rcpp::Rcout << "[partition_levels]: Trying categorical rule on X[" << rule.v_raw <<"]";
      Rcpp::Rcout <<" but could not find levels for this variable." << std::endl;
      Rcpp::stop("[partitoin_levels]: ensure variable is categorical and levels are provided!");
    }
  } else if(avail_levels.size() == 2){
    // only two levels available at the node. we will send the first to the left and the other to the right
    rule.l_vals.insert(*avail_levels.begin());
    rule.r_vals.insert(*avail_levels.rbegin());
  } else if(avail_levels.size() > 2){
    if(tree_pi.edges->at(rule.v_raw).size() > 0){
      // this is a network-structured variable
      graph_partition(rule.l_vals, rule.r_vals, tree_pi.edges->at(rule.v_raw), avail_levels, tree_pi.graph_cut_type, gen);
    } else{
      // v is not network structured, so we partition completely at random
      int rule_counter = 0;
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
        Rcpp::Rcout << " number of avail_levels: " << avail_levels.size() << std::endl;
        Rcpp::stop("[partition_levels]: failed to generate valid categorical split in 1000 attempts"); // this should almost surely not get triggered.
      }
    }
  } else{
    Rcpp::stop("[partition_levels]: no available levels provided...");
  }
}


void compute_nested_theta(std::vector<double> &nest_theta, tree &t, int &nid, int &p_cont, int &p_cat, tree_prior_info &tree_pi)
{
  int p = p_cont + p_cat;
  // initialize prob for all continuous variables to be 1 and categorical variables to 0
  //for(int j = 0; j < p_cont; ++j) nest_theta[j] = 1.0;
  for(int j = 0; j < p_cont; ++j){
    if(tree_pi.theta->at(j) != 0.0) nest_theta[j] = 1.0;
  }
  for(int j = p_cont; j < p; ++j) nest_theta[j] = 0.0;
  
  tree::tree_p nx = t.get_ptr(nid);
  // get categorical variables used at rules at nx's ancestors
  std::vector<int> anc_v;
  nx->get_anc_v_cat(anc_v);
  
  // if anc_v is empty, that means we have not yet split on *any* categorical predictors
  if(anc_v.size() == 0){
    //for(int j = p_cont; j < p; ++j) nest_theta[j] = 1.0;
    // loop over clusters
    for(std::map<int, std::set<int>>::iterator c_it = tree_pi.nest_components->begin(); c_it != tree_pi.nest_components->end(); ++c_it){
      if(c_it->second.size() == 1) nest_theta[*(c_it->second.begin()) + p_cont] = 1.0; // singleton cluster
      else{
        if(tree_pi.nest_v_option == 0){
          // allow splits on any variable from the component
          //for(std::set<int>::iterator it = c_it->second.begin(); it != c_it->second.end(); ++it) nest_theta[*it + p_cont] = 1.0;
          for(std::set<int>::iterator it = c_it->second.begin(); it != c_it->second.end(); ++it){
            if(tree_pi.theta->at(*it + p_cont) != 0.0) nest_theta[*it + p_cont] = 1.0;
          }
          
        } else if(tree_pi.nest_v_option == 1){
          // find the lowest possible resolution variables: these have edges coming in but none going out
          for(std::set<int>::iterator it = c_it->second.begin(); it != c_it->second.end(); ++it){
            edge_map_it in_e_it = tree_pi.nest_in->find(*it);
            edge_map_it out_e_it = tree_pi.nest_out->find(*it);
            if(in_e_it->second.size() > 0 && out_e_it->second.size() == 0 && tree_pi.theta->at(*it + p_cont) != 0.0) nest_theta[*it + p_cont] = 1.0;
          }
        } else if(tree_pi.nest_v_option == 2){
          // find the highest possible resolution variables: these have edges going out but none coming in
          for(std::set<int>::iterator it = c_it->second.begin(); it != c_it->second.end(); ++it){
            edge_map_it in_e_it = tree_pi.nest_in->find(*it);
            edge_map_it out_e_it = tree_pi.nest_out->find(*it);
            if(in_e_it->second.size() == 0 && out_e_it->second.size() > 0 && tree_pi.theta->at(*it + p_cont) != 0.0) nest_theta[*it + p_cont] = 1.0;
          }
        } else if(tree_pi.nest_v_option == 3){
          // allow split on anything but the highest or lowest resolution variable in the cluster
          for(std::set<int>::iterator it = c_it->second.begin(); it != c_it->second.end(); ++it){
            edge_map_it in_e_it = tree_pi.nest_in->find(*it);
            edge_map_it out_e_it = tree_pi.nest_out->find(*it);
            if(in_e_it->second.size() > 0 && out_e_it->second.size() > 0 && tree_pi.theta->at(*it + p_cont) != 0.0) nest_theta[*it + p_cont] = 1.0;
          }
        } else{
          Rcpp::stop("[compute_nested_theta]: nest_v_option must be in {0, 1, 2, 3}");
        } // closes if/else's checking nest_v_option
      } // closes if/else checking if cluster is singleton
    } // closes loop over clusters
  } else{
    // we have previously split on a categorical variables. we need to figure out which clusters were used
    std::map<int,int> cluster_rep; // key: id of the cluster in nest_graph; value: variable used at ancestor
    for(std::vector<int>::iterator v_it = anc_v.begin(); v_it != anc_v.end(); ++v_it){
      for(std::map<int, std::set<int>>::iterator c_it = tree_pi.nest_components->begin(); c_it != tree_pi.nest_components->end(); ++c_it){
        if(c_it->second.count(*v_it) == 1){
          if(cluster_rep.count(c_it->first) == 0){
            // this is the first time a variable from the cluster is present. we need to add *v_it to the set of representatives
            cluster_rep.insert(std::pair<int,int>(c_it->first, *v_it));
            break;
          }
          break;
        }
      }
    }
    // now loop over all the clusters (i.e. keys of nest_components)
    // if the cluster is not in cluster_rep: decide what to do based on nest_v_option
    // otherwise, we check if cluster is singleton: if it is, we set corresponding theta = 1
    // if not a singleton, we set based on nest_v_option
    for(std::map<int, std::set<int>>::iterator c_it = tree_pi.nest_components->begin(); c_it != tree_pi.nest_components->end(); ++c_it){
      if(cluster_rep.count(c_it->first) == 0){
        // did not previously split on variable from this component.
        if(c_it->second.size() == 1 && tree_pi.theta->at(*(c_it->second.begin()) + p_cont) != 0.0) nest_theta[*(c_it->second.begin()) + p_cont] = 1.0; // singleton cluster
        else{
          if(tree_pi.nest_v_option == 0){
            // allow splits on any variable from the component
            for(std::set<int>::iterator it = c_it->second.begin(); it != c_it->second.end(); ++it){
              if(tree_pi.theta->at(*it + p_cont) != 0.0) nest_theta[*it + p_cont] = 1.0;
            }
          } else if(tree_pi.nest_v_option == 1){
            // find the lowest possible resolution variables: these have edges coming in but none going out
            for(std::set<int>::iterator it = c_it->second.begin(); it != c_it->second.end(); ++it){
              edge_map_it in_e_it = tree_pi.nest_in->find(*it);
              edge_map_it out_e_it = tree_pi.nest_out->find(*it);
              if(in_e_it->second.size() > 0 && out_e_it->second.size() == 0 && tree_pi.theta->at(*it + p_cont) != 0.0) nest_theta[*it + p_cont] = 1.0;
            }
          } else if(tree_pi.nest_v_option == 2){
            // find the highest possible resolution variables: these have edges going out but none coming in
            for(std::set<int>::iterator it = c_it->second.begin(); it != c_it->second.end(); ++it){
              edge_map_it in_e_it = tree_pi.nest_in->find(*it);
              edge_map_it out_e_it = tree_pi.nest_out->find(*it);
              if(in_e_it->second.size() == 0 && out_e_it->second.size() > 0 && tree_pi.theta->at(*it + p_cont) != 0.0) nest_theta[*it + p_cont] = 1.0;
            }
          } else if(tree_pi.nest_v_option == 3){
            // allow split on anything but the highest or lowest resolution variable in the cluster
            for(std::set<int>::iterator it = c_it->second.begin(); it != c_it->second.end(); ++it){
              edge_map_it in_e_it = tree_pi.nest_in->find(*it);
              edge_map_it out_e_it = tree_pi.nest_out->find(*it);
              if(in_e_it->second.size() > 0 && out_e_it->second.size() > 0 && tree_pi.theta->at(*it + p_cont) != 0.0) nest_theta[*it + p_cont] = 1.0;
            }
          } else{
            Rcpp::stop("[compute_nested_theta]: nest_v_option must be in {0, 1, 2, 3}");
          } // closes if/else's checking nest_v_option
        }
      } else{
        if(c_it->second.size() == 1) nest_theta[*(c_it->second.begin()) + p_cont] = 1.0; // singleton cluster
        else{
          int v_rep = cluster_rep.find(c_it->first)->second; // most recent variable from cluster
          if(tree_pi.nest_v_option == 0){
            // only allow ourselves to split on v_rep
            if(tree_pi.theta->at(v_rep+p_cont) != 0.0) nest_theta[v_rep + p_cont] = 1.0;
          } else if(tree_pi.nest_v_option == 1){
            // include v_rep and everything that has edge into v_rep
            if(tree_pi.theta->at(v_rep + p_cont) != 0.0) nest_theta[v_rep + p_cont] = 1.0;
            edge_map_it e_it = tree_pi.nest_in->find(v_rep);
            for(std::vector<edge>::iterator it = e_it->second.begin(); it != e_it->second.end(); ++it){
              if(tree_pi.theta->at(it->source + p_cont) != 0.0) nest_theta[it->source + p_cont] = 1.0;
            }
          } else if(tree_pi.nest_v_option == 2){
            // include v_rep and everything that has edge out of v_rep
            nest_theta[v_rep + p_cont] = 1.0;
            edge_map_it e_it = tree_pi.nest_out->find(v_rep);
            for(std::vector<edge>::iterator it = e_it->second.begin(); it != e_it->second.end(); ++it){
              if(tree_pi.theta->at(it->sink + p_cont) != 0.0) nest_theta[it->sink + p_cont] = 1.0;
            }
          } else if(tree_pi.nest_v_option == 3){
            // include v_rep and everything has edge out of v_rep or into v
            nest_theta[v_rep + p_cont] = 1.0;
            edge_map_it e_it = tree_pi.nest_in->find(v_rep);
            for(std::vector<edge>::iterator it = e_it->second.begin(); it != e_it->second.end(); ++it){
              if(tree_pi.theta->at(it->source + p_cont) != 0.0) nest_theta[it->source + p_cont] = 1.0;
            }
            e_it = tree_pi.nest_out->find(v_rep);
            for(std::vector<edge>::iterator it = e_it->second.begin(); it != e_it->second.end(); ++it){
              if(tree_pi.theta->at(it->sink + p_cont) != 0.0) nest_theta[it->sink + p_cont] = 1.0;
            }
          } else{
            Rcpp::stop("[compute_nested_theta]: nest_v_option must be in {0, 1, 2, 3}");
          } // closes if/else's checking nest_v_option
        } // closes if/else checking that the cluster is not a singleton
      } // closes if/else checking if any ancestors split on variable from this cluster
    } // closes loop over the components in the nest graph
  } // closes if/else checking that anc_v is empty or not
}

void draw_rule(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen)
{
  rule.clear();
  int v_raw = 0;
  if(tree_pi.nest_v){
    std::vector<double> nest_theta(di.p, 0.0);
    compute_nested_theta(nest_theta, t, nid, di.p_cont, di.p_cat, tree_pi);
    v_raw = gen.categorical(nest_theta);
  } else{
    v_raw = gen.categorical(tree_pi.theta);
  }
  if(v_raw < di.p_cont){
    rule.is_cat = false;
    rule.v_aa = v_raw;
    draw_aa_cutpoint(rule, t, nid, di, tree_pi, gen);
  } else{
    rule.is_cat = true;
    rule.v_cat = v_raw - di.p_cont;
    std::set<int> avail_levels;
    if(tree_pi.nest_c){
      t.get_ptr(nid)->get_rg_nested_cat(avail_levels, rule.v_cat, tree_pi);
    } else{
      t.get_ptr(nid)->get_rg_cat(avail_levels, rule.v_cat);
    }
    if(avail_levels.size() <= 1) avail_levels = tree_pi.cat_levels->at(rule.v_cat);
    partition_levels(rule, avail_levels, tree_pi, gen);
  }
}
