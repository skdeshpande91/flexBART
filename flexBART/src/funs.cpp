#include "funs.h"

void tree_traversal(suff_stat &ss, tree &t, data_info &di)
{
  double* xx_cont = 0;
  int* xx_cat = 0;
  tree::tree_cp bn;
  int nid;
  ss.clear(); // clear out the sufficient statistic map
  
  tree::npv bnv;
  t.get_bots(bnv);
  suff_stat_it ss_it; // used to look up to which element of ss a particular bottom node corresponds
  
  // add an element to suff stat map for each bottom node
  for(tree::npv_it it = bnv.begin(); it != bnv.end(); ++it){
    nid = (*it)->get_nid(); // bnv is a vector of pointers. it points to elements of bnv so we need (*it) to access the members of these elements
    // add element to sufficient stats map with key = nid, value = empty vector to hold index of observations assigned to this node
    ss.insert(std::pair<int,std::vector<int>>(nid,std::vector<int>()));
  }
  
  // ready to loop over all observations
  for(int i = 0; i < di.n; i++){
    if(di.x_cont != 0) xx_cont = di.x_cont + i * di.p_cont;
    if(di.x_cat != 0) xx_cat = di.x_cat + i * di.p_cat;
    bn = t.get_bn(xx_cont, xx_cat);
    if(bn == 0){
      Rcpp::Rcout << "i = " << i << std::endl;
      t.print();
      Rcpp::stop("[tree_traversal]: could not find bottom node!");
    }
    else{
      nid = bn->get_nid();
      if(ss.count(nid) != 1) Rcpp::stop("[tree_traversal]: bottom node not included in sufficient statistic map!"); // should never be encountered
      else{
        ss_it = ss.find(nid); // iterator now set to element of ss corresponding to the bottom node holding observation i
        ss_it->second.push_back(i);
      } // closes if/else checking that i-th observation's bottom node was in our map
    } // closes if/else checking that i-th observation maps to valid bottom node
  } // closes loop over all observation
}


void fit_ensemble(std::vector<double> &fit, std::vector<tree> &t_vec, data_info &di){
  if(fit.size() != di.n) Rcpp::stop("[fit_ensemble]: size of fit must be equal to di.n!"); // honestly should never get triggered
  double* xx_cont = 0;
  int* xx_cat = 0;
  for(int i = 0; i < di.n; i++){
    if(di.x_cont != 0) xx_cont = di.x_cont + i * di.p_cont;
    if(di.x_cat != 0) xx_cat = di.x_cat + i * di.p_cat;
    fit[i] = 0.0;
    for(int m = 0; m < t_vec.size(); m++) fit[i] += t_vec[m].evaluate(xx_cont, xx_cat); // involves a tree traversal
  }
}

std::string write_tree(tree &t, tree_prior_info &tree_pi, set_str_conversion &set_str)
{
  std::ostringstream os;
  os.precision(32);
  
  tree::cnpv nds;
  rule_t rule;
  t.get_nodes(nds);
  
  for(tree::cnpv_it nd_it = nds.begin(); nd_it != nds.end(); ++nd_it){
    os << (*nd_it)->get_nid() << " ";
    if( (*nd_it)->get_ntype() == 'b'){
      // it's a bottom node
      os << "m " << (*nd_it)->get_mu() << " " << std::endl; // m tells us to expect mu next
    } else if( ((*nd_it)->get_ntype() == 't') && ( !(*nd_it)->l) ){ // because we need to look at left child of a node, make write_tree a friend in tree class
      // it's the top node and it has no children
      os << "m " << (*nd_it)->get_mu() << " " << std::endl; // tree is a stump, m tells us to expect mu next
    } else{
      // we need to print out the rule
      //os << "make rule " << std::endl;
      os << "r "; // node has a decision rule. r tells us to expect a rule next
      rule.clear();
      rule = (*nd_it)->get_rule();
      

      // os << rule.is_cat << " " << rule.is_rc << " ";
      os << rule.is_cat << " ";

      if(!rule.is_cat){
        // axis-aligned rule
        os << rule.c << " " << rule.v_aa;
      } else{
        // categorical rule
        int K = tree_pi.K->at(rule.v_cat); // how many levels
        os << rule.v_cat << " " << K << " ";
        os << set_str.set_to_hex(K, rule.l_vals) << " ";
        os << set_str.set_to_hex(K, rule.r_vals) << " ";
      }
      os << std::endl;
    } // closes if/else checking what type of node we are writing
  } // closes loop over the nodes in t
  
  return os.str();
  
}

void read_tree(tree &t, std::string &tree_string, set_str_conversion &set_str)
{
  std::istringstream tree_ss(tree_string); // an in stringstream of the tree's string representation
  std::string node_string; // string for each individual node in the tree
  
  int nid;
  char stream_type; // either 'm' to indicate that next element in stream is mu or 'r' to indicate a rule follows
  
  double tmp_mu; // holds the value of mu for a leaf node

  //char aa; // '0' or '1' for rule.is_aa
  char cat; // '0' or '1' for rule.is_cat
  rule_t tmp_rule; // temporary rule that gets populated as we read the tree's string/stream
  
  int tmp_v; // for reading in rc weights
  double tmp_phi; // for reading in rc weights
  
  int K; // tells us how many levels there were to the categorical variable
  std::string l_hex; // string representation of the l_vals in a categorical rule
  std::string r_hex; // string representation of the r_vals in a categorical rule
  
  std::map<int, rule_t> decision_nodes;
  std::map<int, double> leaf_nodes;
  
  while(tree_ss){
    std::getline(tree_ss, node_string, '\n');
    if(node_string.size() > 0){
      std::istringstream node_ss(node_string); // in stream for the single node
      node_ss >> nid; // get the node nid
      node_ss >> stream_type;
      
      if(stream_type == 'm'){
        node_ss >> tmp_mu;
        leaf_nodes.insert(std::pair<int,double>(nid, tmp_mu));
      } else if(stream_type == 'r'){
        tmp_rule.clear();
        node_ss >> cat;
        
        if(cat == '0') tmp_rule.is_cat = false;
        else tmp_rule.is_cat = true;

        
        if(!tmp_rule.is_cat){
          // axis-aligned
          node_ss >> tmp_rule.c;
          node_ss >> tmp_rule.v_aa;
        } else{
          // categorical rule
          node_ss >> tmp_rule.v_cat; // get the variable index
          node_ss >> K; // we now know how many levels of the categorical variable there were
          node_ss >> l_hex;
          node_ss >> r_hex;
          
          tmp_rule.l_vals = set_str.hex_to_set(K, l_hex);
          tmp_rule.r_vals = set_str.hex_to_set(K, r_hex);
        }
        decision_nodes.insert(std::pair<int, rule_t>(nid, tmp_rule));
      } // closes if/else checking what type of node we're parsing
    } // closes if checking that we found a valid node in the stream
  } // closes while that parses stream for the tree
  
  // we now have decision_nodes and leaf_nodes and are ready to build up our tree
  t.to_null(); // clear out the tree if there was anything there
  
  // remember std::map is sorted by key.
  // we have always used node id as the key so by iterating over our map, we will *never*
  // attempt to birth from a node that has not already been created.

  for(std::map<int,rule_t>::iterator it = decision_nodes.begin(); it != decision_nodes.end(); ++it){
    t.birth(it->first, it->second); // do the birth.
  }

  tree::npv bnv;
  t.get_bots(bnv); // get the bottom nodes
  std::map<int,double>::iterator leaf_it;

  for(tree::npv_it it = bnv.begin(); it != bnv.end(); ++it){
    leaf_it = leaf_nodes.find( (*it)->get_nid() );
    if(leaf_it == leaf_nodes.end()){
      Rcpp::Rcout << "[read_tree]: we didn't read a leaf node with nid" << (*it)->get_nid() << std::endl;
      Rcpp::stop("mistake in reading in tree");
    } else{
      (*it)->set_mu(leaf_it->second);
    }
  }
  
}



void draw_rule(rule_t &rule, tree &t, int &nid, data_info &di, tree_prior_info &tree_pi, RNG &gen){
  rule.clear(); // clear out the rule
  // we now are ready to draw a decision rule

  int rule_counter = 0; // we are allowed multiple tries to draw a valid random combination or categorical rule
  double c_upper = 1.0; // upper bound for range of cutpoints in axis aligned split
  double c_lower = -1.0; // lower bound for range of cutpoints in axis aligned split
  double tmp_weight = 0.0; // weights of random combination
  double c_max = 1.0; // upper bound for absolute value of cutpoint in random combination split
  tree::tree_p nx = t.get_ptr(nid); // at what node are we proposing this rule.
  
  double unif = gen.uniform();
  if( (!tree_pi.rc_split) || (tree_pi.rc_split && unif > *tree_pi.prob_rc) ){
    // either we were never allowed to try a rc split OR (1) we are allowed to try rc splits but (2) randomly choose a different type of rule
    int v_raw = gen.multinomial(di.p, tree_pi.theta);
    if(v_raw < di.p_cont){
      // continuous variable so it's an axis-aligned splitting rule
      rule.is_aa = true;
      rule.is_cat = false;
      rule.v_aa = v_raw;
      if(tree_pi.unif_cuts[rule.v_aa] == 0){
        // draw the cutpoint from the supplied cutpoints
        c_lower = *(tree_pi.cutpoints->at(rule.v_aa).begin()); // returns smallest element in set
        c_upper = *(tree_pi.cutpoints->at(rule.v_aa).rbegin()); // reverse iterator, returns largest value in set
        nx->get_rg_aa(rule.v_aa, c_lower, c_upper);
        if(c_lower >= c_upper){
          // this is a weird tree and we'll just propose a trivial split
          c_lower = *(tree_pi.cutpoints->at(rule.v_aa).begin());
          c_upper = *(tree_pi.cutpoints->at(rule.v_aa).rbegin());
        }
        std::vector<double> valid_cutpoints;
        if(tree_pi.cutpoints->at(rule.v_aa).count(c_lower) != 1 || tree_pi.cutpoints->at(rule.v_aa).count(c_upper) != 1){
          // c_lower and c_upper were not found in the set of available cutpoints
          Rcpp::Rcout << "[draw_rule]: attempting to select a cutpoint from given set" << std::endl;
          Rcpp::Rcout << "  lower bound is: " << c_lower << " count in set is " << tree_pi.cutpoints->at(rule.v_aa).count(c_lower) << std::endl;
          Rcpp::Rcout << "  upper bound is: " << c_upper << " count in set is " << tree_pi.cutpoints->at(rule.v_aa).count(c_upper) << std::endl;
          //Rcpp::Rcout << "  cutpoints are:";
          //for(std::set<double>::iterator it = tree_pi.cutpoints->at(rule.v_aa).begin(); it != tree_pi.cutpoints->at(rule.v_aa).end(); ++it) Rcpp::Rcout << " " << *it;
          //Rcpp::Rcout << std::endl;
          Rcpp::stop("we should never have a c that is outside the pre-defined set of cutpoints!");
        }
        // we want to draw from the cutpoints exclusive of c_lower & c_upper;
        // i.e. we want to start with the one just after c_lower and just before c_upper
        // std::set::lower_bound: iterator at first element that is not considered to come before
        // std::set::upper_bound: iterator at first element considered to come after
        // if value is not in set, lower_bound and upper_bound give same result
        // if value is in set: lower bound returns the value, upper bound returns the next value
        for(std::set<double>::iterator it = tree_pi.cutpoints->at(rule.v_aa).upper_bound(c_lower); it != tree_pi.cutpoints->at(rule.v_aa).lower_bound(c_upper); ++it){
          valid_cutpoints.push_back(*it);
        }
        int num_cutpoints = valid_cutpoints.size();
        if(num_cutpoints < 1){
          // no valid splits are available; we will just pick something, all of the observations will go to one child anyway...
          valid_cutpoints.clear();
          for(std::set<double>::iterator it = tree_pi.cutpoints->at(rule.v_aa).begin(); it != tree_pi.cutpoints->at(rule.v_aa).end(); ++it){
            valid_cutpoints.push_back(*it);
          }
          num_cutpoints = valid_cutpoints.size();
        }
        // at this point, valid cutpoints is a vector containing the available cutpoints at this node. we pick one uniformly.
        rule.c = valid_cutpoints[floor(gen.uniform() * num_cutpoints)];
        
      } else{
        // draw cutpoints uniformly
        c_upper = 1.0;
        c_lower = -1.0;
        nx->get_rg_aa(rule.v_aa, c_lower, c_upper);
        if(c_lower >= c_upper){
          c_lower = -1.0;
          c_upper = 1.0;
        }
        rule.c = gen.uniform(c_lower, c_upper);
      }
    } else{
      // categorical decision rule
      rule.is_aa = false;
      rule.is_cat = true;
      rule.v_cat = v_raw - di.p_cont;

      std::set<int> avail_levels = tree_pi.cat_levels->at(rule.v_cat); // get the full set of levels for this variable
      nx->get_rg_cat(rule.v_cat, avail_levels); // determine the set of levels available at nx.
      // if there is only one level left for this variable at nx, we will just propose a trivial split
      // and will reset the value of avail_levels to be the full set of all levels for the variable
      if(avail_levels.size() <= 1) avail_levels = tree_pi.cat_levels->at(rule.v_cat);
      
      rule.l_vals.clear();
      rule.r_vals.clear();
      
      if(tree_pi.graph_split[rule.v_cat] == 1 && tree_pi.edges->at(rule.v_cat).size() > 0){
        // if we explicitly say to use the graph to split the variables
        //graph_partition(avail_levels, rule.l_vals, rule.r_vals, tree_pi.edges->at(rule.v_cat), tree_pi.K->at(rule.v_cat), tree_pi.graph_cut_type, gen);
        graph_partition(rule.l_vals, rule.r_vals, tree_pi.edges->at(rule.v_cat), avail_levels, tree_pi.graph_cut_type, gen);
      } else{
        // otherwise we default to splitting the available levels uniformly at random: prob 0.5 to go to each child
        rule_counter = 0;
        double tmp_prob = 0.5;
        if(tree_pi.a_cat > 0  && tree_pi.b_cat > 0) tmp_prob = gen.beta(tree_pi.a_cat, tree_pi.b_cat);
        else tmp_prob = 0.5;
        
        while( ((rule.l_vals.size() == 0) || (rule.r_vals.size() == 0)) && rule_counter < 1000 ){
          rule.l_vals.clear();
          rule.r_vals.clear();
          for(set_it it = avail_levels.begin(); it != avail_levels.end(); ++it){
            //if(gen.uniform() <= 0.5) rule.l_vals.insert(*it);
            if(gen.uniform() <= tmp_prob) rule.l_vals.insert(*it);
            else rule.r_vals.insert(*it);
          }
          ++(rule_counter);
        }
        if(rule_counter == 1000){
          Rcpp::stop("failed to generate valid categorical split in 1000 attempts"); // this should almost surely not get triggered.
        }
      }
      if( (rule.l_vals.size() == 0) || (rule.r_vals.size() == 0) ){
        Rcpp::Rcout << "[draw_rule]: proposed an invalid categorical rule" << std::endl;
        Rcpp::stop("proposed an invalid categorical rule!");
      }
    } // closes if/else determining whether we do an axis-aligned or categorical decision rule
  } else if(tree_pi.rc_split && unif <= *tree_pi.prob_rc){
    // random combination rule
    rule.is_aa = false;
    rule.is_cat = false;

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
    Rcpp::stop("[draw_rule]: unable to draw a rule. check values of rc_split & prob_rc");
  } // closes all if/else determining the type of rule (axis-aligned or categorical) OR random combination
}


void update_theta_u(std::vector<double> &theta, double &u, std::vector<int> &var_count, int &p, double &a_u, double &b_u, RNG &gen)
{
  if(theta.size() != p){
    Rcpp::Rcout << "theta has size " << theta.size() << "  p = " << p << std::endl;
    Rcpp::stop("theta must have size p!");
  } else{
    double tmp_sum = 0.0;
    double tmp_concentration = 0.0;
    double sum_log_theta = 0.0;
    int v_count;
    std::vector<double> tmp_gamma(p, 0.0);
    
    // update theta first
    double u_orig = u;
    for(int j = 0; j < p; j++){
      v_count = var_count[j];
      tmp_concentration = u_orig/(1.0 - u_orig) + (double) v_count;
      tmp_gamma[j] = gen.gamma(tmp_concentration, 1.0);
      tmp_sum += tmp_gamma[j];
    }
    for(int j = 0; j < p; j++){
      theta[j] = tmp_gamma[j]/tmp_sum;
      sum_log_theta += log(theta[j]);
    }
    
    // we're now ready to update u
    double u_prop = gen.beta(a_u,b_u);
    double log_like_prop = (u_prop)/(1.0 - u_prop) * sum_log_theta;
    double log_like_orig = (u_orig)/(1.0 - u_orig) * sum_log_theta;
    
    log_like_prop += lgamma( (double) p * u_prop/(1.0 - u_prop)) - ((double) p) * lgamma(u_prop/(1.0 - u_prop));
    log_like_orig += lgamma( (double) p * u_orig/(1.0 - u_orig)) - ((double) p) * lgamma(u_orig/(1.0 - u_orig));
    double log_accept = log_like_prop - log_like_orig;
    if(gen.log_uniform() <= log_accept) u = u_prop;
    else u = u_orig;
  }
}


