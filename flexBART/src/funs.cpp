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
    // add element to sufficient stats map with key = nid
    // value is a vector of integers of length di.n recording which observations are assigned to this bottom node
    ss.insert(std::pair<int,std::vector<int>>(nid,std::vector<int>()));
  }
  
  /*
  // for debugging purposes
  Rcpp::Rcout << "[tree_traversal]: added elements to suff. stat. map for each bottom node!" << std::endl;
  for(suff_stat_it ss_it = ss.begin(); ss_it != ss.end(); ++ss_it){
    Rcpp::Rcout << "  nid = " << ss_it->first << "  n = " << ss_it->second.size() << std::endl;
  }
  Rcpp::Rcout << "[tree_traversal]: ready to start traversal!" << std::endl;
  */
  
  // ready to loop over all observations
  for(int i = 0; i < di.n; i++){
    if(di.x_cont != 0) xx_cont = di.x_cont + i * di.p_cont;
    if(di.x_cat != 0) xx_cat = di.x_cat + i * di.p_cat;
    bn = t.get_bn(xx_cont, xx_cat);
    if(bn == 0) Rcpp::stop("[tree_traversal]: could not find bottom node!"); // should never be encountered
    else{
      nid = bn->get_nid();
      if(ss.count(nid) != 1) Rcpp::stop("[tree_traversal]: bottom node not included in sufficient statistic map!"); // should never be encountered
      else{
        ss_it = ss.find(nid); // iterator now set to element of ss corresponding to the bottom node holding observation i
        ss_it->second.push_back(i);
      }
    } // closes if/else checking that i-th observation's bottom node was in our map
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

void compute_suff_stat_grow(suff_stat &orig_suff_stat, suff_stat &new_suff_stat, int &nx_nid, rule_t &rule, tree &t, data_info &di)
{
  double* xx_cont;
  int* xx_cat;
  int i;
  double tmp_x;
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
  
  for(int_it it = nx_it->second.begin(); it != nx_it->second.end(); ++it){
    i = *it;
    if(di.x_cont != 0) xx_cont = di.x_cont + i * di.p_cont;
    if(di.x_cat != 0) xx_cat = di.x_cat + i * di.p_cat;
    
    if(!rule.is_cat && !rule.is_rc){
      // axis-aligned rule
      if(xx_cont[rule.v_aa] < rule.c) nxl_it->second.push_back(i);
      else if(xx_cont[rule.v_aa] >= rule.c) nxr_it->second.push_back(i);
      else Rcpp::stop("[compute_ss_grow]: could not assign observation to left or right child");
      
    } else if(!rule.is_cat && rule.is_rc){
      // random combination rule
      tmp_x = 0.0;
      for(rc_it rcit = rule.rc_weight.begin(); rcit != rule.rc_weight.end(); ++rcit) tmp_x += (rcit->second) * xx_cont[rcit->first];
      if(tmp_x < rule.c) nxl_it->second.push_back(i);
      else if(tmp_x >= rule.c) nxr_it->second.push_back(i);
      else Rcpp::stop("[compute_ss_grow]: could not assign observation to left or right child");
    } else if(rule.is_cat && !rule.is_rc){
      // categorical rule
      // we need to see whether i-th observation's value of the categorical pred goes to left or right
      // std::set.count returns 1 if the value is in the set and 0 otherwise
      l_count = rule.l_vals.count(xx_cat[rule.v_cat]);
      r_count = rule.r_vals.count(xx_cat[rule.v_cat]);
      if(l_count == 1 && r_count == 0){
        nxl_it->second.push_back(i);
      } else if(l_count == 0 && r_count == 1){
        nxr_it->second.push_back(i);
      } else if(l_count == 1 && r_count == 1){
        Rcpp::stop("[compute_ss_grow]: observation goes to both left & right child...");
      } else{
        Rcpp::stop("[compute_ss_grow]: observation doesn't go to either left or right child");
      }
    } else{
      // we should never hit this error
      Rcpp::stop("[compute_ss_grow]: cannot resolve the type of decision rule");
    }
  } // closes loop over all entries in nx
  
}

void compute_suff_stat_prune(suff_stat &orig_suff_stat, suff_stat &new_suff_stat, int &nl_nid, int &nr_nid, int &np_nid, tree &t, data_info &di)
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
  /*
  for(int_it it = nl_it->second.begin(); it != nl_it->second.end(); ++it){
    i = *it;
    np_it->second.push_back(i); // could probably get away with *it but let's be safe
  }
  // now let's add the elements from nr_it
  for(int_it it = nr_it->second.begin(); it != nr_it->second.end(); ++it){
    i = *it;
    np_it->second.push_back(i);
  }
  */
}

double compute_lil(suff_stat &ss, int &nid, double &sigma, data_info &di, tree_prior_info &tree_pi)
{
  // reminder posterior of jump mu is N(P^-1 Theta, P^-1)
  //int i;
  if(ss.count(nid) != 1) Rcpp::stop("[compute_lil]: did not find node in suff stat map!");
  suff_stat_it ss_it = ss.find(nid);
  
  double P = 1.0/pow(tree_pi.tau, 2.0) + ( (double) ss_it->second.size())/pow(sigma, 2.0); // precision of jump mu
  double Theta = tree_pi.mu0/pow(tree_pi.tau, 2.0);
  
  for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it) Theta += di.rp[*it]/pow(sigma, 2.0);
  
  /*
  for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
    i = *it;
    Theta += di.rp[i]/pow(sigma,2.0);
  }
   */
  return(-0.5 * log(P) + 0.5 * pow(Theta,2.0) / P);
  
}

void draw_mu(tree &t, suff_stat &ss, double &sigma, data_info &di, tree_prior_info &tree_pi, RNG &gen)
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
      P = 1.0/pow(tree_pi.tau, 2.0) + ( (double) ss_it->second.size())/pow(sigma, 2.0); // precision of jump mu
      Theta = tree_pi.mu0/pow(tree_pi.tau, 2.0);
      for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it) Theta += di.rp[*it]/pow(sigma,2.0);
      /*
      for(int_it it = ss_it->second.begin(); it != ss_it->second.end(); ++it){
        i = *it;
        Theta += di.rp[i]/pow(sigma,2.0);
      }
       */
      post_sd = sqrt(1.0/P);
      post_mean = Theta/P;
      bn->set_mu(gen.normal(post_mean, post_sd));
    }
  }
}

std::string write_tree(tree &t, data_info &di, set_str_conversion &set_str)
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

      os << rule.is_cat << " " << rule.is_rc << " ";
      if(!rule.is_cat && !rule.is_rc){
        // axis-aligned rule
        os << rule.c << " " << rule.v_aa;
      } else if(!rule.is_cat && rule.is_rc){
        os << rule.c << " ";
        for(rc_it rcit = rule.rc_weight.begin(); rcit != rule.rc_weight.end(); ++rcit){
          os << rcit->first << " " << rcit->second << " ";
        }
      } else if(rule.is_cat){
        int K = di.K->at(rule.v_cat); // how many levels
        os << rule.v_cat << " " << K << " ";
        os << set_str.set_to_hex(K, rule.l_vals) << " ";
        os << set_str.set_to_hex(K, rule.r_vals) << " ";
      }
      os << std::endl;
    } // closes if/else checking what type of node we are writing
  } // closes loop over the nodes in t
  
  return os.str();
  
}

void read_tree(tree &t, std::string &tree_string, data_info &di, set_str_conversion &set_str)
{
  std::istringstream tree_ss(tree_string); // an in stringstream of the tree's string representation
  std::string node_string; // string for each individual node in the tree
  
  int nid;
  char stream_type; // either 'm' to indicate that next element in stream is mu or 'r' to indicate a rule follows
  
  double tmp_mu; // holds the value of mu for a leaf node

  
  char cat; // '0' or '1' for rule.is_cat
  char rc; // '0' or '1' for rule.is_rc
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
        node_ss >> rc;
        if(cat == '0') tmp_rule.is_cat = false;
        else tmp_rule.is_cat = true;
        
        if(rc == '0') tmp_rule.is_rc = false;
        else tmp_rule.is_rc = true;
        
        if(!tmp_rule.is_cat && !tmp_rule.is_rc){
          // axis-aligned decision rule
          node_ss >> tmp_rule.c;
          node_ss >> tmp_rule.v_aa;
        } else if(!tmp_rule.is_cat && tmp_rule.is_rc){
          // random combination rule
          node_ss >> tmp_rule.c; // get the cutpoint first
          while(node_ss){
            node_ss >> tmp_v;
            node_ss >> tmp_phi;
            tmp_rule.rc_weight.insert(std::pair<int,double>(tmp_v, tmp_phi));
          }
        } else if(tmp_rule.is_cat && !tmp_rule.is_rc){
          // categorical rule
          node_ss >> tmp_rule.v_cat; // get the variable index
          node_ss >> K; // we now know how many levels of the categorical variable there were
          node_ss >> l_hex;
          node_ss >> r_hex;
          
          if(K != di.K->at(tmp_rule.v_cat)){
            Rcpp::Rcout << "v_cat = " << tmp_rule.v_cat << std::endl;
            Rcpp::Rcout << "Read in K = " << K << " total categorical levels" << std::endl;
            Rcpp::Rcout << "Corresponding entry of di.K has length " << di.K->at(tmp_rule.v_cat) << std::endl;
            Rcpp::stop("mismatch in number of levels recorded in cat_levels & in saved tree!");
          }
          
          tmp_rule.l_vals = set_str.hex_to_set(K, l_hex);
          tmp_rule.r_vals = set_str.hex_to_set(K, r_hex);
        } else{
          Rcpp::Rcout << "cannot resolve the type of rule" << std::endl;
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
  
  // since we're messing with private members of tree, do we need to make this function a friend of tree? couldn't hurt.
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

// depth-first search, used to find connected components of a graph
void dfs(int i, std::vector<bool> &visited, std::vector<int> &comp, int &n, arma::mat &A)
{
  visited[i] = true; // dfs has reached i for the first time so mark it
  comp.push_back(i); // now that i has been marked, add it to the connected component
  for(int ii = 0; ii < n; ii++){
    if( std::fabs(A(i,ii)) > 1e-16 ){ // in case A is a weighted matrix
      // i is connected to ii
      if(!visited[ii]){
        // somehow ii hasn't been visited before, so we need to continue our dfs from there.
        dfs(ii, visited, comp, n, A);
      }
    }
  }
}

// find connected components of a graph
void find_components(std::vector<std::vector<int> > &components, arma::mat &A)
{
  components.clear(); // clear it out
  int n = A.n_rows;
  std::vector<bool> visited(n,false);
  for(int i = 0; i < n; i++){
    if(!visited[i]){
      std::vector<int> new_comp;
      dfs(i, visited, new_comp, n, A);
      components.push_back(new_comp);
    }
  }
}

// find minimum edge weight in Boruvka's algorithm
std::pair<int,int> find_min_edge_weight(std::vector<int> &components, int &n, arma::mat &W)
{
  if(components.size() == n) Rcpp::stop("[add_min_weight_edge]: component already contains all n nodes!");
  
  std::vector<unsigned int> tmp_in;
  std::vector<unsigned int> tmp_out;
  
  std::vector<int>::iterator find_it;
  for(int i = 0; i < n; i++){
    find_it = std::find(components.begin(), components.end(), i);
    if(find_it != components.end()) tmp_in.push_back( (unsigned int) i);
    else tmp_out.push_back( (unsigned int) i);
  }
  arma::uvec in_index(tmp_in);
  arma::uvec out_index(tmp_out);
  
  arma::mat tmp_W = W(in_index, out_index);
  // check that there are actually edges from component to rest of the graph
  if(arma::all(arma::abs(arma::vectorise(tmp_W)) < 1e-16)) Rcpp::stop("tmpW contains all 0's");
  tmp_W.elem(arma::find(arma::abs(tmp_W) < 1e-16)).ones(); // just to be extra safe, convert all 0's into 1's
  
  
  arma::uword min_index = tmp_W.index_min(); //
  arma::uvec min_sub = arma::ind2sub(arma::size(tmp_W), min_index);
  std::pair<int,int> results( in_index(min_sub(0)), out_index(min_sub(1)) );
  return results;
}
// implement's Boruvka's algorithm
arma::mat boruvka(arma::mat &W)
{
  int n = W.n_rows;
  
  std::vector<std::vector<int> > W_components;
  find_components(W_components, W);
  if(W_components.size() > 1){
    Rcpp::stop("W is not connected!");
  }
  
  arma::mat A_mst = arma::zeros<arma::mat>(n,n); // initialize the
  std::vector<std::vector<int> > mst_components;
  find_components(mst_components, A_mst);
  
  int counter = 0; // failsafe which should *never* be triggered since graph is connected
  while( (mst_components.size() > 1) && (counter < 100) ){
    for(std::vector<std::vector<int> >::iterator it = mst_components.begin(); it != mst_components.end(); ++it){
      
      std::pair<int,int> min_edge = find_min_edge_weight(*it,n, W);
      A_mst(min_edge.first, min_edge.second) = 1;
      A_mst(min_edge.second, min_edge.first) = 1;
    }
    find_components(mst_components, A_mst);
    ++counter;
  }
  
  if(mst_components.size() > 1){
    Rcpp::Rcout << "Was not able to find a spanning tree (check connectivity of W!)" << std::endl;
  }
  return A_mst;
}

void get_edge_probs(std::vector<double> &cut_ix_probs, const arma::mat &cut_A, const arma::uvec &mst_index, const int &n)
{
  arma::mat tmp_cut_A = cut_A;
  cut_ix_probs.clear();
  double tmp_sum = 0.0;
  for(int cut_ix = 0; cut_ix < n-1; cut_ix++){
    tmp_cut_A = cut_A;
    arma::uvec cut_sub = arma::ind2sub(arma::size(tmp_cut_A), mst_index(cut_ix));
    if(tmp_cut_A(cut_sub(0), cut_sub(1)) == 0){
      Rcpp::stop("[get_min_component_size]: trying to remove an edge that doesn't exist...");
    } else{
      tmp_cut_A(cut_sub(0), cut_sub(1)) = 0; // delete the edge!
      arma::mat part_mst_A = arma::symmatl(tmp_cut_A); // adjacency matrix of the partitioned tree
      std::vector<std::vector<int> > cut_components;
      find_components(cut_components, part_mst_A);
      cut_ix_probs.push_back( (double) std::min(cut_components[0].size(), cut_components[1].size()));
      tmp_sum += (double) std::min(cut_components[0].size(), cut_components[1].size());
    }
  }
  
  for(int cut_ix = 0; cut_ix < n-1; cut_ix++) cut_ix_probs[cut_ix] /= tmp_sum;
}


// adj_support is only for the lower triangle of the adjacency matrix
void graph_partition(std::set<int> &vals, std::set<int> &l_vals, std::set<int> &r_vals, std::vector<unsigned int> &adj_support, int &K, bool &reweight, RNG &gen)
{
  arma::mat W = arma::zeros<arma::mat>(K,K); // we need to recreate the a weighted version of adjacency matrix
  // note that W is lower triangular
  for(std::vector<unsigned int>::iterator w_it = adj_support.begin(); w_it != adj_support.end(); ++w_it) W(*w_it) = gen.uniform();
  W = arma::symmatl(W); // make W symmetric
  
  // we need to subset  W to just the rows & columns corresponding to vals
  int n = vals.size();
  if(n == 1) Rcpp::stop("[graph_partition]: vals contains only one 1 value; cannot partition it further!");
  
  std::vector<double> cut_ix_probs(n-1, 1.0/( (double) (n-1) ));
  //int cut_ix = gen.multinomial(n-1,cut_ix_probs);
  int cut_ix = 0;
  std::vector<unsigned int> tmp_in;
  for(set_it it = vals.begin(); it != vals.end(); ++it) tmp_in.push_back( (unsigned int) *it);
  arma::uvec index(tmp_in); //
  arma::mat tmp_W = W(index,index);
  
  std::vector<std::vector<int> > components;
  find_components(components, tmp_W); // how many connected components are in W?
  if(components.size() != 1){
    Rcpp::stop("subgraph induced by current set of levels is not connected");
  } else{
    // now that we know levels induces a connected subgraph of the original graph
    // we can run Boruvka's algorithm
    arma::mat A_mst = boruvka(tmp_W);
    arma::mat cut_A = arma::trimatl(A_mst); // get the lower triangle of the adjacency matrix of the MST
    arma::uvec mst_index = arma::find(cut_A); // get the indices of the edges in the MST
    
    if(reweight){
      // probability that edge gets deleted is proportional to the size of the smallest component
      // that results when edge gets deleted.
      // this will lower the chance that we propose a split
      // to avoid splits that create singleton partition cells, what if we compute the size of the smallest component formed by deleting an edge
      get_edge_probs(cut_ix_probs, cut_A, mst_index, n);
    } else{
      // delete an edge uniformly at random
      cut_ix_probs.resize(n-1, 1.0/( (double) (n-1)));
    }
    cut_ix = gen.multinomial(n-1, cut_ix_probs);
  
    arma::uvec cut_sub = arma::ind2sub(arma::size(cut_A), mst_index(cut_ix)); // get the row/column subscripts corresponding to the edge being deleted
    
    if(cut_A(cut_sub(0), cut_sub(1)) != 1){
      Rcpp::stop("Trying to delete an edge that doesn't exist...");
    }
    
    cut_A(cut_sub(0), cut_sub(1)) = 0; // delete the edge!
    arma::mat part_mst_A = arma::symmatl(cut_A); // adjacency matrix of the partitioned tree
    std::vector<std::vector<int> > cut_components;
    find_components(cut_components, part_mst_A); // get the connected components of the partitioned tree!
    
    if(cut_components.size() != 2){
      Rcpp::stop("we have more than 2 connected components after deleting a single edge from a MST...");
    } else{
      l_vals.clear();
      r_vals.clear();
      
      // *it is an integer between 0 & n and indexes the temporary edge labels
      // the i-th element of tmp_in corresponds to the i-th edge label
      // we need to look up the correct level (i.e. value in vals) using tmp_in
      for(int_it it = cut_components[0].begin(); it != cut_components[0].end(); ++it) l_vals.insert( (int) tmp_in[*it]);
      for(int_it it = cut_components[1].begin(); it != cut_components[1].end(); ++it) r_vals.insert( (int) tmp_in[*it]);
    } // closes if/else checking that we have 2 components after deleting an edge
  } // closes if/else checking that levels in vals induced a connected subgraph of the original adjancecy matrix
  
}


/* eventually we will add this functionality back in.
void update_theta_cont(std::vector<double> &theta_cont, std::vector<int> &cont_var_count, int &cont_rule_count, double &a_cont, double &b_cont, int &p_cont, RNG &gen)
{
  if(theta_cont.size() != p_cont) Rcpp::stop("[update_theta_cont]: theta_cont must have size p_cont");
  double a_post = a_cont;
  double b_post = b_cont;
  for(int j = 0; j < p_cont; j++){
    a_post = a_cont + (double) cont_var_count[j];
    b_post = b_cont + (double)(cont_rule_count - cont_var_count[j]);
    theta_cont[j] = gen.beta(a_post, b_post);
  }
  
}

void update_theta_u_cat(std::vector<double> &theta_cat, std::vector<int> &cat_var_count, double &u_cat, double& a_cat, double& b_cat, int &p_cat, RNG &gen)
{
  // stuff for updating theta
  double tmp_sum = 0.0;
  int v_count = 0;
  double tmp_concentration = 0.0;
  std::vector<double> tmp_gamma(p_cat);
  
  // stuff for updating u
  double u_prop = 0.0;
  double u_orig = 0.0;
  double sum_log_theta = 0.0;
  double log_like_prop = 0.0;
  double log_like_orig = 0.0;
  double log_accept = 0.0;
  
  u_orig = u_cat;
  for(int j = 0; j < p_cat; j++){
    v_count = cat_var_count[j];
    tmp_concentration = u_orig/(1.0 - u_orig) + (double) v_count;
    tmp_gamma[j] = gen.gamma(tmp_concentration, 1.0);
    tmp_sum += tmp_gamma[j];
  }
  for(int j = 0; j < p_cat; j++){
    theta_cat[j] = tmp_gamma[j]/tmp_sum;
    sum_log_theta += log(theta_cat[j]);
  }
  
  u_prop = gen.beta(a_cat, b_cat);
  log_like_prop = (u_prop)/(1.0 - u_prop) * sum_log_theta;
  log_like_prop += lgamma( ((double) p_cat) * u_prop/(1.0 - u_prop)) - ((double) p_cat) * lgamma(u_prop/(1.0 - u_prop));
  
  log_like_orig = (u_orig)/(1.0 - u_orig) * sum_log_theta;
  log_like_orig += lgamma( ((double) p_cat) * u_orig/(1.0 - u_orig)) - ( (double) p_cat) * lgamma(u_orig/(1.0 - u_orig));
  
  log_accept = log_like_prop - log_like_orig;
  if(log_accept >= 0.0) log_accept = 0.0;
  if(gen.log_uniform() <= log_accept) u_cat = u_prop;
  else u_cat = u_orig;
}
*/


