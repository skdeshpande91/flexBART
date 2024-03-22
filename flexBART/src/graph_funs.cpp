#include "graph_funs.h"

// our graph algorithms defined here

// sometimes it's convenient to have a redundant (i.e. symmetric) representation of the edges
// each edge represented in two elements of the map with keys correspond to the "source" and "sink"
// will call this specifically in boruvka's and in the LERW style code
void parse_edge_mat(std::vector<edge> &edges, Rcpp::NumericMatrix &edge_mat, int &n_vertex)
{
  int n_edges = edge_mat.rows();
  edges.clear();
  if(edge_mat.cols() == 3){
    for(int i = 0; i < n_edges; i++){
      edges.push_back(edge( (int) edge_mat(i,0), (int) edge_mat(i,1), edge_mat(i,2)));
    }
  } else if(edge_mat.cols() == 2){
    for(int i = 0; i < n_edges; i++){
      edges.push_back(edge( (int) edge_mat(i,0), (int) edge_mat(i,1), 1.0));
    }
  } else{
    Rcpp::stop("[parse_edge_mat]: The matrix edge_mat must have 2 columns (unweighted graph) or 3 columns (weighted graph)");
  }
}

// takes in the List of edge_mat's
void parse_graphs(std::vector<std::vector<edge>> &edges, int &p_cat, std::vector<int> &K, Rcpp::List &tmp_edge_mats, Rcpp::LogicalVector &graph_split)
{
  edges.clear();
  edges.resize(p_cat, std::vector<edge>());
  if(tmp_edge_mats.size() == p_cat){
    for(int j = 0; j < p_cat; j++){
      if(graph_split(j) == 1){
        Rcpp::NumericMatrix edge_mat = Rcpp::as<Rcpp::NumericMatrix>(tmp_edge_mats[j]);
        parse_edge_mat(edges[j], edge_mat, K[j]);
      } else{
        // do nothing
      }
    }
  } else{
    Rcpp::Rcout << "[parse_graphs]: detected " << p_cat << " categorical variables";
    Rcpp::Rcout << " edge_mat_list has length " << tmp_edge_mats.size() << std::endl;
    Rcpp::stop("edge_mat_list must have length equal to p_cat!");
  }
}



void build_symmetric_edge_map(edge_map &emap, std::vector<edge> &edges, std::set<int> &vertices)
{
  emap.clear();
  
  // create an element with a key for each vertex in the specified set
  for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it) emap.insert(std::pair<int, std::vector<edge>>(*v_it, std::vector<edge>()));
  
  for(edge_vec_it it = edges.begin(); it != edges.end(); ++it){
    
    if(emap.count(it->source) != 1 || emap.count(it->sink) != 1){
      Rcpp::Rcout << "[build_symmetric_edge_map]: Supplied set includes an edge with one vertex outside set of supplied vertices!" << std::endl;
      Rcpp::Rcout << "  edge is from " << it->source << " to " << it->sink << std::endl;
      Rcpp::Rcout << "  Supplied vertices:";
      for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it) Rcpp::Rcout << " " << *v_it;
      Rcpp::Rcout << std::endl;
    } else{
      emap.find(it->source)->second.push_back(edge(it->source, it->sink, it->weight));
      emap.find(it->sink)->second.push_back(edge(it->sink, it->source, it->weight));
    }
  } // finish looping over all of the edges
}


// in the main BART loop, we will often have a set of categories that form a subset of the vertices of a network
// we will have to partition the subgraph induced by this set of vertices
// get_induced_edges loops over all edges and pulls out those whose source & sink are in the vertex_subset
// we need to pull out the edges
std::vector<edge> get_induced_edges(std::vector<edge> &edges, std::set<int> &vertex_subset)
{
  std::vector<edge> subset_edges;
  // loop over the edges
  for(edge_vec_it it = edges.begin(); it != edges.end(); ++it){
    // check that both the source and the sink belong to the vertex_set
    if(vertex_subset.count(it->source) == 1 && vertex_subset.count(it->sink) == 1){
      subset_edges.push_back(edge(it->source, it->sink, it->weight));
    }
  }
  return subset_edges;
}

// stuff for finding connected components:
//  depth-first search: we maintain a map with keys corresponding to vertices and boolean values; this records which vertices our DFS has visited
//  dfs should only ever be called using a *symmetric* edge_map; we need a key for every vertex
void dfs(int v, std::map<int, bool> &visited, std::vector<int> &comp, edge_map &emap)
{
  // this is the first time we have visited vertex labelled v
  std::map<int,bool>::iterator v_it = visited.find(v);
  v_it->second = true; // mark that we have visited vertex v
  comp.push_back(v); // now that v has been marked, we can add it to the connected component

  // now find the element of our map containing all verties whose source is v
  if(emap.count(v) != 1){
    // if we are calling dfs correctly this should *never* be hit
    Rcpp::Rcout << "[dfs]: at vertex" << v << std::endl;
    Rcpp::Rcout << "keys in emap are:";
    for(edge_map_it tmp_it = emap.begin(); tmp_it != emap.end(); ++tmp_it) Rcpp::Rcout << " " << tmp_it->first << std::endl;
    Rcpp::stop("something is wrong with our edge map");
  } else{
    edge_map_it v_edges_it = emap.find(v);
    for(edge_vec_it it = v_edges_it->second.begin(); it != v_edges_it->second.end(); ++it){
      // it points to an edge leaving v, let us see if the sink has been visited
      int vv = it->sink;
      std::map<int,bool>::iterator vv_it = visited.find(vv);
      if(vv_it == visited.end()){
        // we were unable to find an entry for vv in visited.
        // if we are calling dfs correctly, this should *never* be hit
        Rcpp::Rcout<< "[dfs]: Traversing edges from " << v << " to " << vv << std::endl;
        Rcpp::Rcout << "  did not find entry for " << vv << " in visited" << std::endl;
        Rcpp::Rcout << "  keys of visited:";
        for(std::map<int,bool>::iterator visit_it = visited.begin(); visit_it != visited.end(); ++visit_it) Rcpp::Rcout << " " << visit_it->first;
        Rcpp::Rcout << std::endl;
        Rcpp::stop("something is wrong with our edge map or visited!");
      } else{
        if(!vv_it->second){
          // we have not yet visited the vertex vv and so we continue our dfs from there
          dfs(vv, visited, comp, emap);
        }
      }
    } // closes loop over edges incident to v
  } // closes if/else checking that we have an entry for v in our edge_map
}


void find_components(std::vector<std::vector<int> > &components, std::vector<edge> &edges, std::set<int> &vertices)
{
  components.clear();
  
  if(edges.size() == 0){
    // no edges: every vertex is is its own component
    //Rcpp::Rcout << "[find_components]: no edges" << std::endl;
    for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it) components.push_back(std::vector<int>(1, *v_it));
  } else{
    std::map<int, bool> visited;
    // initially mark every vertex as unvisited
    for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it) visited.insert(std::pair<int,bool>(*v_it, false));
    edge_map emap;
    build_symmetric_edge_map(emap, edges, vertices);
    
    // loop over the vertices and run dfs if we haven't already visited that vertex
    for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it){
      //Rcpp::Rcout << " starting from " << *v_it << std::endl;
      if(!visited.find(*v_it)->second){
        // we have not yet visited vertex labelled *v_it, so it represents a new component
        std::vector<int> new_comp;
        dfs(*v_it, visited, new_comp, emap);
        // by the time the dfs finishes, we have visited everything in the current connected component
        components.push_back(new_comp); // add current connected component to our vector of components
      }
    }
  }
}

void get_unique_edges(std::vector<edge> &edges)
{
  std::vector<edge> unik_edges;
  if(edges.size() > 0){
    unik_edges.push_back(edges[0]);
    for(std::vector<edge>::iterator it = edges.begin(); it != edges.end(); ++it){
      bool add_edge = true;
      for(std::vector<edge>::iterator uit = unik_edges.begin(); uit != unik_edges.end(); ++uit){
        if( (it->source == uit->source && it->sink == uit->sink) || (it->sink == uit->source && it->source == uit->sink)){
          // *it (or its reverse) is already in our set of unique edges. DO NOT keep comparing *it to elements of unik_edges & DO NOT add it to unik_edges
          //Rcpp::Rcout << "duplicate edge found!" << std::endl;
          add_edge = false;
          break;
        }
      }
      if(add_edge) unik_edges.push_back(*it);
    }
  }
  edges.clear();
  for(std::vector<edge>::iterator it = unik_edges.begin(); it != unik_edges.end(); ++it) edges.push_back(edge(it->source, it->sink, it->weight));
}

arma::mat get_adjacency_matrix(edge_map emap){
  // adjacency matrix uses the indices 0 --> n_vertex-1
  // we need a mapping from our actual vertices (given by the keys of emap)
  // to the indices of the adjacency matrix
  
  int n_vertex = emap.size();
  std::map<int, int> vertex_index_map;
  int v_counter = 0;
  for(std::map<int, std::vector<edge>>::iterator v_it = emap.begin(); v_it != emap.end(); ++v_it){
    vertex_index_map.insert(std::pair<int,int>(v_it->first, v_counter));
    ++v_counter;
  }
  arma::mat A = arma::zeros<arma::mat>(n_vertex, n_vertex);
  
  int source_index;
  int sink_index;
  
  for(std::map<int, std::vector<edge>>::iterator v_it = emap.begin(); v_it != emap.end(); ++v_it){
    source_index = vertex_index_map.find(v_it->first)->second;
    for(std::vector<edge>::iterator e_it = v_it->second.begin(); e_it != v_it->second.end(); ++e_it){
      if(source_index != vertex_index_map.find(e_it->source)->second){
        Rcpp::stop("Something wrong with edge indexing!");
      }
      sink_index = vertex_index_map.find(e_it->sink)->second;
      A(source_index, sink_index) = e_it->weight;
    }
  }
  return(A);
}

arma::mat floydwarshall(std::vector<edge> &edges, std::set<int> &vertices)
{
  // get symmetric edge map
  edge_map emap;
  build_symmetric_edge_map(emap, edges, vertices);
  
  // at this point indices run from 0 to (n_vertex-1)
  // we need lookup tables
  std::map<int,int> vertex_index_map;
  std::map<int,int> index_vertex_map;
  int v_counter = 0;
  for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it){
    vertex_index_map.insert(std::pair<int,int>(*v_it, v_counter));
    index_vertex_map.insert(std::pair<int,int>(v_counter, *v_it));
    ++v_counter;
  }
  
  int n_vertex = vertices.size(); // number of vertices
  //arma::mat A = get_adjacency_matrix(emap);
  
  
  
  arma::mat D = get_adjacency_matrix(emap);
  D.elem(arma::find(D == 0)).fill(pow(n_vertex, 2.0));
  D.diag().zeros();
  
  for(int k = 0; k < n_vertex; ++k){
    //Rcpp::Rcout << "Starting k = " << k << std::endl;
    for(int i = 0; i < n_vertex; ++i){
      for(int j = 0; j < n_vertex; ++j){
        //Rcpp::Rcout << "i = " << i << "j = " << j << "D(i,j) = " << D(i,j);
        //Rcpp::Rcout << " val = " << D(i,k) + D(k,j) << std::endl;
        if(D(i,j) > D(i,k) + D(k,j)){
          //Rcpp::Rcout << "Found shorter path!" << std::endl;
          D(i,j) = D(i,k) + D(k,j);
        }
        
      }
      //Rcpp::Rcout << std::endl;
    }
  }
  return(D);
  
}

void boruvka(std::vector<edge> &mst_edges, std::vector<edge> &edges, std::set<int> &vertices)
{
  
  mst_edges.clear();
  
  edge_map emap;
  build_symmetric_edge_map(emap, edges, vertices);

  std::vector<std::vector<int> > mst_components;
  find_components(mst_components, mst_edges, vertices);
  int n_vertex = vertices.size();
  int counter = 0;
  // Boruvka has worst-case complexity of O(log(n)), so we add a considerable buffer
  while( mst_components.size() > 1 && counter < 2*n_vertex ){
    std::vector<edge> new_edges; // the new edges that will get added to the MST
    //Rcpp::Rcout << " Starting Round " << counter << " of Boruvka. Edges are:" << std::endl;
    //for(std::vector<edge>::iterator mst_eit = mst_edges.begin(); mst_eit != mst_edges.end(); ++mst_eit){
    //  Rcpp::Rcout << mst_eit->source << " to " << mst_eit->sink << std::endl;
    //}
    for(std::vector<std::vector<int> >::iterator comp_it = mst_components.begin(); comp_it != mst_components.end(); ++comp_it){
      // looping over all of the existing components in the current MST
      int min_source = 0;
      int min_sink = 0;
      double min_weight = 1.0;
      for(std::vector<int>::iterator v_it = comp_it->begin(); v_it != comp_it->end(); ++v_it){
        // v_it points to a particular vertex in the component *comp_it.
        // We need to loop over all incident edges to *v_it
        // we check whether (A) the "sink" is outside the component and (B) whether it has smallest weight
        
        //Rcpp::Rcout << "    Visiting vertex " << *v_it << std::endl;
        
        edge_map_it ve_it = emap.find(*v_it); // ve_it points to the vector of edges incident to *v_it
        for(std::vector<edge>::iterator e_it = ve_it->second.begin(); e_it != ve_it->second.end(); ++e_it){
          // e_it points to a specific edge
          // we check whether (A) the "sink" is outside the component and (B) whether it has smallest weight
          //Rcpp::Rcout << "   checking edge from " << e_it->source << " to " << e_it->sink << " count = " << std::count(comp_it->begin(), comp_it->end(), e_it->sink) << std::endl;
          if(std::count(comp_it->begin(), comp_it->end(), e_it->sink) == 0 && e_it->weight < min_weight){
            // found a new minimum edge weight leaving the component!
            min_source = e_it->source; // this had better be *v_it
            min_sink = e_it->sink;
            min_weight = e_it->weight;
          }
        } // closes loop over edge incident to vertex *v_it
      } // closes loop over vertices in component *comp_it
      new_edges.push_back(edge(min_source, min_sink, min_weight));
    } // closes loop over components
    
    // at this point, we've finished looping over all of the vertices in a particular component and we have found the edge
    // that leaves the component and has minimum weight. we dump that edge into our running collection of edges in the MSt
    
    for(std::vector<edge>::iterator it = new_edges.begin(); it != new_edges.end(); ++it){
      mst_edges.push_back(*it);
    }
    get_unique_edges(mst_edges); // we have may duplicate edges so kill them off here.
    find_components(mst_components, mst_edges, vertices); // re-compute the number of components
    counter++;
  } // closes main while loop
  
  if(mst_components.size() > 1){
    Rcpp::Rcout << "[boruvka]: after " << counter << " rounds, we have not yet formed a single connected MST" << std::endl;
    Rcpp::Rcout << "  returning an empty set of edges" << std::endl;
    mst_edges.clear();
  }
}


void wilson(std::vector<edge> &mst_edges, std::vector<edge> &edges, std::set<int> &vertices, RNG &gen)
{
  int n_vertex = vertices.size();
  edge_map emap;
  build_symmetric_edge_map(emap, edges, vertices);
  
  
  std::map<int, bool> in_tree; // key is vertex label, value is boolean of whether vertex is in tree
  std::vector<int> possible_roots; // we have to start the tree from somewhere

  for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it){
    in_tree.insert(std::pair<int, bool>(*v_it, false)); // each vertex is not tree
    possible_roots.push_back(*v_it);
  }
  int root = possible_roots[floor(gen.uniform() * n_vertex)];
  in_tree.find(root)->second = true; // mark the root as a member of the tree
  
  
  bool tree_full = true; // true when all vertices are members of the tree
  
  // when we start, we know most values in in_tree are false so this loop isn't super expensive
  for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it){
    if(!in_tree.find(*v_it)->second){
      tree_full = false;
      break;
    }
  }
  
  possible_roots.clear();
  // from now on, possible_roots will store the potential roots for the LERW to the tree
  int outer_counter = 0; // number of rounds in which we try to walk to the tree
  
  while(!tree_full && outer_counter <= n_vertex){
    // each time through the loop, we choose a new vertex and walk from it to the tree
    // so long as the inner loop terminates (i.e. we successfully walk to the tree)
    // we should never have to run the outer loop more than n_vertex times
    
    possible_roots.clear();
    for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it){
      if(!in_tree.find(*v_it)->second) possible_roots.push_back(*v_it);
    }
    int next_state = possible_roots[floor(gen.uniform() * possible_roots.size())]; // next_state is start of the walk
    int old_state;
    
    std::vector<int> lew;
    std::vector<edge> lew_edges;
    // as we do the loop erasure, we have to keep track of which vertices have been visited
    // key is vertex label and value is boolean for whether we visited or not
    std::map<int, bool> visited;
    for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it){
      visited.insert(std::pair<int,bool>(*v_it, false));
    }
    
    lew.push_back(next_state); // only used for testing
    visited.find(next_state)->second = true; // mark the root as having been visited
    
    
    int inner_counter = 0;
    int edge_index;
    
    // next few lines only for testing purposes
    //Rcpp::Rcout << "Starting Round " << outer_counter << "! Tree contains:";
    //for(std::map<int,bool>::iterator tr_it = in_tree.begin(); tr_it != in_tree.end(); ++tr_it){
    //  if(tr_it->second) Rcpp::Rcout << " " << tr_it->first;
    //}
    //Rcpp::Rcout << std::endl;
    // done with the printing used to test
    
    while(!in_tree.find(next_state)->second && inner_counter < (int) 10 * pow(n_vertex,3.0)){
      old_state = next_state;
      edge_map_it em_it = emap.find(old_state); // em_it now points to vector of edges incident to old_state
      if(em_it->second.size() == 0){
        //Rcpp::Rcout << "Random walk is a vertex that has no incident edges. Stopping now" << std::endl;
        mst_edges.clear(); // clear out mst_edges
        tree_full = true;
        break; // break out of inner loop
      } else{
        std::vector<double> weights;
        for(std::vector<edge>::iterator e_it = em_it->second.begin(); e_it != em_it->second.end(); ++e_it){
          weights.push_back(e_it->weight);
        }
        edge_index = gen.categorical(weights);
        
        next_state = em_it->second[edge_index].sink;
        if(!visited.find(next_state)->second){
          // we have not yet visited the vertex next_state
          lew.push_back(next_state); // only for testing purposes
          lew_edges.push_back( em_it->second[edge_index] ); // add the edge to next_state
          visited.find(next_state)->second = true; // mark next_state has having been visited
        } else{
          // we have already visited next_state
          std::vector<int>::iterator prev_it = lew.end()-1; // prev_it points to the last current state in the LERW
          std::vector<edge>::iterator prev_eit = lew_edges.end() - 1; // points to last edge in the LERW
          for(;prev_it != lew.begin(); --prev_it){
            if(*prev_it == next_state) break;
          }
          // prev_it should point to the last instance of next_state in lew
          for(;prev_eit != lew_edges.begin(); --prev_eit){
            if(prev_eit->sink == next_state) break;
          }
          // prev_eit should now point to the last edge in lew that goes into next_state
          
          if(prev_it == lew.begin() && lew[0] != next_state){
            // if we've reached the beginning of lew but haven't seen next_state something has gone quite wrong
            Rcpp::Rcout << "[wilson]: LERW trying to re-visit " << next_state << " but could not find last visit" << std::endl;
            Rcpp::Rcout << "  LERW currently contains:";
            for(std::vector<int>::iterator lew_it = lew.begin(); lew_it != lew.end(); ++lew_it) Rcpp::Rcout << " " << *lew_it;
            Rcpp::Rcout << std::endl;
            Rcpp::stop("could not perform loop-erasure");
            
          } else if(prev_it == lew.begin() && lew[0] == next_state){
            // we are attempting to revisit the the starting point of the LERW
            // we reset lew to contain just next_state and empty lew_edges
            lew.clear();
            lew.push_back(next_state);
            lew_edges.clear();
            
            // we need to update visited: in this case, only next_state should be marked as visited
            for(std::map<int, bool>::iterator visit_it = visited.begin(); visit_it != visited.end(); ++visit_it){
              if(visit_it->first != next_state) visit_it->second = false;
              else visit_it->second = true;
            }
          } else{
            // we have found the last time the LERW has visited next state and it wasn't at the beginning
            lew.erase(prev_it, lew.end()); // erase also kills the last visit to next_state...
            lew.push_back(next_state); // ...so we add that back in!
            
            // prev_eit points to the last edge that goes INTO next_state
            // this is not an edge in the cycle being popped
            int prev_source = prev_eit->source;
            int prev_sink = prev_eit->sink; // had better be equal to next_state!
            double prev_weight = prev_eit->weight;
            
            lew_edges.erase(prev_eit, lew_edges.end()); // not really sure what prev_eit will point to after the erase...
            lew_edges.push_back(edge(prev_source, prev_sink, prev_weight)); // ... so we make a new copy of the edge to be safe
            
            // we must update visited.
            // looping over lew is not sufficient because we may have removed a visited vertex while popping a cycle
            for(std::map<int, bool>::iterator visit_it = visited.begin(); visit_it != visited.end(); ++visit_it){
              if(std::count(lew.begin(), lew.end(), visit_it->first) == 1) visit_it->second = true;
              else if(std::count(lew.begin(), lew.end(), visit_it->first) == 0) visit_it->second = false;
              else{
                Rcpp::Rcout << "[wilson]: After popping cycle, our LERW visits a vertex twice!" << std::endl;
                Rcpp::Rcout << "  LERW currently contains:";
                for(std::vector<int>::iterator lew_it = lew.begin(); lew_it != lew.end(); ++lew_it) Rcpp::Rcout << " " << *lew_it;
                Rcpp::Rcout << std::endl;
                Rcpp::stop("LERW somehome contains a loop...");
              }
            } // closes loop updating visited after cycle popping
          } // closes if/else checking when we visited next_state
        } // closes if/else checking whether we have visited next_state (and popping cycle if so)
      } // closes if/else checking that old_state has incident edges
      inner_counter++;
    } // closes inner while loop
    
    if(inner_counter == (int) 10*pow(n_vertex,3.0)){
      Rcpp::Rcout << "[wilson]: inner loop did not hit tree in " << 10*pow(n_vertex,3.0) <<  " steps. stopping now" << std::endl;
      Rcpp::Rcout << "  LERW currently contains :";
      for(std::vector<int>::iterator lew_it = lew.begin(); lew_it != lew.end(); ++lew_it) Rcpp::Rcout << " " << *lew_it;
      Rcpp::Rcout << std::endl;
      Rcpp::stop("may need to increase number of steps allowed for inner loop!");
    } else{
      // dump edges from lew into mst_edges
      for(std::vector<edge>::iterator e_it = lew_edges.begin(); e_it != lew_edges.end(); ++e_it) mst_edges.push_back(*e_it);
      // update in_tree
      for(std::vector<int>::iterator lew_it = lew.begin(); lew_it != lew.end(); ++lew_it) in_tree.find(*lew_it)->second = true;
      
      tree_full = true;
      for(std::map<int, bool>::iterator tr_it = in_tree.begin(); tr_it != in_tree.end(); ++tr_it){
        if(!tr_it->second){
          tree_full = false;
          break;
        }
      }
      
      // some printing only used for testing
      //Rcpp::Rcout << "Ending Round " << outer_counter << "! Tree contains:";
      //for(std::map<int,bool>::iterator tr_it = in_tree.begin(); tr_it != in_tree.end(); ++tr_it){
      //  if(tr_it->second) Rcpp::Rcout << " " << tr_it->first;
      //}
      //Rcpp::Rcout << std::endl;
    }
    outer_counter++;
  } // closes outer loop
  
  if(!tree_full){
    Rcpp::Rcout << "[wilson]: unable to reach every vertex in the graph. returning an empty set of edges" << std::endl;
    mst_edges.clear();
  }
}

void delete_unif_edge(std::set<int> &l_vals, std::set<int> &r_vals, std::vector<edge> &edges, std::set<int> &vertices, RNG &gen){
  
  //std::vector<double> cut_probs(;
  //for(std::vector<edge>::iterator e_it = edges.begin(); e_it != edges.end(); ++e_it){
  //  cut_probs.push_back(e_it->weight);
  //}
  //int cut_index = gen.categorical(cut_probs);
  int cut_index = floor(gen.uniform() * edges.size()); // pick index of edge to be deleted uniformly
  
  std::vector<edge> cut_edges;
  for(int e = 0; e < edges.size(); e++){
    if(e != cut_index) cut_edges.push_back(edges[e]);
  }
  std::vector<std::vector<int> > cut_components;
  find_components(cut_components, cut_edges, vertices);
  if(cut_components.size() != 2){
    Rcpp::Rcout << "[delete_unif_edge]: Attempted to partition connected subgraph by deleting an edge" << std::endl;
    Rcpp::Rcout << "   Resulting graph does not have 2 components. supplied graph might not have been a tree" << std::endl;
    Rcpp::Rcout << "  Returning empty l_vals and r_vals!" << std::endl;
    //Rcpp::stop("error in cutting edge from spanning tree");
    // return empty l_vals and r_vals and handle this exception in a higher function
    l_vals.clear();
    r_vals.clear();
  } else{
    l_vals.clear();
    r_vals.clear();
    for(std::vector<int>::iterator it = cut_components[0].begin(); it != cut_components[0].end(); ++it) l_vals.insert(*it);
    for(std::vector<int>::iterator it = cut_components[1].begin(); it != cut_components[1].end(); ++it) r_vals.insert(*it);
  } // closes if/else checking that we have 2 connected components after deleting edge from spanning tree
}

void signcheck_split(std::set<int> &l_vals, std::set<int> &r_vals,
                std::vector<edge> &edges, std::set<int> &vertices)
{
  // first we need our edge map
  edge_map emap;
  build_symmetric_edge_map(emap, edges, vertices);
  
  std::map<int,int> vertex_index_map;
  std::map<int,int> index_vertex_map;
  int v_counter = 0;
  for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it){
    vertex_index_map.insert(std::pair<int,int>(*v_it, v_counter));
    index_vertex_map.insert(std::pair<int,int>(v_counter, *v_it));
    ++v_counter;
  }
  
  int n_vertex = vertices.size();
  
  // get adjacency matrix
  arma::mat A = get_adjacency_matrix(emap);
  arma::mat D = arma::zeros<arma::mat>(n_vertex, n_vertex);
  for(int i = 0; i < n_vertex; ++i) D(i,i) = arma::accu(A.row(i));
  arma::mat L = D - A;
  
  arma::vec eigval;
  arma::mat eigvec;
  bool eigen = arma::eig_sym(eigval, eigvec, L);
    
  if(eigen){
    arma::uvec lval_index = arma::find(eigvec.col(1) < 0);
    l_vals.clear();
    for(int ix = 0; ix < lval_index.size(); ++ix){
      l_vals.insert(index_vertex_map.find(lval_index(ix))->second);
    }
    for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it){
      if(l_vals.count(*v_it) == 0) r_vals.insert(*v_it);
    }
   
  } else{
    Rcpp::Rcout << "[signcheck_split]: Eigendecomposition failed. Returning empty l_vals and r_vals" << std::endl;
    l_vals.clear();
    r_vals.clear();
  }
  
}
/*
 cut_type = 0: randomly pick one of the processes
 cut_type = 1: wilson + delete uniform edge
 cut_type = 2: wilson + partition based on sign of 2nd eigenvector of Laplacian
 cut_type = 3: just do cut based on sign of 2nd eigenvector Laplacian
 cut_type = 4: hotspot
 */
void graph_partition(std::set<int> &l_vals, std::set<int> &r_vals, std::vector<edge> &orig_edges,std::set<int> &avail_levels, int &cut_type, RNG &gen)
{
  // get the edges for the induced subgraph
  std::vector<edge> edges = get_induced_edges(orig_edges, avail_levels);
  
  // check if the induced subgraph is connected
  std::vector<std::vector<int> > components;
  find_components(components, edges, avail_levels);
  if(components.size() == 0){
    // this should never be reached
    Rcpp::stop("[graph_partition]: graph has no components...");
  } else if(components.size() > 1){
    // graph is not connected. assign whole components to left and right children uniformly at random
    l_vals.clear();
    r_vals.clear();
    int rule_counter = 0;
    while( ((l_vals.size() == 0 || r_vals.size() == 0)) && rule_counter < 1000 ){
      l_vals.clear();
      r_vals.clear();
      for(int comp_ix = 0; comp_ix < components.size(); comp_ix++){
        if(gen.uniform() <= 0.5){
          // send everything in this component to the left child
          for(int_it it = components[comp_ix].begin(); it != components[comp_ix].end(); ++it) l_vals.insert(*it);
        } else{
          // send everything in this component to the right child
          for(int_it it = components[comp_ix].begin(); it != components[comp_ix].end(); ++it) r_vals.insert(*it);
        }
      }
      ++rule_counter;
    }
    if(rule_counter == 1000){
      Rcpp::stop("[graph partition]: graph disconnected. failed to generate a non-trivial partiton of components in 1000 attempts!");
    }
  } else{
    // induced subgraph is connected and we try to partition it
    int tmp_cut_type = cut_type;
    while(tmp_cut_type == 0){
      tmp_cut_type = floor(gen.uniform() * 4.0) + 1; // uniform on the set 1, 2, 3, 4
    }
    if(tmp_cut_type == 1 || tmp_cut_type == 2){
      std::vector<edge> ust_edges;
      wilson(ust_edges, edges, avail_levels, gen);
      if(tmp_cut_type == 1) delete_unif_edge(l_vals, r_vals, ust_edges, avail_levels, gen);
      else signcheck_split(l_vals, r_vals, ust_edges, avail_levels);
    } else if(tmp_cut_type == 3){
      signcheck_split(l_vals, r_vals, edges, avail_levels);
    } else if(tmp_cut_type == 4){
      hotspot(l_vals, r_vals, edges, avail_levels, gen);
    } else{
      // should never hit this
      Rcpp::Rcout << "[graph_partition]: tmp_cut_type = " << tmp_cut_type << std::endl;
      Rcpp::stop("tmp_cut_type should never exceed 4!");
    } // closes if/else checking tmp_cut_type
    
    if(l_vals.size() == 0 || r_vals.size() == 0){
      Rcpp::Rcout << "[graph_partition]: one of l_vals or r_vals contains no elements. something is wrong" << std::endl;
      Rcpp::stop("could not generate nontrivial graph partition!");
    }
    
  } // closes if/else checking whether the graph induced by avail_levels is connected or not
}

void hotspot(std::set<int> &l_vals, std::set<int> &r_vals, std::vector<edge> &edges, std::set<int> &vertices, RNG &gen, bool debug)
{
  // first we need our edge map
  edge_map emap;
  build_symmetric_edge_map(emap, edges, vertices);
  
  
  //std::map<int,int> vertex_index_map;
  std::map<int,int> index_vertex_map;
  int v_counter = 0;
  for(std::map<int,std::vector<edge>>::iterator v_it = emap.begin(); v_it != emap.end(); ++v_it){
    //vertex_index_map.insert(std::pair<int,int>(v_it->first, v_counter));
    index_vertex_map.insert(std::pair<int,int>(v_counter, v_it->first));
    ++v_counter;
  }
  
  int n_vertex = vertices.size();
  arma::mat D = floydwarshall(edges, vertices);
  
  int seed_index = floor(gen.uniform() * n_vertex);
  int max_dist = pow(n_vertex, 2);
  int radius = 1;
  std::vector<double> radius_probs;
  
  int attempt = 0;
  int max_attempt = 10;
  bool trivial_split = true;
  
  while( (attempt < max_attempt) && trivial_split){
    if(debug) Rcpp::Rcout << "[hotspot]: Starting attempt = " << attempt << std::endl;
    seed_index = floor(gen.uniform() * n_vertex); // seed of the hotspot
    if(debug) Rcpp::Rcout << "  seed = " << seed_index << std::endl;
    max_dist = D.col(seed_index).max(); // maximum distance in graph from hotspot seed
    
    if(max_dist == 1){
      // everything is at distance 1 from hotspot seed
      // we set radius = 1 so that l_vals contains all vertices & r_vals contains none
      // note: this is a trivial split.
      trivial_split = true;
      radius = 1;
    } else if(max_dist == 2){
      // there is a vertex at distance 2 from hotspot seed
      // we will set radius = 1 so that
      //   l_vals contains seed + its neighbors
      //   r_vals contains at least one other vertex
      // this results in a non-trivial splitting rule so we can set trivial_split = false
      trivial_split = false;
      radius = 1;
    } else{
      // we have options for the radius.
      // P(radius = r) = 0.5^r for r = 1, 2, ..., max_dist - 2; residual prob for r = max_dist - 1
      trivial_split = false;
      radius_probs.clear();
      radius_probs.resize(max_dist-1);
      radius_probs[max_dist-2] = 1.0;
      for(int r = 0; r < max_dist-2; ++r){
        radius_probs[r] = pow(0.5,r+1);
        radius_probs[max_dist-2] -= pow(0.5, r+1);
      }
      if(debug){
        if(debug) Rcpp::Rcout << "  radius probabilities:";
        for(int r = 0; r < max_dist-1; ++r) Rcpp::Rcout << " " << radius_probs[r];
        Rcpp::Rcout << std::endl;
      }
      radius = gen.categorical(radius_probs) + 1; // without +1, we will get radius of 0 sometimes...
    }
    if(debug) Rcpp::Rcout << "  radius = " << radius << std::endl;
    
    l_vals.clear();
    r_vals.clear();
    
    arma::uvec l_index = arma::find(D.col(seed_index) <= radius);
    for(int ix = 0; ix < l_index.size(); ++ix){
      l_vals.insert(index_vertex_map.find(l_index(ix))->second);
    }
    for(std::set<int>::iterator v_it = vertices.begin(); v_it != vertices.end(); ++v_it){
      if(l_vals.count(*v_it) == 0) r_vals.insert(*v_it);
    }
    
    ++attempt;
  }
}

