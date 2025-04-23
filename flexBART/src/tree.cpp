#include <string>
#include <vector>
#include <map>
#include "tree.h"

//--------------------------------------------------
// constructors

tree::tree(){
  mu = 0.0;
  
  rule.is_cat = false;
  rule.v_aa = 0;
  rule.c = 0.0;
  
  rule.v_cat = 0;
  rule.l_vals = std::set<int>();
  rule.r_vals = std::set<int>();
  
  p = 0;
  l = 0;
  r = 0;
}

// for the destructor: cut back to one node
void tree::to_null()
{
  size_t ts = get_treesize();
  while(ts > 1){
    npv nv;
    get_nogs(nv); // get nodes with no grandchildren
    for(size_t i = 0; i < nv.size(); i++){
      delete nv[i]->l;
      delete nv[i]->r;
      nv[i]->l=0;
      nv[i]->r=0;
    }
    ts = get_treesize();
  }
  mu = 0.0;
  rule.clear();
  p = 0;
  l = 0;
  r = 0;
}

// print
void tree::print(bool pc) const // pc is flag to print children
{
  size_t id = get_nid();
  
  if(pc && (get_ntype() == 't')) Rcpp::Rcout << "tree size:" << get_treesize() << std::endl;
  Rcpp::Rcout << "id." << id;
  if(get_ntype() == 'b') Rcpp::Rcout << "  mu: " << mu << std::endl;
  else if(get_ntype() == 't' && get_treesize() == 1){
    // tree is just the top node
    Rcpp::Rcout << "  mu: " << mu << std::endl;
  } else { // internal node or nnog or top node
    if(!rule.is_cat){
      Rcpp::Rcout << "  axis-aligned split on X_cont[," << rule.v_aa+1 << "] at c = " << rule.c << std::endl;
    } else{
      // categorical decision rule
      Rcpp::Rcout << "  split on categorical X_cat[," << rule.v_cat+1 << "]" << std::endl;
      Rcpp::Rcout << "    left levels: ";
      for(set_it it = rule.l_vals.begin(); it != rule.l_vals.end(); ++it) Rcpp::Rcout << " " << *it;
      Rcpp::Rcout << std::endl;
      
      Rcpp::Rcout << "    right levels: ";
      for(set_it it = rule.r_vals.begin(); it != rule.r_vals.end(); ++it) Rcpp::Rcout << " " << *it;
      Rcpp::Rcout << std::endl;
    }
  }
  
  if(pc){
    if(l){
      l->print(pc);
      r->print(pc);
    }
  }
}



//--------------------------------------------------
//operators

// this overloads = if we say tree1 = tree2

tree& tree::operator=(const tree& rhs)
{
   if(&rhs != this) {
      to_null(); //kill left hand side (this)
      cp(this,&rhs); //copy right hand side to left hand side
   }
   return *this;
}

// tree-level gets

// get a vector of pointers to the bottom nodes
void tree::get_bots(npv& bv)
{
  if(l) { //have children
    l->get_bots(bv);
    r->get_bots(bv);
  } else bv.push_back(this);
}

//get a vector of pointers to the no grandchildren nodes
void tree::get_nogs(npv& nv)
{
  if(l) { //have children
    if((l->l) || (r->l)) {  //have grandchildren
      if(l->l) l->get_nogs(nv);
      if(r->l) r->get_nogs(nv);
    } else nv.push_back(this);
  }
}

bool tree::is_nog() const
{
  bool isnog=true;
  if(l) {
    if(l->l || r->l) isnog=false; //one of the children has children.
  } else isnog=false; //no children
  return isnog;
}

//get a vector of pointers to *ALL* nodes
void tree::get_nodes(npv& v)
{
  v.push_back(this);
  if(l) {
    l->get_nodes(v);
    r->get_nodes(v);
  }
}

void tree::get_nodes(cnpv& v)  const
{
  v.push_back(this);
  if(l) {
    l->get_nodes(v);
    r->get_nodes(v);
  }
}


// get the size of the tree (i.e. number of nodes)
int tree::get_treesize() const
{
   if(!l) return 1;  //if bottom node, tree size is 1
   else return (1+l->get_treesize()+r->get_treesize());
}
//number of nodes with no
int tree::get_nnogs() const
{
  if(!l) return 0; // this is a bottom node
  if(l->l || r->l) return(l->get_nnogs() + r->get_nnogs()); // this has at least one grandchild
  else return 1; // this is a nog node
}
// count number of bottom nodes
int tree::get_nbots() const
{
  if(!l) return 1; // this is a bottom node
  else return(l->get_nbots() + r->get_nbots());
}

// get depth of tree
int tree::get_depth() const
{
   if(!p) return 0; //no parents
   else return (1+p->get_depth());
}

// get the node id
int tree::get_nid() const
//recursion up the tree
{
   if(!p) return 1; //if you don't have a parent, you are the top
   if(this==p->l) return 2*(p->get_nid()); //if you are a left child
   else return 2*(p->get_nid())+1; //else you are a right child
}

// get the node type
char tree::get_ntype() const
{
   //t:top, b:bottom, n:no grandchildren, i:internal
  if(p == 0) return 't';
  if(l == 0) return 'b';
  if( (l->l != 0) && (r->l != 0)) return 'n'; // no grandchildren
  return 'i';
}

tree::tree_p tree::get_ptr(int nid)
{
  if(this->get_nid() == nid) return this; //found it
  if(!l) return 0; //no children, did not find it
  tree_p lp = l->get_ptr(nid);
  if(lp) return lp; //found on left
  tree_p rp = r->get_ptr(nid);
  if(rp) return rp; //found on right
  return 0; //never found it
}

tree::tree_cp tree::get_bn(double* x_cont, int* x_cat)
{
  if(!l) return this; // node has no left child, so it must be a leaf
  int l_count = 0;
  int r_count = 0;
  
  if(!rule.is_cat){
    // axis-aligned rule
    if(x_cont != 0){
      if(x_cont[rule.v_aa] < rule.c) return l->get_bn(x_cont, x_cat);
      else if(x_cont[rule.v_aa] >= rule.c) return r->get_bn(x_cont, x_cat);
      else{
        Rcpp::Rcout << "[get_bn]: could not resolve continuous decision rule: " <<std::endl;
        Rcpp::Rcout << "    x = " << x_cont[rule.v_aa] << " and c = " << rule.c << std::endl;
        return 0;
      }
    } else{
      Rcpp::Rcout << "[get_bn]: encountered continuous decision rule but no continuous predictors were supplied!" << std::endl;
      return 0;
    }
  } else{
    // categorical decision rule
    if(x_cat != 0){
      l_count = rule.l_vals.count(x_cat[rule.v_cat]);
      r_count = rule.r_vals.count(x_cat[rule.v_cat]);
      
      if(l_count == 1 && r_count == 0) return l->get_bn(x_cont, x_cat);
      else if(l_count == 0 && r_count == 1) return r->get_bn(x_cont, x_cat);
      else if(l_count == 0 && r_count == 0){
        Rcpp::Rcout << "[get_bn]: could not find value of categorical predictor in either left or right cutset!" << std::endl;
        return 0;
      } else{ // l_count == 1 & r_count == 1
        Rcpp::Rcout << "[get_bn]: value of categorical predictor in both left & right cutset!" << std::endl;
        return 0;
      }
    } else{
      Rcpp::Rcout << "[get_bn]: encountered categorical decision rule but no categorical predictors were supplied!" << std::endl;
      return 0;
    }
  }
}


double tree::evaluate(double* x_cont, int* x_cat)
{
  tree::tree_cp bn = get_bn(x_cont, x_cat);
  if(bn == 0) return std::nan(""); // when we don't have a valid bottom node, return nan
  else return bn->get_mu();
}

// birth
void tree::birth(int nid, rule_t rule)
{
  tree_p np = get_ptr(nid); // get the pointer to the node being split
  if(!np){
    Rcpp::Rcout << "Trying birth at node w/ id " << nid << " but pointer is invalid!" << std::endl;
    Rcpp::stop("[birth]: invalid pointer to node!");
  }
  if(np->l) Rcpp::stop("[birth]: trying to split a node that already has children!");
  
  tree_p l = new tree; // initialize the new tree that will be the left child of nid
  l->mu = 0.0; // we will overwrite this later
  l->rule.clear();


  
  tree_p r = new tree; // initialize the new tree that will be the left child of nid
  r->mu = 0.0; // we will overwrite this later
  r->rule.clear();
  
  np->l = l;
  np->r = r;
  
  np->rule.is_cat = rule.is_cat;
  np->rule.v_aa = rule.v_aa;
  np->rule.c = rule.c;
  np->rule.v_cat = rule.v_cat;
  np->rule.l_vals = rule.l_vals;
  np->rule.r_vals = rule.r_vals;
  np->mu = 0.0; // we will overwrite this later
  l->p = np;
  r->p = np;
  
}

// perform a death
void tree::death(int nid)
{
  tree_p nb = get_ptr(nid);
  if(!nb) Rcpp::stop("[death]: missing pointer for nid!");
  
  if(nb->is_nog()){
    delete nb->l;
    delete nb->r;
    
    nb->l = 0; // nb now has no children so set corresponding pointers to 0
    nb->r = 0; // nb now has no children so set corresponding pointers to 0
    nb->mu = 0.0; // this will be over-written when we update the mu's
    // reset the rule values
    nb->rule.clear();
  } else Rcpp::stop("[death]: cannot perform death move on a node with grandchildren.");
}

void tree::get_rg_aa(int &v, double &c_lower, double &c_upper){
  if(p){
    if(!p->rule.is_cat && p->rule.v_aa == v){
      // my parent does an axis-aligned split on v
      if(this == p->l){
        // this is the left child of its parent
        // this's parent splits on v so any observation in this should be less p->rule.c
        // p->rule.c establishes an upper bound on the valid cutpoints for v in this
        // If c_upper is already smaller than p_rule.c then we don't want to update c_upper
        if(p->rule.c < c_upper) c_upper = p->rule.c;
      } else if(this == p->r){
        // this is the right child of its parent
        // this's parent splits on v so any observation in this should be greater than p->rule.c
        // p->rule.c is a lower bound on the valid cutpoints for v in this
        // If c_lower is already greater than p->rule.c then we don't want to update c_lower
        if(p->rule.c > c_lower) c_lower = p->rule.c;
      } else Rcpp::stop("[get_rg_aa]: this was not equal to either left or right child of its parent..."); // we should *never* hit this spot
    }
    p->get_rg_aa(v, c_lower, c_upper); // actually recurse up the tree
  }
}

void tree::get_rg_cat(int &v, std::set<int> &levels){
// recruse up the tree:
//    1. we will initialize levels to be the set of all levels for x_cat[v]
//    2. if this has a parent that splits on v, check if this is p->l or p->r.
//        * if this is the left child of a node that splits on v, then replace levels with p->rule.l_vals & set recurse = false
//
  bool recurse = true;
  if(p){
    // this has a parent.
    if(p->rule.is_cat && p->rule.v_cat == v){
      // parent of this splits on v
      if(this == p->l){
        levels.clear();
        for(set_it it = p->rule.l_vals.begin(); it != p->rule.l_vals.end(); ++it) levels.insert(*it);
        recurse = false; // no need to continue recursing up the tree
      } else if(this == p->r){
        levels.clear();
        for(set_it it = p->rule.r_vals.begin(); it != p->rule.r_vals.end(); ++it) levels.insert(*it);
        recurse = false; // no need to continue recursing up the tree
      } else Rcpp::stop("[get_rg_cat]: this was not equal to either left or right child of its parent!");
    }
    if(recurse) p->get_rg_cat(v, levels);
  }
}
void tree::get_anc_v_cat(std::vector<int> &anc_v)
{
  if(p){
    if(p->rule.is_cat){
      anc_v.push_back(p->rule.v_cat);
      p->get_anc_v_cat(anc_v);
    }
  }
}

void tree:get_rg_nested_cat(std::set<int> &levels, int &v, tree_prior_info &tree_pi)
{
  bool recurse = true;
  std::set<int>* parent_vals;
  if(p){
    // this has a parent
    if(p->rule.is_cat){
      // parent splits on categorical predictor
      // get parent values
      if(this == p->l) parent_vals = &p->rule.l_vals;
      else if(this == p->r) parent_vals = &p->rule.r_vals;
      else Rcpp::stop("[get_rg_nested_cat]: this was not left nor right child of its parent!");
      
      if(levels.size() == 0){
        // this->p is the first ancestor at which we split on this cluster. We will populate levels for the first time here
        // 3 possibilities: 1. Parent splits on v; 2. Parent splits on higher res than v (e.g., blockgroup vs tract); 3. Parent splits on lower res than v (tract vs blockgroup)
        // When we recurse up the tree, we don't want to over-include things. So let's add some checks to see whether levels is empty or not
        
        if(p->rule.v_cat == v){
          // We inherit the values from the appropriate branch. At this point, we can stop the recursion
          for(std::set<int>::iterator it = parent_vals->begin(); it != parent_vals->end(); ++it) levels.insert(*it);
          recurse = false;
        } else{
          // we now loop over the hi_lo maps to find the ones involving v and p->rule.v_cat
          
          std::vector<hi_lo_map>::iterator nest_it = tree_pi.nesting->begin();
          for(; nest_it != tree_pi.nesting->end(); ++nest_it){
            if(nest_it->hi == v && nest_it->lo == p->rule.v_cat){
              // p splits on lower resolution than v; we will inherit all the high resolution values contained in appropriate branch
              // now we loop over lo-res values in appropriate branch and get each of the hi-res values
              for(std::set<int>::iterator lo_it = parent_vals->begin(); lo_it != parent_vals->end(); ++lo_it){
                // The lo-res value *lo_it is an available at the node. We need to find the corresponding hi-resolution values
                std::map<int, std::set<int>>::iterator hi_val_it = nest_it->map.find(*lo_it);
                if(hi_val_it != nest_it->map.end()){
                  for(std::set<int>::iterator it = hi_val_it->second.begin(); it != hi_val_it->second.end(); ++it) levels.insert(*it);
                }
              } // closes loop over the lo-res values available from parent
              recurse = false; // no need to continue recursion past a lower-resolution ancestor
              break; // no need to continue looping over elements of nesting.
            }
            if(nest_it->hi == p->rule.v_cat && nest_it->lo == v){
              // p splits on higher resolution than v;
              // loop over hi-res values from parent and add the lo-res value containing it to levels.
              std::map<int, std::set<int>>::iterator lo_it = nest_it->map.begin();
              for(std::set<int>::iterator hi_it = parent_vals->begin(); hi_it != parent_vals->end(); ++hi_it){
                lo_it = nest_it->map.begin();
                for(; lo_it != nest_it->map.end(); ++lo_it){
                  if(lo_it->second.count(*hi_it) == 1){
                    // lo_it->first is the lo-res value containing *hi_it
                    // we need to dump lo_it->first into the levels
                    // we can also break out of this loop
                    levels.insert(lo_it->first);
                    break;
                  }
                } // closes loop looking for the lo-res value containing current hi-res value
              } // closes loop over available hi-res values from parent
              recurse = true; // we need to continue recursion to make sure we don't over-include
              break;
            }  // closes if/else checking whether p splits on higher or lower resolution than v
          } // closes loop over hi-lo maps in the family
        } // closes if/else checking if p splits on v or not
      } else{
        // levels is not empty. earlier in the recursion, we encountered a split
        // on a variable connected to v in nest_graph
        // we cannot simply inherit values
        // if p->rule.v_cat = v, then we need to remove anything from levels not contained in the appropriate branch
        // and if p->rule.v_cat is higher resolution than v, we also need to
        // remove anything from levels not contained in appropriate branch
        
        std::set<int> parent_levels;

        if(p->rule.v_cat == v){
          for(std::set<int>::iterator it = parent_vals->begin(); it != parent_vals->end(); ++it) parent_levels.insert(*it);
          recurse = false; // after this point we need not continue recursion up the tree
        } else{
          std::vector<hi_lo_map>::iterator nest_it = tree_pi.nesting->begin();
          for(; nest_it != tree_pi.nesting->end(); ++nest_it){
            if(nest_it->hi == p->rule.v_cat && nest_it->lo == v){
              // loop over hi-res values in parent_vals and add them to temporary parent_levels
              std::map<int, std::set<int>>::iterator lo_it = nest_it->map.begin();
              for(std::set<int>::iterator hi_it = parent_vals->begin(); hi_it != parent_vals->end(); ++hi_it){
                lo_it = nest_it->map.begin();
                for(; lo_it != nest_it->map.end(); ++lo_it){
                  if(lo_it->second.count(*hi_it) == 1){
                    // lo_it->first is the lo-res value containing *hi_it
                    // we need to dump lo_it->first into the levels
                    // we can also break out of this loop
                    parent_levels.insert(lo_it->first);
                    break;
                  }
                } // closes loop looking for the lo-res value containing current hi-res value
              } // closes loop over available hi-res values from parent
              recurse = true; // we should still continue recursion to make sure we don't over-include levels
              break;
            } // closes if checking whehter parent split on higher-resolution than v in the branch
            if(nest_it->hi == v && nest_it->lo == p->rule.v_cat){
              // we have split on a lower-resolution before and we can just inherit these values.
              // this bit is redundant but it will allow us to do the intersection by hand
              for(std::set<int>::iterator lo_it = parent_vals->begin(); lo_it != parent_vals->end(); ++lo_it){
                // The lo-res value *lo_it is an available at the node. We need to find the corresponding hi-resolution values
                std::map<int, std::set<int>>::iterator hi_val_it = nest_it->map.find(*lo_it);
                if(hi_val_it != nest_it->map.end()){
                  for(std::set<int>::iterator it = hi_val_it->second.begin(); it != hi_val_it->second.end(); ++it) parent_levels.insert(*it);
                }
              } // closes loop over the lo-res values available from parent
              recurse = false; // no need to continue recursion past a lower-resolution ancestor
              break;
            } // closes if checking whether parent split on lower-resolution than v in the branch
          } // closes loop over the hi_lo maps in the family looking for the one with v and p->rule.v_cat
        } // closes if/else checking if we split on v or a high-resolution variable
        
        // at this point, levels may contain more than what we want. we need to include only those values in both levels & parent_levels.
        std::set<int> new_levels;
        for(std::set<int>::iterator it = levels.begin(); it != levels.end(); ++it){
          if(parent_levels.count(*it) == 1) new_levels.insert(*it);
        }
        // now overwrite levels to include everything in new_levels.
        levels.clear();
        for(std::set<int>::iterator it = new_levels.begin(); it != new_levels.end(); ++it) levels.insert(*it);
        
        
      } // closes if/else checking if levels is currently empty.
      
      
      
      std::vector<hi_lo_map>::iterator nest_it = tree_pi.nesting->begin();
      for(; nest_it != tree_pi.nesting->end(); ++nest_it){
        if(nest_it->hi == v && nest_it->lo == p->rule.v_cat){
          // p splits on lower resolution than v; we will inherit all the high resolution values contained in appropriate branch
          // now we loop over lo-res values in appropriate branch and get each of the hi-res values
          for(std::set<int>::iterator lo_it = parent_vals->begin(); lo_it != parent_vals->end(); ++lo_it){
            // The lo-res value *lo_it is an available at the node. We need to find the corresponding hi-resolution values
            std::map<int, std::set<int>>::iterator hi_val_it = nest_it->map.find(*lo_it);
            if(hi_val_it != nest_it->map.end()){
              for(std::set<int>::iterator it = hi_val_it->second.begin(); it != hi_val_it->second.end(); ++it) levels.insert(*it);
            }
          } // closes loop over the lo-res values available from parent
          recurse = false; // no need to continue recursion past a lower-resolution ancestor
          break; // no need to continue looping over elements of nesting.
        } else if(nest_it->hi == p->rule.v_cat && nest_it->lo == v){
          
        }
      }
      
      
      
    }
  }
}


//private functions

//copy tree o to tree n
void tree::cp(tree_p n, tree_cp o)
//assume n has no children (so we don't have to kill them)
//recursion down
{
  if(n->l) Rcpp::stop("[cp]:tree n has children.");
  // if we haven't stopped by now, it's valid to continue to copying

  n->mu = o->mu;
  // not 100% sure if it's valid to do n->rule = o->rule
  // but better to be safe and deliberately copy all members
  n->rule.is_cat = o->rule.is_cat;
  n->rule.v_aa = o->rule.v_aa;
  n->rule.c = o->rule.c;
  n->rule.v_cat = o->rule.v_cat;
  n->rule.l_vals = o->rule.l_vals;
  n->rule.r_vals = o->rule.r_vals;
  
  if(o->l){
    // if o has children
    n->l = new tree; // create new tree for n's left child
    (n->l)->p = n; // assign parent of n's left child as n
    cp(n->l, o->l); // recurse for left child
    n->r = new tree;
    (n->r)->p = n;
    cp(n->r, o->r);
  }
}

