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

