library(igraph)
library(RColorBrewer)
library(scales)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("flexBART/src/rflexBART.cpp")

K <- 10
g <- erdos.renyi.game(K, 0.1, type = "gnp")
counter <- 0
while( (!is_connected(g)) & counter < 100){
  g <- erdos.renyi.game(K, 0.1, type = "gnp")
  counter <- counter + 1
}
if(!is_connected(g)) stop("graph isn't connected")

my_layout <- layout.circle(g)
plot(g, layout = my_layout)
A <- as_adjacency_matrix(g, sparse = FALSE, names = FALSE)
A_lower <- A
A_lower[upper.tri(A_lower)] <- 0


cat_levels_list <- list(c(0:(K-1)))
adj_support_list <- list(which(A_lower != 0) - 1) # C++ is 1-indexed

X_cont <- matrix(1, nrow = 1, ncol = 1)
X_cat <- matrix(as.integer(0:(K-1)), nrow = K, ncol = 1)

#
M <- 50000

# first let's ignore adjacency entirely
# when we split, we partition the vertices of the induced subgraph
# into two sets uniformly at random (conditionally on not producing a trivial partition)
no_adj_trees <- .draw_ensemble(tX_cont = t(X_cont),
                         tX_cat = t(X_cat),
                         cat_levels_list = cat_levels_list,
                         adj_support_list = adj_support_list,
                         mst_split = FALSE, mst_reweight = FALSE,
                         alpha = 0.95, beta = 2, mu0 = 0, tau = 1,
                         prob_aa = 0, prob_rc = 0, M = M, verbose = TRUE, print_every = M/10)

# now let's respect adjacency: we split by 
# (1) drawing uniform weights for edges of induced subgraph
# (2) computing MST for this weighted subgraph
# (3) delete one edge from the MST uniformly at random
unif_edge_trees <- .draw_ensemble(tX_cont = t(X_cont),
                            tX_cat = t(X_cat),
                            cat_levels_list = cat_levels_list,
                            adj_support_list = adj_support_list,
                            mst_split = TRUE, mst_reweight = FALSE,
                            alpha = 0.95, beta = 2, mu0 = 0, tau = 1,
                            prob_aa = 0, prob_rc = 0, M = M, verbose = TRUE, print_every = M/10)
# same as above but instead of deleting an edge uniformly at random
# delete edge with prob proportional to the size of the *smallest* cluster
# that is produced when edge is deleted
# this biases us away from creating singletons
size_biased_trees <- .draw_ensemble(tX_cont = t(X_cont),
                              tX_cat = t(X_cat),
                              cat_levels_list = cat_levels_list,
                              adj_support_list = adj_support_list,
                              mst_split = TRUE, mst_reweight = TRUE,
                              alpha = 0.95, beta = 2, mu0 = 0, tau = 1,
                              prob_aa = 0, prob_rc = 0, M = M, verbose = TRUE, print_every = M/10)

# .draw_ensemble returns a list. one element, named tree_fits, is a matrix of size n x M,
# recording the fit of each tree, where n is number of rows in X_cat
# 

get_kernel <- function(n, tree_fits){
  # tree_fits is a n x M matrix containing fit of individual trees (columns index trees)
  # if observations i and ii are in the same cluster, they will have the same value in a column
  # so mean(tree_fits[i,] == tree_fits[ii,]) is approx. proportion of times k and kk are clustered together
  
  co_cluster <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(ii in 1:i){
      co_cluster[i,ii] <- mean(tree_fits[i,] == tree_fits[ii,])
      co_cluster[ii,i] <- co_cluster[i,ii]
    }
  }
  return(co_cluster)
}

no_adj_kernel <- get_kernel(K, no_adj_trees$tree_fits)
unif_edge_kernel <- get_kernel(K, unif_edge_trees$tree_fits)
size_biased_kernel <- get_kernel(K, size_biased_trees$tree_fits)

View(round(no_adj_kernel, digits = 3))
View(round(unif_edge_kernel, digits = 3))
View(round(size_biased_kernel, digits = 3))

# graph laplacian 
