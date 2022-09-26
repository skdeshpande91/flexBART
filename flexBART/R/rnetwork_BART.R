rnetwork_BART <- function(M, vertex_id, 
                          X_cont = matrix(0, nrow = 1, ncol = 1),
                          unif_cuts = rep(TRUE, times = ncol(X_cont)),
                          cutpoints_list = NULL,
                          A = matrix(0, nrow = 1, ncol = 1),
                          graph_split = TRUE,
                          graph_cut_type = 0,
                          alpha = 0.95, beta = 2, mu0 = 0, tau = 1/sqrt(M),
                          verbose = TRUE, print_every = floor(M/10))
{
  if(!graph_cut_type %in% c(0,1)) stop("graph_cut_type must be 0 or 1.")
  
  if(nrow(A) == 1) stop("Adjacency matrix A must contain > 1 row (i.e. graph must have > 1 vertex)")
  n_vertex <- nrow(A)
  cat_levels_list <- list(0:(n_vertex-1)) # remember C++ is 0-indexed
  X_cat <- matrix(as.integer(vertex_id - 1), ncol = 1)
 
  
  # get a matrix that lists the edges
  # this is is similar to what would obtained from igraph::get.data.frame(graph, what = "edges")
  tmp_indices <- which(A !=0, arr.ind = TRUE)
  edge_mat <- tmp_indices[tmp_indices[,1] < tmp_indices[,2],] # get only the lower triangle of A
  colnames(edge_mat) <- c("from", "to")
  edge_mat_list <- list(edge_mat-1) # remember C++ is 0-indexed
  
  tree_draws <- .drawEnsemble(tX_cont = t(X_cont),
                              tX_cat = t(X_cat),
                              unif_cuts = unif_cuts,
                              cutpoints_list = cutpoints_list,
                              cat_levels_list = cat_levels_list,
                              edge_mat_list = edge_mat_list,
                              graph_split = graph_split,
                              graph_cut_type = 0,
                              perc_rounds = 0, perc_threshold = 0.5,
                              rc_split = FALSE, prob_rc = 0, a_rc = 1, b_rc = 1,
                              alpha = alpha, beta = beta,
                              mu0 = mu0, tau = tau, M = M,
                              verbose = verbose, print_every = print_every)
  return(tree_draws)
}