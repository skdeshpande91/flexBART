# function to make a graph partition
make_partition <- function(g, K, min_size = 5, max_attempts = 1000){
  if(K == 1) stop("K must be greater than 1")
  m <- length(E(g)) # number of edges
  tmp_g <- g
  E(tmp_g)$weights <- runif(m, min = 0, max = 1)
  st <- mst(tmp_g, weights = E(tmp_g)$weights, algorithm = "prim")
  A_st <- as_adjacency_matrix(graph = st, type = "lower", sparse = FALSE)
  
  edge_list <- which(A_st != 0)
  removal_index <- sample(1:length(edge_list), size = K-1, replace = FALSE) # cut K-1 edges to get K clusters
  A_st[edge_list[removal_index]] <- 0
  new_g <- graph_from_adjacency_matrix(A_st, mode = "lower")
  comps <- components(new_g)
  
  attempts <- 1
  while(min(table(comps$membership)) < min_size & attempts <= max_attempts){
    E(tmp_g)$weights <- runif(m, min = 0, max = 1)
    st <- mst(tmp_g, weights = E(tmp_g)$weights, algorithm = "prim")
    A_st <- as_adjacency_matrix(graph = st, type = "lower", sparse = FALSE)
    
    edge_list <- which(A_st != 0)
    removal_index <- sample(1:length(edge_list), size = K-1, replace = FALSE) # cut K-1 edges to get K clusters
    A_st[edge_list[removal_index]] <- 0
    new_g <- graph_from_adjacency_matrix(A_st, mode = "lower")
    comps <- components(new_g)
    attempts <- attempts + 1
  }
  if(attempts == max_attempts){
    print(paste("Could not create partition with minimum cluster size >= ", min_size))
    return(NULL)
  } else return(comps$membership) # labels the different components
}