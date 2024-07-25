Y_all <- 
  rep(mu, each = n_obs) + 
  sigma * rnorm(n_vertex*n_obs, mean = 0, sd = 1)
vertex_id_all <- rep(1:n_vertex, each = n_obs)

bart_df_all <- data.frame(vertex = factor(vertex_id_all))
bart_mm_all <- 
  dbarts::makeModelMatrixFromDataFrame(bart_df_all, drop = FALSE)

if(!is.na(embed_dim)){
  ase_mm_all <-
    dbarts::makeModelMatrixFromDataFrame(data.frame(X_ase[vertex_id_all,1:embed_dim]), drop = FALSE)
} else{
  ase_mm_all <-
    dbarts::makeModelMatrixFromDataFrame(data.frame(X_ase[vertex_id_all,1]), drop = FALSE)
}


test_vertices <- sample(1:n_vertex, size = floor(0.1 * n_vertex), replace = FALSE)
train_vertices <- (1:n_vertex)[-test_vertices]
train_index <- which(vertex_id_all %in% train_vertices)

Y_train <- Y_all[train_index]
vertex_id_train <- vertex_id_all[train_index]
X_cat_train <- matrix(vertex_id_train-1, ncol = 1)
bart_df_train <- bart_df_all[train_index,]
bart_mm_train <- bart_mm_all[train_index,]
ase_mm_train <- ase_mm_all[train_index,]


# For testing, make a prediction at every vertex

vertex_id_test <- 1:n_vertex
X_cat_test <- matrix(vertex_id_test-1, ncol = 1)
bart_df_test <- data.frame(vertex = vertex_id_test)
bart_mm_test <- unique(bart_mm_all)
ase_mm_test <- unique(ase_mm_all)





