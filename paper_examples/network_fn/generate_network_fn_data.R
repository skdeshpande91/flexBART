vertex_id_all <- rep(1:n_vertex, each = n_obs)
X_cont_all <- matrix(runif(n_vertex*n_obs*p_cont,min = -1,max = 1),
                     nrow = n_vertex*n_obs,
                     ncol = p_cont)
X_cat_all <- matrix(vertex_id_all - 1, ncol = 1)

bart_df_all <- data.frame(X_cont_all, vertex = factor(vertex_id_all))
bart_mm_all <- 
  dbarts::makeModelMatrixFromDataFrame(bart_df_all, drop = FALSE)
if(!is.na(embed_dim)){
  ase_mm_all <- 
    dbarts::makeModelMatrixFromDataFrame(data.frame(X_cont_all, X_ase[vertex_id_all,1:embed_dim]), drop = FALSE)
} else{
  ase_mm_all <- 
    dbarts::makeModelMatrixFromDataFrame(data.frame(X_cont_all), drop = FALSE)
}


mu_all <- mu_true(X_cont_all, vertex_id_all)

Y_all <- mu_all + sigma * rnorm(n_vertex*n_obs, mean = 0, sd = 1)


test_vertices <- sample(1:n_vertex, size = floor(0.1 * n_vertex), replace = FALSE)
train_vertices <- (1:n_vertex)[-test_vertices]
train_index <- which(vertex_id_all %in% train_vertices)


Y_train <- Y_all[train_index]
vertex_id_train <- vertex_id_all[train_index]
X_cont_train <- X_cont_all[train_index, ]
X_cat_train <- matrix(X_cat_all[train_index, ], ncol = 1)
bart_df_train <- bart_mm_all[train_index,]
bart_mm_train <- bart_mm_all[train_index,]
ase_mm_train <- ase_mm_all[train_index,]

mu_train <- mu_all[train_index]

n_test <- 500
#n_test <- 100
#n_test <- 250
vertex_id_test <- rep(1:n_vertex, each = n_test)

X_cont_test <- matrix(runif(n_vertex * n_test * p_cont, min = -1, max = 1), ncol = p_cont)
X_cat_test <- matrix(vertex_id_test - 1, ncol = 1)
bart_df_test <- data.frame(X_cont_test, vertex = factor(vertex_id_test))

bart_mm_test <- 
  dbarts::makeModelMatrixFromDataFrame(bart_df_test, drop = FALSE)
if(!is.na(embed_dim)){
  ase_mm_test <- 
    dbarts::makeModelMatrixFromDataFrame(data.frame(X_cont_test, X_ase[vertex_id_test,1:embed_dim]), drop = FALSE)
} else{
  ase_mm_test <- 
    dbarts::makeModelMatrixFromDataFrame(data.frame(X_cont_test), drop = FALSE)
}

mu_test <- mu_true(X_cont_test, vertex_id_test)

test_index_1 <- which(vertex_id_test %in% train_vertices)
test_index_2 <- which(vertex_id_test %in% test_vertices)


mu_test1 <- mu_test[test_index_1]
mu_test2 <- mu_test[test_index_2]



