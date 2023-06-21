network_BART <- function(Y_train,
                         vertex_id_train, # observation i belongs to node vertex_id_train[i].
                         X_cont_train = matrix(0, nrow = 1, ncol = 1),
                         vertex_id_test = NULL,
                         X_cont_test = matrix(0, nrow = 1, ncol = 1),
                         unif_cuts = rep(TRUE, times = ncol(X_cont_train)),
                         cutpoints_list = NULL,
                         A = matrix(0, nrow = 1, ncol = 1),
                         graph_split = TRUE,
                         graph_cut_type = 1,
                         sparse = FALSE,
                         M = 200,
                         nd = 1000, burn = 1000, thin = 1,
                         save_samples = TRUE,
                         save_trees = TRUE, verbose = TRUE, print_every = floor( (nd*thin + burn))/10)
{
  if(!graph_cut_type %in% c(0,1,2,3)) stop("graph_cut_type must be 0,1,2, or 3.")
  
  # pre-processing specific for network stuff
  
  if(nrow(A) == 1) stop("Adjacency matrix A must contain > 1 row (i.e. graph must have > 1 vertex)")
  n_vertex <- nrow(A)
  cat_levels_list <- list(0:(n_vertex-1)) # remember C++ is 0-indexed
  X_cat_train <- matrix(as.integer(vertex_id_train - 1), ncol = 1, dimnames = list(c(), "vertex"))
  if(!is.null(vertex_id_test)){
    X_cat_test <- matrix(as.integer(vertex_id_test - 1), ncol = 1)
  } else{
    X_cat_test <- matrix(0, nrow = 1, ncol = 1)
  }
  
  # get a matrix that lists the edges
  # this is is similar to what would obtained from igraph::get.data.frame(graph, what = "edges")
  tmp_indices <- which(A !=0, arr.ind = TRUE)
  edge_mat <- tmp_indices[tmp_indices[,1] < tmp_indices[,2],] # get only the lower triangle of A
  colnames(edge_mat) <- c("from", "to")
  edge_mat_list <- list(edge_mat-1) # remember C++ is 0-indexed
  
  # the usual pre-processing
  y_mean <- mean(Y_train)
  y_sd <- stats::sd(Y_train)
  std_Y_train <- (Y_train - y_mean)/y_sd # standardize the output
  tau <- (max(std_Y_train) - min(std_Y_train))/(2 * 2 * sqrt(M)) # CGM10 prior sd on all leaf parameters
  nu <- 3
  lambda <- stats::qchisq(0.1, df = nu)/nu
  
  p_cont <- 0
  p_cat <- 0
  cont_names <- c()
  cat_names <- c()
  
  if(length(X_cont_train) > 1){
    p_cont <- ncol(X_cont_train)
    if(is.null(colnames(X_cont_train))){
      cont_names <- paste0("X", 1:p_cont)
    } else{
      cont_names <- colnames(X_cont_train)
    }
  } else{
    cont_names <- c()
  }
  if(length(X_cat_train) > 1){
    p_cat <- ncol(X_cat_train)
    if(is.null(colnames(X_cat_train))){
      cat_names <- paste0("X", (p_cont+1):(p_cont+p_cat))
    } else{
      cat_names <- colnames(X_cat_train)
    }
  } else{
    cat_names <- c()
  }
  pred_names <- c(cont_names, cat_names)
  
  fit <- .flexBART_fit(Y_train = std_Y_train,
                       tX_cont_train = t(X_cont_train),
                       tX_cat_train = t(X_cat_train),
                       tX_cont_test = t(X_cont_test),
                       tX_cat_test = t(X_cat_test),
                       unif_cuts = unif_cuts,
                       cutpoints_list = cutpoints_list,
                       cat_levels_list = cat_levels_list,
                       edge_mat_list = edge_mat_list,
                       graph_split = graph_split,
                       graph_cut_type = graph_cut_type,
                       a_cat = 0, b_cat = 0,
                       rc_split = FALSE, prob_rc = 0, a_rc = 1, b_rc = 1,
                       sparse = sparse, a_u = 0.5, b_u = 1,
                       mu0 = 0, tau = tau, 
                       lambda = lambda, nu = nu,
                       M = M, nd = nd, burn = burn, thin = thin,
                       save_samples = save_samples,
                       save_trees = save_trees, verbose = verbose, 
                       print_every = print_every)

  yhat_train_mean <- y_mean + y_sd * fit$fit_train_mean
  if(save_samples){
    yhat_train_samples <- y_mean + y_sd * fit$fit_train
  }
  if(!is.null(fit$fit_test_mean)){
    yhat_test_mean <- y_mean + y_sd * fit$fit_test_mean
    if(save_samples){
      yhat_test_samples <- y_mean + y_sd * fit$fit_test
    }
  }
  
  
  sigma_samples <- y_sd * fit$sigma
  
  results <- list()
  results[["y_mean"]] <- y_mean
  results[["y_sd"]] <- y_sd
  results[["yhat.train.mean"]] <- yhat_train_mean
  if(save_samples) results[["yhat.train"]] <- yhat_train_samples
  if(!is.null(fit$fit_test_mean)){
    results[["yhat.test.mean"]] <- yhat_test_mean
    if(save_samples) results[["yhat.test"]] <- yhat_test_samples
  }
  results[["sigma"]] <- y_sd * fit$sigma
  
  varcounts <- fit$var_count
  if(length(pred_names) != ncol(varcounts)){
    warning("There was an issue tracking variable names. Not naming columns of varcounts object")
  } else{
    colnames(varcounts) <- pred_names
  }
  results[["varcounts"]] <- varcounts
  
  if(save_trees) results[["trees"]] <- fit$trees
  return(results)
}