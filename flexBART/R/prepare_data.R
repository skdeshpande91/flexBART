prepare_data <- function(train_data,
                         outcome_name,
                         cov_ensm, 
                         test_data = NULL,
                         probit = FALSE,
                         ...)
{
  
  ###############################
  # Check if user supplied certain arguments
  ###############################
  usr_args <- list(...)
  usr_names <- names(usr_args) # equivalent to `...names()`
  
  # padding for re-scaling continuous covariates
  # default is 0.2 SD's
  pad <- NULL
  if("pad" %in% usr_names) pad <- usr_args[["pad"]]
  
  # for detecting whether a continuous variable is actually ordinal
  # we look at the number of unique differences. if it is below (resp. above)
  # threshold, we say that the variable is ordinal (resp. continuous)
  n_unik_diffs <- NULL
  if("n_unik_diffs" %in% usr_names) n_unik_diffs <- usr_args[["n_unik_diffs"]]
  
  # are there adjacency matrices describing structure of categorial levels
  #adjacency_list <- NULL
  tmp_adj <- pmatch(usr_names, table = "adjacency_list", duplicates.ok = FALSE)
  if(any(!is.na(tmp_adj))){
    ix <- which(!is.na(tmp_adj))
    tmp_name <- usr_names[[ix]]
    adjacency_list <- usr_args[[tmp_name]]
  } else{
    adjacency_list <- NULL
  }
  
  ###############################
  # Get the problem dimensions
  ###############################
  covariate_names <- rownames(cov_ensm)
  R <- ncol(cov_ensm) # how many ensembles
  p <- nrow(cov_ensm) # how many total covariates (across all ensembles)
  
  n_train <- nrow(train_data)
  n_test <- 0
  if(!is.null(test_data)){
    stopifnot(is.data.frame(test_data))
    if(!all(covariate_names %in% colnames(test_data))){
      cat("[prepare_data]: Following covariates used for training not found in supplied testing data: \n")
      cat(covariate_names[!covariate_names %in% colnames(test_data)], "\n")
      stop("[prepare_data]: test_data must contain all training covariates")
    }
    cat("[prepare_data]: Using both training & testing data to get covariate information \n")
    if(p == 1){
      cov_data <- 
        data.frame(c(train_data[,covariate_names[1]], test_data[,covariate_names[1]]))
      colnames(cov_data) <- covariate_names[1]
    } else{
      cov_data <- 
        rbind(train_data[,covariate_names], test_data[,covariate_names])
    }
    n_test <- nrow(test_data)
  } else{
    if(p == 1){
      cov_data <- data.frame(train_data[,covariate_names[1]])
      colnames(cov_data) <- covariate_names[1]
    } else{
      cov_data <- train_data[,covariate_names]
    }
  }
  dinfo <- 
    get_covariate_info(cov_data = cov_data, pad = pad, n_unik_diffs = n_unik_diffs)
  dinfo$outcome_name <- outcome_name
  
  ###############################
  # Create training object
  ###############################
  trinfo <- new_flexBART_training()
  ###############################
  # Standardize outcome & check for missingness
  ###############################
  y <- train_data[,outcome_name]
  if(any(is.na(y))){
    message(paste("[prepare_data]: Detected missing values in", outcome_name))
    stop("[prepare_data]: flexBART does not yet support missing outcomes. 
         Please re-run after removing observations w/ missing outcomes.")
  }
  if(!probit){
    trinfo$y_mean <- mean(y)
    trinfo$y_sd <- sd(y)
    trinfo$std_Y <- (y - trinfo$y_mean)/trinfo$y_sd 
  } else{
    if(!is.integer(y)){
      stop("For probit regression ", outcome_name, " must be an integer")
    }
    if(!all(y %in% c(0L,1L))){
      stop("For probit regression ", outcome_name, " must take values in 0L or 1L")
    }
    trinfo$std_Y <- y
  }

  
  ###############################
  # Build X_cont_train (if applicable)
  ###############################
  if(dinfo$p_cont > 0){
    # loop over all the continuous covariates and re-scale as needed
    trinfo$X_cont <- 
      matrix(nrow = n_train, ncol = dinfo$p_cont,
             dimnames = list(c(), dinfo$cont_names))
    tmp_cutpoints_list <- list()
    for(j in dinfo$cont_names){
      tmp <- 
        convert_continuous(x = train_data[,j],
                           x_min = dinfo$x_min[j],
                           x_max = dinfo$x_max[j],
                           discrete = is.null(dinfo$x_sd[j]))
      trinfo$X_cont[,j] <- tmp$std_x
      if(is.null(tmp$cutpoints)){
        tmp_cutpoints_list[j] <- list(NULL)
      } else{
        tmp_cutpoints_list[[j]] <- tmp$cutpoints
      }
    } # closes loop over all continuous covariates
    
    # check if there are any cutpoints
    if(any(sapply(tmp_cutpoints_list, FUN = length) != 0)){
      trinfo$cutpoints <- tmp_cutpoints_list
    }
  } # closes if checking for continuous covariates
  
  ###############################
  # Build X_cat_train (if applicable)
  ###############################
  if(dinfo$p_cat > 0){
    trinfo$X_cat <- 
      matrix(NA, nrow = n_train, ncol = dinfo$p_cat,
             dimnames = list(c(), dinfo$cat_names))
    trinfo$cat_levels_list <- list()
    
    for(j in dinfo$cat_names){
      trinfo$cat_levels_list[[j]] <- 
        dinfo$cat_mapping_list[[j]][,"integer_coding"]
      trinfo$X_cat[,j] <- 
        convert_categorical(x = train_data[,j], 
                            name = j, 
                            mapping = dinfo$cat_mapping_list[[j]])
    } # closes loop over categorical covariates
    
    # now parse adjacency structures
    if(!is.null(adjacency_list)){
      if(!all(names(adjacency_list) %in% dinfo$cat_names)){
        stop("[preprocess]: all names of adjacency_list must be categorical variable names")
      }
      trinfo$edge_mat_list <- 
        parse_adjacency(adjacency_list, dinfo)
    } # closes if checking if there are any network-structured predictors
  } # closes if checking if there are categorical predictors
  
  ###############################
  # Build X_cont_test and X_cat_test (if applicable)
  ###############################
  teinfo <- new_flexBART_testing()
  if(!is.null(test_data)){
    if(dinfo$p_cont > 0){
      teinfo$X_cont <- 
        matrix(NA, nrow = n_test, ncol = dinfo$p_cont,
               dimnames = list(c(), dinfo$cont_names))
      for(j in dinfo$cont_names){
        if(is.null(trinfo$cutpoints[[j]])){
          x_max <- dinfo$x_max[j]
          x_min <- dinfo$x_min[j]
          if(any(test_data[,j] > x_max) | any(test_data[,j] < x_min)){
            warning(paste("[prepare_data]: found test set value for", j, 
                          "outside range found in training. Consider increased padding"))
          }
        } else{
          # j is a gridded/discrete ordinal variable
          x_max <- max(trinfo$cutpoints[[j]])
          x_min <- min(trinfo$cutpoints[[j]])
          if(any(test_data[,j] > x_max) | any(test_data[,j] < x_min)){
            warning(paste("[prepare_data]: found test set value for", j, 
                          "outside range of cutpoints used in training.",
                          "Consider setting cutpoints_list manually"))
          }
        }
        teinfo$X_cont[,j] <- 
          convert_continuous(x = test_data[,j],
                             x_min = dinfo$x_min[j],
                             x_max = dinfo$x_max[j],
                             discrete = is.null(dinfo$x_sd[j]))$std_x
      } # closes loop over continuous predictors
    } # closes if checking if there are continuous predictors
    if(dinfo$p_cat > 0){
      teinfo$X_cat <- 
        matrix(nrow = n_test, ncol = dinfo$p_cat,
               dimnames = list(c(), dinfo$cat_names))
      for(j in dinfo$cat_names){
        teinfo$X_cat[,j] <- 
          convert_categorical(x = test_data[,j], name = j,
                              mapping = dinfo$cat_mapping_list[[j]])
      } # closes loop over categorical covariates j
    } # closes if checking there are categorical covariates
  } # closes if checking if there is testing data
  
  ###############################
  # Detect nesting structure
  ###############################
  if(dinfo$p_cat > 0){
    if(!is.null(test_data)){
      cat("[preprocess]: Using training & testing to detect nesting structure among categorical variables\n")
      tmp_nest <- parse_nesting(rbind(trinfo$X_cat, teinfo$X_cat), dinfo)
    } else{
      tmp_nest <- parse_nesting(trinfo$X_cat, dinfo)
    }
    if(!is.null(tmp_nest$nest_list)){
      trinfo$nest_list <- tmp_nest$nest_list
    }
  }

  ###############################
  # Build Z_train
  ###############################
  R <- ncol(cov_ensm)
  trinfo$Z <- 
    matrix(1, nrow = n_train, ncol = R,
           dimnames = list(c(), colnames(cov_ensm)))
  if(n_test > 0){
    teinfo$Z <- 
      matrix(1, nrow = n_test, ncol = R,
             dimnames = list(c(), colnames(cov_ensm)))
  }
  
  if(any(!is.na(colnames(cov_ensm)))){
    # There are non-intercept terms
    nonint_names <- colnames(cov_ensm)[!is.na(colnames(cov_ensm))]
    for(j in nonint_names){
      trinfo$Z[,j] <- train_data[,j]
      if(n_test > 0) teinfo$Z[,j] <- test_data[,j]
    }
  }
  
  ########################################
  # Standardize Z
  ########################################
  
  trinfo$R <- R
  trinfo$z_mean <- rep(NA, times = R)
  trinfo$z_sd <- rep(NA, times = R)
  
  for(r in 1:R){
    if(all(trinfo$Z[,r] == 1)){
      trinfo$z_mean[r] <- 0
      trinfo$z_sd[r] <- 1
    } else{
      trinfo$z_mean[r] <- mean(trinfo$Z[,r])
      trinfo$z_sd[r] <- sd(trinfo$Z[,r])
      trinfo$Z[,r] <- (trinfo$Z[,r] - trinfo$z_mean[r])/trinfo$z_sd[r]
    }
    # remember to standardize Z_test
    if(n_test > 0){
      teinfo$Z[,r] <- (teinfo$Z[,r] - trinfo$z_mean[r])/trinfo$z_sd[r]
    }
  }
  
  ######################################
  # Figure out which beta_j's are identified
  # maybe pass a message to user that things are not identified so they should
  # tread cautiously
  ######################################
  
  trinfo$z_col_id <- get_z_col_id(trinfo$Z)
  if(n_test > 0){
    if(!identical(trinfo$z_col_id, get_z_col_id(teinfo$Z))){
      stop("Z_test must have same patterns of identical columns as Z_train!")
    }
  }
  ########################################
  # Validate that training, testing, and data
  # info objects are built correctly
  ########################################
  
  
  return(list(data_info = dinfo,
              training_info = trinfo,
              testing_info = teinfo))
}