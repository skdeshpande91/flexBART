new_flexBART_training <- function()
{
  out <- list()
  out["std_Y"] <- list(NULL)
  out["y_mean"] <- list(NULL)
  out["y_sd"] <- list(NULL)
  # in the C++ code, we use a 1x1 matrix as a null value for design matrices
  # this avoids having to cast nullable arguments first
  out["X_cont"] <- matrix(0, nrow = 1, ncol = 1) 
  out["X_cat"] <- matrix(0, nrow = 1, ncol = 1) 
  out["Z"] <- matrix(0, nrow = 1, ncol = 1)
  
  out["cutpoints"] <- list(NULL)
  out["cat_levels_list"] <- list(NULL)
  out["edge_mat_list"] <- list(NULL)
  out["nest_list"] <- list(NULL)
  
  out["R"] <- list(NULL)
  out["z_mean"] <- list(NULL)
  out["z_sd"] <- list(NULL)
  out["z_col_id"] <- list(NULL)
  
  structure(out, class = "flexBART_training")
}

validate_flexBART_training <- function(trinfo)
{
  exp_names <- c("std_Y", "y_mean", "y_sd", "X_cont", "X_cat", "Z",
                 "cutpoints",
                 "cat_levels_list", "edge_mat_list",
                 "nest_list",
                 "R", "z_mean", "z_sd", "z_col_id")
  x_uncl <- unclass(trinfo)
  if(!identical(sort(names(x_uncl)), sort(exp_names))){
    stop("[validate_flexBART_training]: trinfo does not have valid names")
  }
  
  # check that std_Y is not null and has mean 0 and variance 1
  if(is.null(trinfo$std_Y)){
    stop("[validate_flexBART_training]: std_Y is null")
  } else{
    if( abs(mean(trinfo$std_Y)) > 1e-12 | abs(sd(trinfo$std_Y - 1)) > 1e-12 ){
      stop("[validate_flexBART_training]: internal standardization of Y failed.")
    }
  }
  
  if(ncol(trinfo$Z) != trinfo$R){
    message("[validate_flexBART_training]: supplied ncol(trinfo$Z) = ", ncol(trinfo$Z), "but R = ", trinfo$R)
    stop("[validate_flexBART_training]: Z must have exactly R columns.")
  }
  if(length(trinfo$z_mean) != trinfo$R){
    message("[validate_flexBART_training]: supplied z_mean has length = ", length(trinfo$z_mean), "but R = ", trinfo$R)
    stop("[validate_flexBART_training]: z_mean must be length R.")
  }
  if(length(trinfo$z_sd) != trinfo$R){
    message("[validate_flexBART_training]: supplied z_sd = ", length(trinfo$z_sd), "but R = ", trinfo$R)
    stop("[validate_flexBART_training]: z_sd must be length R.")
  }
  
  
  for(r in 1:R){
    if(is.na(trinfo$z_mean[r]) != is.na(trinfo$z_sd[r])){
      message("[validate_flexBART_training]: r = ", r, "z_mean[r] =", z_mean[r], " z_sd[r] = ", z_sd[r])
      stop("[validate_flexBART_training]: r-th element of z_mean & z_sd must be both NA or both non-NA")
    }
    
    if(is.na(trinfo$z_mean[r])){
      if(any(abs(trinfo$Z[,r] - 1) > 1e-12)){
        stop("[validate_flexBART_training]: Expected column", r, " of Z to be all 1's but found non-1 entries")
      }
    } else{
      if(abs(mean(trinfo$Z[,r])) > 1e-12){
        message("[validate_flexBART_training]: r = ", r, " mean of Z[,r] = ", mean(trinfo$Z[,r]))
        stop("[validate_flexBART_training]: expected column of Z to have mean 0")
      }
      if(abs(sd(trinfo$Z[,r]) - 1) > 1e-12){
        message("[validate_flexBART_training]: r = ", r, " mean of Z[,r] = ", sd(trinfo$Z[,r]))
        stop("[validate_flexBART_training]: expected column of Z to have sd 1")
      }
    }
  }
  
  if(!is.null(trinfo$cutpoints)){
    if(any(sapply(trinfo$cutpoints, FUN = length) == 1)){
      stop("Only one cutpoint supplied for some variables. Should be NULL or should provide >= 2 cutpoints ")
    }
  }
  
}


new_flexBART_testing <- function()
{
  out <- list()
  out["X_cont"] <- matrix(0, nrow = 1, ncol = 1)
  out["X_cat"] <- matrix(0, nrow = 1, ncol = 1)
  out["Z"] <- matrix(0, nrow = 1, ncol = 1)
  structure(out, class = "flexBART_testing")
}

validate_flexBART_testing <- function(teinfo){
  exp_names <- c("X_cont", "X_cat", "Z")
  x_uncl <- unclass(teinfo)
  if(!identical(sort(names(x_uncl)), sort(exp_names))){
    stop("[validate_flexBART_testing]: teinfo does not have valid names")
  }
}

new_flexBART_data_info <- function()
{
  out <- list()
  out["outcome_name"] <- list(NULL)
  
  out["p_cont"] <- 0
  out["p_cat"] <- 0
  out["p"] <- 0
  out["cont_names"] <- list(NULL)
  out["cat_names"] <- list(NULL)
  
  out["x_min"] <- list(NULL)
  out["x_max"] <- list(NULL)
  out["x_sd"] <- list(NULL)
  out["cat_mapping_list"] <- list(NULL)
  structure(out, class = "flexBART_data_info")
}

validate_flexBART_data_info <- function(dinfo)
{
  exp_names <- c("outcome_name",
                 "p_cont", "p_cat", "p", "cont_names", "cat_names",
                 "x_min", "x_max", "x_sd",
                 "cat_mapping_list")
  x_uncl <- unclass(dinfo)
  if(!identical(sort(names(x_uncl)), sort(exp_names))){
    stop("[validate_flexBART_data_info]: dinfo does not have valid names")
  }
  # check that length(cont_names) == p_cont
  # check that length(cat_names) == p_cat
  # check that cat_mapping_list has exactly p_cat elements
}


new_flexBART_fit <- function()
{
  out <- list()
  out["dinfo"] <- list(NULL)
  out["trees"] <- list(NULL)
  out["is.probit"] <- list(NULL)
  out["scaling_info"] <- list(NULL)
  out["M"] <- list(NULL)
  out["cov_ensm"] <- list(NULL)

  structure(out, class = "flexBART_fit")
}
validate_flexBART_fit <- function(fit)
{
  exp_names <- c("dinfo", "trees", "is.probit",
                 "scaling_info", "M", "cov_ensm")
  x_uncl <- unclass(fit)
  if(!identical(sort(names(x_uncl)), sort(exp_names))){
    stop("[validate_flexBART_fit]: fit does not have valid names")
  }
  
  # for probit, we don't really need y_mean or y_sd to be set and until we get
  # multiple ensemble probit working, no need to 
}