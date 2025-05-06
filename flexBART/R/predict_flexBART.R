predict.flexBART <- function(object, newdata,
                             verbose = FALSE, print_every = 50)
{

  # validate fit; it should be class flexBART_fit
  if(class(object) != "flexBART_fit"){
    stop("object must be of class 'flexBART_fit'. If using output of flexBART or probit_flexBART, pass only the element named 'fit' ")
  }
  validate_flexBART_fit(object)
  
  n <- nrow(newdata)
  cov_ensm <- object[["cov_ensm"]]
  R <- ncol(cov_ensm)
  ###############################
  # Build X_cont
  ###############################
  if(object$dinfo$p_cont > 0){
    X_cont <- matrix(nrow = n, ncol = object$dinfo$p_cont,
                     dimnames = list(c(), object$dinfo$cont_names))
    for(j in object$dinfo$cont_names){
      x_max <- object$dinfo$x_max[j]
      x_min <- object$dinfo$x_min[j]
      if(any(newdata[,j] > x_max) || any(newdata[,j] < x_min)){
        warning(paste0("[predict_flexBART]: found value for", j, "outside range of training data."))
      }
      X_cont[,j] <- 
        convert_continuous(x = newdata[,j],
                           x_min = x_min,
                           x_max = x_max,
                           discrete = is.null(object$dinfo$x_sd[j]))$std_x
    }
  } else{
    X_cont <- matrix(0, nrow = 1, ncol = 1)
  }
  ###############################
  # Build X_cat
  ###############################
  if(object$dinfo$p_cat > 0){
    X_cat <- matrix(nrow = n, ncol = object$dinfo$p_cat,
                    dimnames = list(c(), object$dinfo$cat_names))
    for(j in object$dinfo$cat_names){
      X_cat[,j] <- 
        convert_categorical(x = newdata[,j], name = j,
                            mapping = object$dinfo$cat_mapping_list[[j]])
    }
  } else{
    X_cat <- matrix(0L, nrow = 1, ncol = 1)
  }
  
  if(R == 1){
    tmp <- 
      .single_ensm_predict(tree_draws = object[["trees"]],
                           tX_cont = t(X_cont),
                           tX_cat = t(X_cat),
                           M = object[["M"]][1],
                           probit = object[["is.probit"]],
                           verbose = verbose,
                           print_every = print_every)
    if(!object[["is.probit"]]) output <- object$scaling_info$y_mean + object$scaling_info$y_sd * tmp
    else output <- tmp
  } else{
    ###############################
    # Build Z
    ###############################
    Z <- matrix(1, nrow = n, ncol = R,
                dimnames = list(c(), colnames(cov_ensm)))
    if(any(!is.na(colnames(cov_ensm)))){
      # There are non-intercept terms
      nonint_names <- colnames(cov_ensm)[!is.na(colnames(cov_ensm))]
      for(j in nonint_names){
        Z[,j] <- newdata[,j]
      }
    }
    ########################################
    # Standardize Z
    ########################################
    for(r in 1:R){
      if(!all(Z[,r] == 1)){
        Z[,r] <- (Z[,r] - object$scaling_info$z_mean[r])/object$scaling_info$z_sd[r]
      }
    }
    
    tmp <-
      .multi_ensm_predict(tree_draws = object[["trees"]],
                          tZ = t(Z),
                          tX_cont = t(X_cont),
                          tX_cat = t(X_cat),
                          M_vec = object[["M"]],
                          verbose = verbose, print_every = print_every)
    yhat <- 
      object$scaling_info$y_mean + 
      object$scaling_info$y_sd * tmp[["fit"]]
    
    beta_samples <-  
      rescale_beta(tmp[["raw_beta"]], 
                   y_mean = object$scaling_info$y_mean,
                   y_sd = object$scaling_info$y_sd,
                   z_mean = object$scaling_info$z_mean, 
                   z_sd = object$scaling_info$z_sd,
                   z_col_id = object$scaling_info$z_col_id)
    output <- list()
    output[["yhat"]] <- yhat
    output[["beta"]] <- beta_samples
    output[["raw_beta"]] <- tmp[["raw_beta"]]
    #output[["yhat_raw"]] <- tmp[["fit"]]
    #output[["Z"]] <- Z
  } # closes if/else checking whether it is single ensemble or multiple ensembles
  return(output)
}
  