###############################
# Object to hold hyperparameters
###############################
new_flexBART_hyper <- function()
{
  out <- list()
  out["M_vec"] <- list(NULL)
  out["alpha_vec"] <- list(NULL)
  out["beta_vec"] <- list(NULL)
  out["mu0_vec"] <- list(NULL)
  out["tau_vec"] <- list(NULL)
  out["graph_cut_type"] <- 2
  out["nest_v"] <- FALSE
  out["nest_v_option"] <- 3
  out["nest_c"] <- FALSE
  out["sigest"] <- 1
  out["sigquant"] <- 0.9
  out["nu"] <- 3
  out["lambda"] <- qchisq(1-out$sigquant, df = out$nu)/out$nu * out$sigest * out$sigest
  out["sparse"] <- TRUE
  out["a_u"] <- 0.5
  out["b_u"] <- 1
  structure(out, class = "flexBART_hyper")
}

###############################
# Check that hyperparameters are valid
###############################
validate_flexBART_hyper <- function(hyper)
{
  exp_names <- 
    c("M_vec", "alpha_vec", "beta_vec", 
      "mu0_vec", "tau_vec", 
      "graph_cut_type", 
      "nest_v", "nest_v_option", "nest_c",
      "sigest", "sigquant",
      "nu", "lambda", 
      "sparse", "a_u", "b_u")
  
  if(!identical(sort(names(hyper)), sort(exp_names))){
    cat("[validate_flexBART_hyper]: supplied hyper names: \n", sort(names(hyper)), "\n")
    cat("[validate_flexBART_hyper]: expected hyper names: \n", sort(exp_names), "\n")
    stop("[validate_flexBART_hyper]: hyper does not have valid names")
  }
  
  ###############################
  # Check that vectors of tree prior parameters
  # have the same length
  ###############################
  tmp_length <- 
    sapply(hyper[c("M_vec", "alpha_vec", "beta_vec", "mu0_vec", "tau_vec")], FUN = length)
  if(length(unique(tmp_length)) != 1){
    stop("[validate_flexBART_hyper]: M_vec, alpha_vec, beta_vec, mu0_vec, and tau_vec must all have the same length!")
  }
  
  ###############################
  # Check that number of trees in each ensemble
  # are positive integers
  ###############################
  if(!(is.integer(hyper$M_vec) & all(hyper$M_vec > 0))){
    message("[validate_flexBART_hyper]: supplied M_vec =")
    print(hyper$M_vec)
    stop("[validate_flexBART_hyper]: All entries in M_vec must be positive integers")
  }
  
  ###############################
  # Check positive vector valued hyperparameters
  ###############################
  pos_vecs <- c("alpha_vec", "beta_vec", "tau_vec")
  for(param in pos_vecs){
    if(any(hyper[[param]] <= 0)){
      message("[validate_flexBART_hyper]: supplied ", param, "=", hyper[[param]])
      stop(paste("[validate_flexBART_hyper]: All entries in", param, "must be positive"))
    }
  }
  
  ###############################
  # Check positive valued hyperparameters
  ###############################
  pos_params <- 
    c("sigest", "sigquant", "nu", "lambda", "a_u", "b_u")
  for(param in pos_params){
    if(any(hyper[[param]] <= 0)){
      message("[validate_flexBART_hyper]: supplied ", param, "=", hyper[[param]])
      stop(paste("[validate_flexBART_hyper]:", param, "must be positive"))
    }
  }
  
  ###############################
  # Check logical hyperparameters
  ###############################
  log_params <- 
    c("nest_v", "nest_c", "sparse")
  for(param in log_params){
    if(!is.logical(hyper[[param]])){
      message(paste("[validate_flexBART_hyper]: supplied ", param, "=", hyper[[param]]))
      stop(paste("[validate_flexBART_hyper]:", param, "must be logical"))
    }
  }
  
  ###############################
  # Check additional parameters
  ###############################
  
  if(!hyper$graph_cut_type %in% c(1L, 2L, 3L, 4L)){
    message(paste("[validate_flexBART_hyper]: supplied graph_cut_type =", hyper$graph_cut_type))
    stop("[validate_flexBART_hyper]: graph_cut_type must be 1L, 2L, 3L, or 4L")
  }
  
  if(!hyper$nest_v_option %in% c(0L, 1L, 2L, 3L)){
    message(paste("[validate_flexBART_hyper]: supplied graph_cut_type =", hyper$graph_cut_type))
    stop("[validate_flexBART_hyper]: graph_cut_type must be 0L, 1L, 2L, or 3L")
  }
  
}

###############################
# Actually parse the hyperparameters
###############################

parse_hyper <- function(R, y_range,...){
  ###############################
  # Create an empty flexBART_hypers object
  ###############################
  hyper <- new_flexBART_hyper()
  target_names <- names(hyper)
  ###############################
  # Grab any user-specified arguments
  ###############################
  usr_args <- list(...)
  usr_names <- names(usr_args) # equivalent to `...names()`
  
  ###############################
  # Check whether user supplied anything for M_vec
  ###############################
  tmp_M <- pmatch(usr_names, table = "M_vec", duplicates.ok = FALSE)
  if(all(is.na(tmp_M))){
    # user doesn't provide any arguments for M
    hyper$M_vec <- rep(50L, times = R)
  } else{
    ix <- which(!is.na(tmp_M))
    tmp_name <- usr_names[[ix]]
    usr_M <- usr_args[[tmp_name]]
    if(!is.integer(usr_M)) usr_M <- as.integer(usr_M)
    
    if(length(usr_M) == 1) hyper$M_vec <- rep(usr_M, times = R)
    else if(length(usr_M) == R) hyper$M_vec <- usr_M
    else{
      message(paste("[parse_hyper]: matched M_vec to argument", tmp_name, "="))
      print(usr_M)
      message(paste("[parse_hyper]:expected length", R, ", the number of ensembles"))
      stop(paste("[parse_hyper]: supplied argument must have length 1 or", R))
    }
    rm(ix, tmp_name, usr_M)
  }
  rm(tmp_M)
  ###############################
  # Check whether user supplied anything for alpha_vec
  ###############################
  tmp_alpha <- pmatch(usr_names, table = "alpha_vec", duplicates.ok = FALSE)
  if(all(is.na(tmp_alpha))){
    hyper$alpha_vec <- rep(0.95, times = R)
  } else{
    ix <- which(!is.na(tmp_alpha))
    tmp_name <- usr_names[[ix]]
    usr_alpha <- usr_args[[tmp_name]]
    if(length(usr_alpha) == 1) hyper$alpha_vec <- rep(usr_alpha, times = R)
    else if(length(usr_alpha) == R) hyper$alpha_vec <- usr_alpha
    else{
      message(paste("[parse_hyper]: matched alpha_vec to supplied argument", tmp_name))
      message("supplied argument =")
      print(usr_alpha)
      message(paste("[parse_hyper]:expected length", R, ", the number of ensembles"))
      stop(paste("[parse_hyper]: supplied argument must have length 1 or", R))
    }
    rm(ix, tmp_name, usr_alpha)
  }
  rm(tmp_alpha)
  ###############################
  # Check whether user supplied anything for beta_vec
  ###############################
  tmp_beta <- pmatch(usr_names, table = "beta_vec", duplicates.ok = FALSE)
  if(all(is.na(tmp_beta))){
    hyper$beta_vec <- rep(2, times = R)
  } else{
    ix <- which(!is.na(tmp_beta))
    tmp_name <- usr_names[[ix]]
    usr_beta <- usr_args[[tmp_name]]
    if(length(usr_beta) == 1) hyper$beta_vec <- rep(usr_beta, times = R)
    else if(length(usr_beta) == R) hyper$beta_vec <- usr_beta
    else{
      message(paste("[parse_hyper]: matched beta_vec to supplied argument", tmp_name))
      message("supplied argument =")
      print(usr_beta)
      message(paste("[parse_hyper]:expected length", R, ", the number of ensembles"))
      stop(paste("[parse_hyper]: supplied argument must have length 1 or", R))
    }
    rm(ix, tmp_name, usr_beta)
  }
  rm(tmp_beta)
  ###############################
  # Check whether user supplied anything for mu0_vec
  ###############################
  tmp_mu <- pmatch(usr_names, table = "mu0_vec", duplicates.ok = FALSE)
  if(all(is.na(tmp_mu))){
    hyper$mu0_vec <- rep(0, times = R)
  } else{
    ix <- which(!is.na(tmp_mu))
    tmp_name <- usr_names[[ix]]
    usr_mu <- usr_args[[tmp_name]]
    if(length(usr_mu) == 1) hyper$mu0_vec <- rep(usr_mu, times = R)
    else if(length(usr_mu) == R) hyper$mu_vec <- usr_mu
    else{
      message(paste("[parse_hyper]: matched mu0_vec to supplied argument", tmp_name))
      message("supplied argument =")
      print(usr_mu)
      message(paste("[parse_hyper]:expected length", R, ", the number of ensembles"))
      stop(paste("[parse_hyper]: supplied argument must have length 1 or", R))
    }
    rm(ix, tmp_name, usr_mu)
  }
  rm(tmp_mu)
  ###############################
  # Check whether user supplied anything for tau_vec
  ###############################
  tmp_tau <- pmatch(usr_names, table = "tau_vec", duplicates.ok = FALSE)
  if(all(is.na(tmp_tau))){
    if(R == 1){
      hyper$tau_vec <- y_range/(2 * 2 * sqrt(hyper$M_vec[1]))
    } else{
      hyper$tau_vec <- 0.5/sqrt(hyper$M_vec)
    }
  } else{
    ix <- which(!is.na(tmp_tau))
    tmp_name <- usr_names[[ix]]
    usr_tau <- usr_args[[which(is.na(tmp_tau))]]
    if(length(usr_tau) == 1) hyper$tau_vec <- rep(usr_tau, times = R)
    else if(length(usr_tau) == R) hyper$tau_vec <- usr_tau
    else{
      message(paste("[parse_hyper]: matched tau_vec to supplied argument", tmp_name))
      message("supplied argument =")
      print(usr_tau)
      message(paste("[parse_hyper]:expected length", R, ", the number of ensembles"))
      stop(paste("[parse_hyper]: supplied argument must have length 1 or", R))
    }
    rm(ix, tmp_name, usr_tau)
  }
  rm(tmp_tau)

  ###############################
  # Check whether user supplied anything for sigma prior
  ###############################
  if("nu" %in% usr_names) nu <- usr_args[["nu"]]
  else nu <- 3
  
  if(nu < 0){
    message(paste("[parse_hyper]: supplied nu =", nu))
    stop("[parse_hyper]: nu must be positive!")
  }
  
  if("sigest" %in% usr_names) hyper$sigest <- usr_args[["sigest"]]
  else hyper$sigest <- 1
  if(hyper$sigest < 0){
    message(paste0("[parse_hyper]: supplied sigest =", hyper$sigest))
    stop("[parse_hyper]: sig_est must be positive (and, ideally, less than 1; see Details)")
  }
  
  if("sigquant" %in% usr_names) hyper$sigquant <- usr_args[["sigquant"]]
  else hyper$sigquant <- 0.9
  if(abs(hyper$sigquant-0.5) > 0.5){
    message(paste0("[parse_hyper]: supplied sigquant =", hyper$sigquant))
    stop("[parse_hyper]: sigquant must be between 0 and 1")
  }

  if("lambda" %in% usr_names) hyper$lambda <- usr_args[["lambda"]]
  else hyper$lambda <- qchisq(1-hyper$sigquant, df = hyper$nu)/hyper$nu * hyper$sigest^2

  ###############################
  # Check whether user supplied anything else
  ###############################
  param_names <- 
    c("graph_cut_type", "nest_v", "nest_v_option", "nest_c",
      "sparse", "a_u", "b_u")
  for(param in param_names){
    if(param %in% usr_names){
      # override with user supplied value
      hyper[[param]] <- usr_args[[param]]
    }
  }
  ###############################
  # Validate the hyperparameters
  ###############################
  validate_flexBART_hyper(hyper)
  return(hyper)
}
