###############################
# Object to hold control parameters
###############################
new_flexBART_control <- function()
{
  out <- list()
  out["nd"] <- 1000L
  out["burn"] <- 1000L
  out["thin"] <- 1L
  out["n.chains"] <- 4L
  out["n.cores"] <- 1L
  out["save_samples"] <- TRUE
  out["save_trees"] <- TRUE
  out["verbose"] <- TRUE
  out["print_every"] <- as.integer(floor( (out$nd*out$thin + out$burn)/10))
  structure(out, class = "flexBART_control")
}

###############################
# Check that the control parameters are valid
###############################
validate_flexBART_control <- function(cntrl)
{
  exp_names <- 
    c("nd", "burn", "thin", 
      "n.chains", "n.cores",
      "save_samples", "save_trees", 
      "verbose", "print_every")
  if(!identical(names(cntrl), exp_names)){
    cat("[validate_flexBART_control]: supplied cntrl names:\n", names(cntrl), "\n")
    cat("[validate_flexBART_control]: expected cntrl names:\n", exp_names, "\n")
    stop("[validate_flexBART_control]: cntrl does not have valid names")
  }
  if(any(sapply(cntrl, FUN = is.null))){
    message("[validate_flexBART_control]: Null control parameters detected:")
    print(sapply(cntrl, FUN = is.null))
    stop("All control parameters must be non-null!")
  }
  ###############################
  # Check that positive integers parameters are positive & integers
  ###############################
  pos_params <- 
    c("nd", "burn", "thin", "n.chains", "n.cores", "print_every")
  for(param in pos_params){
    if(!(is.integer(cntrl[[param]]) & cntrl[[param]] > 0)){
      message(paste("[validate_flexBART_control]: supplied ", param, "=", cntrl[[param]]))
      stop(paste("[validate_flexBART_control]:", param, "must be a positive integer"))
    }
  }
  
  ###############################
  # Check logical parameters are logical
  ###############################
  log_params <- c("save_samples", "save_trees", "verbose")
  for(param in log_params){
    if(!is.logical(cntrl[[param]])){
      message(paste("[validate_flexBART_control]: supplied ", param, "=", cntrl[[param]]))
      stop(paste("[validate_flexBART_control]:", param, "must be a logical"))
    }
  }
}


parse_controls <- function(...){
  
  ###############################
  # Create an empty flexBART_control object
  ###############################
  cntrl <- new_flexBART_control()
  
  ###############################
  # Check if user passed any control parameters as arguments
  ###############################
  target_names <- names(cntrl)
    
  usr_args <- list(...)
  usr_names <- names(usr_args) # equivalent to `...names()`
  if(any(usr_names %in% target_names)){
    # user has passed some control parameters
    tmp_params <- usr_names[usr_names %in% target_names]
    for(params in tmp_params){
      if(params %in% c("save_samples", "save_trees", "verbose")) cntrl[[params]] <- usr_args[[params]]
      else cntrl[[params]] <- as.integer(usr_args[[params]])
    }
  }
  ###############################
  # Additional checks
  ###############################
  if(cntrl$thin > 1){
    message(paste0("[validate_flexBART_control]: supplied thin = ", cntrl$thin))
    warning("Thinning is not recommended. See documentation.")
  }
  
  if(cntrl$n.cores > 1){
    message(paste("[validate_flexBART_control]: supplied n.cores = ", cntrl$n.cores))
    warning("Parallel execution is not yet supported. Ignoring this argument.")
  }
  
  ###############################
  # Validate the controls
  ###############################
  validate_flexBART_control(cntrl)
  return(cntrl)
}
