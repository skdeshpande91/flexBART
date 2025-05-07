probit_flexBART <- function(formula, 
                            train_data,
                            test_data = NULL,...)
{
  
  ###############################
  # Capture additional arguments
  ###############################
  usr_args <- list(...)
  usr_names <- names(usr_args)
  
  ###############################
  # Parse the formula
  ###############################
  if(!is(formula, "formula")){
    frmla <- stats::formula(formula)
  } else{
    frmla <- formula
  }
  tmp_form <- parse_formula(frmla, train_data)
  outcome_name <- tmp_form$outcome_name
  cov_ensm <- tmp_form$cov_ensm
  
  if(ncol(cov_ensm) > 1) stop("probit_flexBART does not yet support multiple ensembles!") 
  
  ###############################
  # Prepare the data to be passed to 
  # actual sampler
  ###############################
  tmp_data <- 
    prepare_data(train_data = train_data,
                 outcome_name = outcome_name, 
                 cov_ensm = cov_ensm, 
                 test_data = test_data,
                 probit = TRUE,...)
  # It will be useful to have problem dimensions readily accessible
  R <- tmp_data$training_info$R
  n_train <- length(tmp_data$training_info$std_Y)
  p_cont <- tmp_data$data_info$p_cont
  p_cat <- tmp_data$data_info$p_cat
  p <- tmp_data$data_info$p
  n_test <- 0
  if(length(tmp_data$testing_info$Z) > 1) n_test <- nrow(tmp_data$testing_info$Z)
  
  
  ###############################
  # Parse hyperparameters
  ###############################
  y_mean <- mean(tmp_data$training_info$std_Y)
  if(is.null(tmp_data$training_info$nest_list)){
    # no nesting structure detected
    nest_v <- FALSE
    nest_v_option <- 3 # ignored
    nest_c <- FALSE
  } else{
    # we found nesting structure
    # if user didn't provide nest_v, nest_v_option, or nest_c, 
    # we need to use default value
    # should warn the user:
    if(!"nest_v" %in% usr_names){
      warning("[flexBART]: nesting structure detected but no nest_v argument specified. Defaulting to nest_v=TRUE")
      nest_v <- TRUE
    }
    if(! "nest_v_option" %in% usr_names){
      warning("[flexBART]: nesting structure detected but no nest_v argument specified. Defaulting to nest_v_option=3")
      nest_v_option <- 3
    }
    if(! "nest_c" %in% usr_names){
      warning("[flexBART]: nesting structure detected but no nest_c argument specified. Defaulting to nest_c=TRUE")
      nest_c <- TRUE
    }
  }
  
  hyper <- 
    parse_hyper_probit(R = R, y_mean = y_mean,
                nest_v = nest_v, nest_v_option, nest_c = nest_c, 
                ...)
  
  ###############################
  # Set control parameters
  ###############################  
  control <- parse_controls(...)
  
  if(control$verbose){
    if(!is.null(tmp_data$training_info$edge_mat_list)){
      cat("[flexBART]: graph_cut_type = ", hyper$graph_cut_type, "\n")
    }
    if(!is.null(tmp_data$training_info$nest_list)){
      cat("[flexBART]: nest_v = ", hyper$nest_v)
      cat(" nest_v_option = ", hyper$nest_v_option)
      cat(" nest_c = ", hyper$nest_c, "\n")
    }
  }
  
  ###############################
  # Create containers for storing things
  ###############################  
  total_draws <- control$nd * control$thin + control$burn
  total_samples <- control$nd * control$n.chains
  
  # Container for sigma samples:
  # all_sigma could be useful for assessing convergence
  # sigma_samples will get passed to predict to do posterior predictive sampling
  all_sigma <- array(NA, dim = c(total_draws, control$n.chains))
  sigma_samples <- rep(NA, times = total_samples)
  
  # Containers for posterior mean of total fit & each beta
  prob_train_mean <-rep(0, times = n_train)
  if(n_test > 0){
    prob_test_mean <- rep(0, times = n_test)
  }
  # Containers for posterior samples
  if(control$save_samples){
    prob_train_samples <- array(NA, dim = c(total_samples, n_train))
    if(n_test > 0){
      prob_test_samples <- array(NA, dim = c(total_samples, n_test))
    }
  }
  varcounts_samples <- 
    array(NA, dim = c(total_samples, p, R), 
          dimnames = list(c(), 
                          c(tmp_data$data_info$cont_names, 
                            tmp_data$data_info$cat_names), c()))
  # container for timing
  timing <- rep(NA, times = control$n.chains)
  if(control$verbose){
    cat("n_train =", n_train, "n_test =", n_test, "\n")
    cat("R =", R, "p_cont =", p_cont, "p_cat =", p_cat, "\n")
    cat("Number of trees: ", hyper$M_vec, "\n")
  }
  if(control$save_trees){
    tree_list <- list()
  }
  
  
  for(chain_num in 1:control$n.chains){
    if(control$verbose){
      cat("Starting chain", chain_num, "at", 
          as.character(round(Sys.time())), "\n")
    }
    start_index <- (chain_num-1)*control$nd + 1
    end_index <- chain_num*control$nd
    
    tmp_time <-
      system.time(
        fit <-
          .single_fit_probit(Y_train = tmp_data$training_info$std_Y,
                             cov_ensm = cov_ensm,
                             tX_cont_train = t(tmp_data$training_info$X_cont),
                             tX_cat_train = t(tmp_data$training_info$X_cat),
                             tX_cont_test = t(tmp_data$testing_info$X_cont),
                             tX_cat_test = t(tmp_data$testing_info$X_cat),
                             cutpoints_list = tmp_data$training_info$cutpoints,
                             cat_levels_list = tmp_data$training_info$cat_levels_list,
                             edge_mat_list = tmp_data$training_info$edge_mat_list,
                             nest_list = tmp_data$training_info$nest_list,
                             graph_cut_type = hyper$graph_cut_type,
                             sparse = hyper$sparse, 
                             a_u = hyper$a_u, 
                             b_u = hyper$b_u,
                             nest_v = hyper$nest_v,
                             nest_v_option = hyper$nest_v_option,
                             nest_c = hyper$nest_c,
                             M = hyper$M_vec[1],
                             alpha = hyper$alpha_vec[1],
                             beta = hyper$beta_vec[1],
                             mu0 = hyper$mu0_vec[1],
                             tau = hyper$tau_vec[1],
                             nd = control$nd, 
                             burn = control$burn, 
                             thin = control$thin,
                             save_samples = control$save_samples, 
                             save_trees = control$save_trees,
                             verbose = control$verbose, 
                             print_every = control$print_every))
    
    prob_train_mean <- prob_train_mean + fit$fit_train_mean/control$n.chains
    
    if(n_test > 0){
      prob_test_mean <- 
        prob_test_mean + fit$fit_test_mean/control$n.chains
    }
    if(control$save_samples){
      prob_train_samples[start_index:end_index,] <- fit$fit_train
      if(n_test > 0){
        prob_test_samples[start_index:end_index,] <- fit$fit_test
      }
    }
    varcounts_samples[start_index:end_index,,] <- fit$var_count
    if(control$save_trees){
      tree_list <- c(tree_list, fit$trees)
    }
    timing[chain_num] <- tmp_time["elapsed"]
    if(control$verbose){
      cat("Ending chain", chain_num, "at", as.character(round(Sys.time())), "\n")
    }
    
  }
  results <- list()
  results[["dinfo"]] <- tmp_data$data_info
  if(control$save_trees) results[["trees"]] <- tree_list
  results[["M"]] <- hyper$M_vec
  results[["cov_ensm"]] <- cov_ensm
  results[["is.probit"]] <- TRUE

  results[["prob.train.mean"]] <- prob_train_mean
  if(control$save_samples) results[["prob.train"]] <- prob_train_samples
  if(n_test > 0){
    results[["prob.test.mean"]] <- prob_test_mean
    if(control$save_samples) results[["prob.test"]] <- prob_test_samples
  }
  results[["varcounts"]] <- varcounts_samples
  results[["timing"]] <- timing
  if(control$save_trees) results[["trees"]] <- tree_list 
  class(results) <- c(class(results), "flexBART")
  return(results)
}
