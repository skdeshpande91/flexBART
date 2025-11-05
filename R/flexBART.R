flexBART <- function(formula, 
                     train_data,
                     test_data = NULL,
                     initialize_sigma = TRUE, ...)
{
  ###############################
  # Capture additional arguments
  ###############################
  usr_args <- list(...)
  usr_names <- names(usr_args)

  ###############################
  # Parse family argument
  ###############################
  if(!("family" %in% usr_names)){
    # user didn't provide family and link arguments
    # default to standard BART
    family <- "gaussian"
    link <- "identity"
  } else if (inherits(usr_args[["family"]], "family")){
    family <- usr_args[["family"]][["family"]]
    link <- usr_args[["family"]][["link"]]
  } else {
    cat(paste("supplied family is class ", class(usr_args[["family"]]), " \n"))
    stop("family must be a family object")
  }

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
  heteroskedastic <- tmp_form$heteroskedastic
  if(heteroskedastic){
    if (!(family == "gaussian" && link == "identity")){
      stop("[flexBART]: heteroskedasticity is only supported for gaussian family with identity link")
    }
    cov_var <- tmp_form$cov_var
  }
  
  ###############################
  # Prepare the data to be passed to 
  # actual sampler
  ###############################
  tmp_data <- 
    prepare_data(train_data = train_data,
                 test_data = test_data,
                 outcome_name = outcome_name, 
                 cov_ensm = cov_ensm, 
                 ...)
  cat("Made it past prepare_data\n")
  # It will be useful to have problem dimensions readily accessible
  R <- tmp_data$training_info$R
  n_train <- length(tmp_data$training_info$std_Y)
  p_cont <- tmp_data$data_info$p_cont
  p_cat <- tmp_data$data_info$p_cat
  p <- tmp_data$data_info$p
  n_test <- 0
  offset <- tmp_data$training_info$offset
  if(length(tmp_data$testing_info$Z) > 1) n_test <- nrow(tmp_data$testing_info$Z)
  
  ###############################
  # Parse hyperparameters
  ###############################
  
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
    
    if("nest_v" %in% usr_names) nest_v <- usr_args[["nest_v"]]
    else{
      warning("[flexBART]: nesting structure detected but no nest_v argument specified. Defaulting to nest_v=TRUE\n")
      nest_v <- TRUE
    }
    
    if("nest_v_option" %in% usr_names) nest_v_option <- usr_args[["nest_v_option"]]
    else{
      if(nest_v) warning("[rflexBART]: nesting structure detected but no nest_v_option argument specified. Defaulting to nest_v_option=3\n")
      nest_v_option <- 3 # may need to change the default
    }
    
    if("nest_c" %in% usr_names) nest_c <- usr_args[["nest_c"]]
    else{
      if(nest_v){
        warning("[flexBART]: nesting structure detected & nest_v = TRUE but no nest_c argument specified. Default to nest_c = TRUE\n")
        nest_c <- TRUE
      } else{
        warning("[flexBART]: nesting structure detected but nest_v = FALSE and no nest_c argument specified. Default to nest_c = FALSE\n")
        nest_c <- FALSE
      }
    }
  }

  if (family == "gaussian" && link == "identity"){
    y_range <- max(tmp_data$training_info$std_Y) - min(tmp_data$training_info$std_Y)
    if("sigest" %in% usr_names){
      # user has supplied an initial estimate of sigma
      # check that it is less than 1; if it isn't, then we need to divide by sd(y)
      sigest <- usr_args[["sigest"]]
      if(sigest < 0){
        stop(paste("[flexBART]: supplied sigest =", sigest, ". Estimate of residual sd must be positive!"))
      }
      if(sigest > 1){
        cat(paste("supplied sigest = ", sigest, "greater than 1 \n"))
        message(paste("Internally, flexBART operates on standardized outcome scale. Dividing by outcome sd \n"))
        sigest <- sigest/tmp_data$training_info$y_sd
      }
    } else{
      if(initialize_sigma){
        if(p_cont == 1 & p_cat == 0){
          cat("no initial estimate of sigma provided. Initializing using OLS\n")
        } else{
          cat("no initial estimate of sigma provided. Initialize using LASSO\n")
        }
        sigest <- 
          get_sigma(tmp_data$training_info, tmp_data$data_info)
      } else{
        cat("no initial estimate of sigma provided.\n")
        sigest <- 1
      }
    }
    if (heteroskedastic){
      # stop("`sigma()` ensemble is not supported!")
      hyper <- 
        parse_hyper_heteroskedastic(R = R + 1,
                    y_range = y_range,
                    nest_v = nest_v, nest_v_option = nest_v_option, nest_c = nest_c, 
                    sigest = sigest, ...)
    } else{
      hyper <- 
        parse_hyper(R = R,
                    y_range = y_range,
                    nest_v = nest_v, nest_v_option = nest_v_option, nest_c = nest_c, 
                    sigest = sigest, ...)
    }
  } else if (family == "binomial" && link == "logit"){
    y_mean <- mean(tmp_data$training_info$std_Y)
    hyper <- 
      parse_hyper_logit(R = R, y_mean = y_mean,
                  nest_v = nest_v, 
                  nest_v_option = nest_v_option, 
                  nest_c = nest_c, 
                  ...)
  } else if (family == "binomial" && link == "probit"){
    y_mean <- mean(tmp_data$training_info$std_Y)
    hyper <- 
      parse_hyper_probit(R = R, y_mean = y_mean,
                  nest_v = nest_v, 
                  nest_v_option = nest_v_option, 
                  nest_c = nest_c, 
                  ...)
  # } else if (family == "poisson" && link == "log"){
  #   y_mean <- mean(tmp_data$training_info$std_Y)
  #   hyper <- 
  #     parse_hyper_poisson(R = R, y_mean = y_mean,
  #                 nest_v = nest_v, 
  #                 nest_v_option = nest_v_option, 
  #                 nest_c = nest_c, 
  #                 ...)
  } else {
    cat(paste("supplied family = ", family, " and link = ", link, " \n"))
    stop("Unsupported family and link combination!")
  }
  
  ###############################
  # Set control parameters
  ###############################  
  control <- parse_controls(...)
  
  if(control$verbose){
    if (family == "gaussian" && link == "identity") {
      cat("Initial sigma (after standardization) =", 
                  round(hyper$sigest, digits = 6), "\n")
    }
    if(!is.null(tmp_data$training_info$edge_mat_list)){
      cat("graph_cut_type = ", hyper$graph_cut_type, "\n")
    }
    if(!is.null(tmp_data$training_info$nest_list)){
      cat("nest_v = ", hyper$nest_v)
      cat(" nest_v_option = ", hyper$nest_v_option)
      cat(" nest_c = ", hyper$nest_c, "\n")
    }
    cat("n.chains = ", control$n.chains, "\n")
  }
  
  ###############################
  # Create containers for storing things
  ###############################  
  total_draws <- control$nd * control$thin + control$burn
  total_samples <- control$nd * control$n.chains
  
  if (family == "gaussian" && link == "identity") {
    if(heteroskedastic){
      # containers for sigma samples
      sigma_train_mean <- rep(0, times = n_train)
      if (control$save_samples) sigma_train_samples <- array(NA, dim = c(total_draws, n_train, control$n.chains))
      if (n_test > 0){
        sigma_test_mean <- rep(0, times = n_test)
        if (control$save_samples) sigma_test_samples <- array(NA, dim = c(total_draws - control$burn, n_test, control$n.chains))
      }
    } else{
      # Container for sigma samples:
      # all_sigma could be useful for assessing convergence
      # sigma_samples will get passed to predict to do posterior predictive sampling
      all_sigma <- array(NA, dim = c(total_draws, control$n.chains))
      sigma_samples <- rep(NA, times = total_samples)
    }
    # Containers for posterior mean of total fit & each beta
    yhat_train_mean <-rep(0, times = n_train)
    if(R > 1) raw_beta_train_mean <- array(0, dim = c(n_train, R))
    if(n_test > 0){
      yhat_test_mean <- rep(0, times = n_test)
      if(R > 1) raw_beta_test_mean <- array(0, dim = c(n_test, R))
    }
    # Containers for posterior samples
    if(control$save_samples){
      yhat_train_samples <- array(NA, dim = c(total_samples, n_train))
      if(R > 1) raw_beta_train_samples <- array(NA, dim =c(total_samples, n_train, R))
      if(n_test > 0){
        yhat_test_samples <- array(NA, dim = c(total_samples, n_test))
        if(R > 1) raw_beta_test_samples <- array(NA, dim = c(total_samples, n_test, R))
      }
    }
  } else if (family == "binomial") {
    # Containers for posterior mean of total fit & each beta
    prob_train_mean <-rep(0, times = n_train)
    if(R > 1) raw_beta_train_mean <- array(0, dim = c(n_train, R))
    if(n_test > 0){
      prob_test_mean <- rep(0, times = n_test)
      if(R > 1) raw_beta_test_mean <- array(0, dim = c(n_test, R))
    }
    # Containers for posterior samples
    if(control$save_samples){
      prob_train_samples <- array(NA, dim = c(total_samples, n_train))
      if(R > 1) raw_beta_train_samples <- array(NA, dim =c(total_samples, n_train, R))
      if(n_test > 0){
        prob_test_samples <- array(NA, dim = c(total_samples, n_test))
        if(R > 1) raw_beta_test_samples <- array(NA, dim = c(total_samples, n_test, R))
      }
    }
  } else if (family == "poisson" && link == "log"){
    # Containers for posterior mean of total fit & each beta
    yhat_train_mean <-rep(0, times = n_train)
    if(R > 1) raw_beta_train_mean <- array(0, dim = c(n_train, R))
    if(n_test > 0){
      yhat_test_mean <- rep(0, times = n_test)
      if(R > 1) raw_beta_test_mean <- array(0, dim = c(n_test, R))
    }
    # Containers for posterior samples
    if(control$save_samples){
      yhat_train_samples <- array(NA, dim = c(total_samples, n_train))
      if(R > 1) raw_beta_train_samples <- array(NA, dim =c(total_samples, n_train, R))
      if(n_test > 0){
        yhat_test_samples <- array(NA, dim = c(total_samples, n_test))
        if(R > 1) raw_beta_test_samples <- array(NA, dim = c(total_samples, n_test, R))
      }
    }
  } else {
    cat(paste("supplied family = ", family, " and link = ", link, " \n"))
    stop("Unsupported family and link combination!")
  }

  if (heteroskedastic){
  varcounts_samples <- 
    array(NA, dim = c(total_samples, p, R+1), 

          dimnames = list(c(), 
                          c(tmp_data$data_info$cont_names, 
                            tmp_data$data_info$cat_names), c()))
  } else{
    varcounts_samples <- 
    array(NA, dim = c(total_samples, p, R), 

          dimnames = list(c(), 
                          c(tmp_data$data_info$cont_names, 
                            tmp_data$data_info$cat_names), c()))
  }
  # container for timing
  timing <- rep(NA, times = control$n.chains)
  if(control$verbose){
    cat("n_train =", n_train, "n_test =", n_test, "\n")
    cat("R =", R, "p_cont =", p_cont, "p_cat =", p_cat, "\n")
    cat("Number of trees: ", hyper$M_vec, "\n")
    if (family == "gaussian" && link == "identity") {
      cat("Implied marginal priors:\n")
      for(r in 1:R){
        if(!is.na(colnames(cov_ensm))[r]){
          cat("  Effect of", colnames(cov_ensm)[r], ":")
        } else{
          cat("  Intercept: ")
        }
        cat("N(", round(hyper$mu0_vec[r]*hyper$M_vec[r], digits = 6), ", ", 
            round(hyper$tau_vec[r]^2*hyper$M_vec[r], digits = 6),
            ") \n", sep="")
      }
    }
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
    if (family == "gaussian" && link == "identity") {
      if (heteroskedastic){
        cat("Heteroskedasticity detected.\n")
        if (R == 1){
          tmp_time <- 
            system.time(
            fit <-
              ._single_fit_heteroskedastic(Y_train = tmp_data$training_info$std_Y,
                                          cov_ensm = cov_ensm,
                                          cov_var = cov_var,
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
                                          M_vec = hyper$M_vec,
                                          alpha_vec = hyper$alpha_vec,
                                          beta_vec = hyper$beta_vec,
                                          mu0_vec = hyper$mu0_vec,
                                          tau_vec = hyper$tau_vec,
                                          nd = control$nd, 
                                          burn = control$burn, 
                                          thin = control$thin,
                                          max_iter = hyper$max_iter,
                                          save_samples = control$save_samples, 
                                          save_trees = control$save_trees,
                                          verbose = control$verbose, 
                                          print_every = control$print_every)
            )
        } else{
          tmp_time <- 
            system.time(
            fit <-
              ._multi_fit_heteroskedastic(Y_train = tmp_data$training_info$std_Y,
                                          cov_ensm = cov_ensm,
                                          cov_var = cov_var,
                                          tZ_train = t(tmp_data$training_info$Z),
                                          tX_cont_train = t(tmp_data$training_info$X_cont),
                                          tX_cat_train = t(tmp_data$training_info$X_cat),
                                          tZ_test = t(tmp_data$testing_info$Z),
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
                                          M_vec = hyper$M_vec,
                                          alpha_vec = hyper$alpha_vec,
                                          beta_vec = hyper$beta_vec,
                                          mu0_vec = hyper$mu0_vec,
                                          tau_vec = hyper$tau_vec,
                                          nd = control$nd, 
                                          burn = control$burn, 
                                          thin = control$thin,
                                          max_iter = hyper$max_iter,
                                          save_samples = control$save_samples, 
                                          save_trees = control$save_trees,
                                          verbose = control$verbose, 
                                          print_every = control$print_every)
          )
          raw_beta_train_mean <- raw_beta_train_mean + fit$beta_train_mean/control$n.chains
          if(n_test > 0){
            raw_beta_test_mean <- 
              raw_beta_test_mean + fit$beta_test_mean/control$n.chains
          }
          if(control$save_samples){
            raw_beta_train_samples[start_index:end_index,,] <- fit$beta_train
            if(n_test > 0){
              raw_beta_test_samples[start_index:end_index,,] <- fit$beta_test
            }
          }
        }
        sigma_train_mean <- sigma_train_mean + fit$sigma_train_mean/control$n.chains
        if(control$save_samples) sigma_train_samples[,,chain_num] <- fit$sigma_train
        if(n_test > 0){
          sigma_test_mean <- sigma_test_mean + fit$sigma_test_mean/control$n.chains
          if(control$save_samples) sigma_test_samples[,,chain_num] <- fit$sigma_test
        }
        
        yhat_train_mean <- yhat_train_mean + fit$fit_train_mean/control$n.chains
        
        if(n_test > 0){
          yhat_test_mean <- 
            yhat_test_mean + fit$fit_test_mean/control$n.chains
        }
        if(control$save_samples){
          yhat_train_samples[start_index:end_index,] <- fit$fit_train
          if(n_test > 0){
            yhat_test_samples[start_index:end_index,] <- fit$fit_test
          }
        }
        varcounts_samples[start_index:end_index,,] <- fit$var_count
        if(control$save_trees){
          tree_list <- c(tree_list, fit$trees)
        }
      } else{
        if(R == 1){
          tmp_time <- 
            system.time(
              fit <- 
                ._single_fit(Y_train = tmp_data$training_info$std_Y,
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
                            sigest = hyper$sigest,
                            nu = hyper$nu,
                            lambda = hyper$lambda,
                            nd = control$nd, 
                            burn = control$burn, 
                            thin = control$thin,
                            save_samples = control$save_samples, 
                            save_trees = control$save_trees,
                            verbose = control$verbose, 
                            print_every = control$print_every))
          
        } else{
          tmp_time <-
            system.time(
              fit <- 
                ._multi_fit(Y_train = tmp_data$training_info$std_Y,
                            cov_ensm = cov_ensm,
                            tZ_train = t(tmp_data$training_info$Z),
                            tX_cont_train = t(tmp_data$training_info$X_cont),
                            tX_cat_train = t(tmp_data$training_info$X_cat),
                            tZ_test = t(tmp_data$testing_info$Z),
                            tX_cont_test = t(tmp_data$testing_info$X_cont),
                            tX_cat_test = t(tmp_data$testing_info$X_cat),
                            cutpoints_list = tmp_data$training_info$cutpoints,
                            cat_levels_list = tmp_data$training_info$cat_levels_list,
                            edge_mat_list = tmp_data$training_info$edge_mat_list,
                            nest_list = tmp_data$training_info$nest_list,
                            graph_cut_type = hyper$graph_cut_type,
                            sparse = hyper$sparse, 
                            a_u = hyper$a_u, b_u = hyper$b_u,
                            nest_v = hyper$nest_v,
                            nest_v_option = hyper$nest_v_option,
                            nest_c = hyper$nest_c,
                            M_vec = hyper$M_vec,
                            alpha_vec = hyper$alpha_vec, 
                            beta_vec = hyper$beta_vec,
                            mu0_vec = hyper$mu0_vec, 
                            tau_vec = hyper$tau_vec,
                            sigest = hyper$sigest,
                            nu = hyper$nu,lambda = hyper$lambda, 
                            nd = control$nd, 
                            burn = control$burn, 
                            thin = control$thin,
                            save_samples = control$save_samples, 
                            save_trees = control$save_trees,
                            verbose = control$verbose,
                            print_every = control$print_every))
          raw_beta_train_mean <- raw_beta_train_mean + fit$beta_train_mean/control$n.chains
          if(n_test > 0){
            raw_beta_test_mean <- 
              raw_beta_test_mean + fit$beta_test_mean/control$n.chains
          }
          if(control$save_samples){
            raw_beta_train_samples[start_index:end_index,,] <- fit$beta_train
            if(n_test > 0){
              raw_beta_test_samples[start_index:end_index,,] <- fit$beta_test
            }
          }
        }
        
        all_sigma[,chain_num] <- fit$sigma
        sigma_samples[start_index:end_index] <- fit$sigma[-(1:control$burn)]
        
        yhat_train_mean <- yhat_train_mean + fit$fit_train_mean/control$n.chains
        
        if(n_test > 0){
          yhat_test_mean <- 
            yhat_test_mean + fit$fit_test_mean/control$n.chains
        }
        if(control$save_samples){
          yhat_train_samples[start_index:end_index,] <- fit$fit_train
          if(n_test > 0){
            yhat_test_samples[start_index:end_index,] <- fit$fit_test
          }
        }
        varcounts_samples[start_index:end_index,,] <- fit$var_count
        if(control$save_trees){
          tree_list <- c(tree_list, fit$trees)
        }
      }
    } else if (family == "binomial"){
      if (link == "logit"){
        if(R == 1){
          tmp_time <-
            system.time(
              fit <-
                ._single_fit_logit(Y_train = tmp_data$training_info$std_Y,
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
                                  max_iter = hyper$max_iter,
                                  save_samples = control$save_samples, 
                                  save_trees = control$save_trees,
                                  verbose = control$verbose, 
                                  print_every = control$print_every))
        } 
        else{
          tmp_time <-
            system.time(
              fit <-
                ._multi_fit_logit(Y_train = tmp_data$training_info$std_Y,
                                  cov_ensm = cov_ensm,
                                  tZ_train = t(tmp_data$training_info$Z),
                                  tX_cont_train = t(tmp_data$training_info$X_cont),
                                  tX_cat_train = t(tmp_data$training_info$X_cat),
                                  tZ_test = t(tmp_data$testing_info$Z),
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
                                  M_vec = hyper$M_vec,
                                  alpha_vec = hyper$alpha_vec,
                                  beta_vec = hyper$beta_vec,
                                  mu0_vec = hyper$mu0_vec,
                                  tau_vec = hyper$tau_vec,
                                  nd = control$nd, 
                                  burn = control$burn, 
                                  thin = control$thin,
                                  max_iter = hyper$max_iter,
                                  save_samples = control$save_samples, 
                                  save_trees = control$save_trees,
                                  verbose = control$verbose, 
                                  print_every = control$print_every))
          raw_beta_train_mean <- raw_beta_train_mean + fit$beta_train_mean/control$n.chains
          if(n_test > 0){
            raw_beta_test_mean <- 
              raw_beta_test_mean + fit$beta_test_mean/control$n.chains
          }
          if(control$save_samples){
            raw_beta_train_samples[start_index:end_index,,] <- fit$beta_train
            if(n_test > 0){
              raw_beta_test_samples[start_index:end_index,,] <- fit$beta_test
            }
          }
        } # closes if/else checking how many ensembles there are
      } else if (link == "probit"){
        if(R == 1){
          tmp_time <-
            system.time(
              fit <-
                ._single_fit_probit(Y_train = tmp_data$training_info$std_Y,
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
        } else{
          tmp_time <-
            system.time(
              fit <-
                ._multi_fit_probit(Y_train = tmp_data$training_info$std_Y,
                                  cov_ensm = cov_ensm,
                                  tZ_train = t(tmp_data$training_info$Z),
                                  tX_cont_train = t(tmp_data$training_info$X_cont),
                                  tX_cat_train = t(tmp_data$training_info$X_cat),
                                  tZ_test = t(tmp_data$testing_info$Z),
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
                                  M_vec = hyper$M_vec,
                                  alpha_vec = hyper$alpha_vec,
                                  beta_vec = hyper$beta_vec,
                                  mu0_vec = hyper$mu0_vec,
                                  tau_vec = hyper$tau_vec,
                                  nd = control$nd, 
                                  burn = control$burn, 
                                  thin = control$thin,
                                  save_samples = control$save_samples, 
                                  save_trees = control$save_trees,
                                  verbose = control$verbose, 
                                  print_every = control$print_every))
          raw_beta_train_mean <- raw_beta_train_mean + fit$beta_train_mean/control$n.chains
          if(n_test > 0){
            raw_beta_test_mean <- 
              raw_beta_test_mean + fit$beta_test_mean/control$n.chains
          }
          if(control$save_samples){
            raw_beta_train_samples[start_index:end_index,,] <- fit$beta_train
            if(n_test > 0){
              raw_beta_test_samples[start_index:end_index,,] <- fit$beta_test
            }
          }
        } # closes if/else checking how many ensembles there are
      } # closes if/else checking binomial link
      prob_train_mean <- prob_train_mean + fit$fit_train_mean/control$n.chains
      
      if(n_test > 0){
        prob_test_mean <- 
          prob_test_mean + 
          fit$fit_test_mean/control$n.chains
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
    } else if (family == "poisson" && link == "log"){
      if(R == 1){
          tmp_time <-
            system.time(
              fit <-
                ._single_fit_poisson(Y_train = tmp_data$training_info$std_Y,
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
                                     max_iter = hyper$max_iter,
                                     save_samples = control$save_samples, 
                                     save_trees = control$save_trees,
                                     verbose = control$verbose, 
                                     print_every = control$print_every))
        } else{
          tmp_time <-
            system.time(
              fit <-
                ._multi_fit_poisson(Y_train = tmp_data$training_info$std_Y,
                                    cov_ensm = cov_ensm,
                                    tZ_train = t(tmp_data$training_info$Z),
                                    tX_cont_train = t(tmp_data$training_info$X_cont),
                                    tX_cat_train = t(tmp_data$training_info$X_cat),
                                    tZ_test = t(tmp_data$testing_info$Z),
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
                                    M_vec = hyper$M_vec,
                                    alpha_vec = hyper$alpha_vec,
                                    beta_vec = hyper$beta_vec,
                                    mu0_vec = hyper$mu0_vec,
                                    tau_vec = hyper$tau_vec,
                                    nd = control$nd, 
                                    burn = control$burn, 
                                    thin = control$thin,
                                    max_iter = hyper$max_iter,
                                    save_samples = control$save_samples, 
                                    save_trees = control$save_trees,
                                    verbose = control$verbose, 
                                    print_every = control$print_every))
          raw_beta_train_mean <- raw_beta_train_mean + fit$beta_train_mean/control$n.chains
          if(n_test > 0){
            raw_beta_test_mean <- 
              raw_beta_test_mean + fit$beta_test_mean/control$n.chains
          }
          if(control$save_samples){
            raw_beta_train_samples[start_index:end_index,,] <- fit$beta_train
            if(n_test > 0){
              raw_beta_test_samples[start_index:end_index,,] <- fit$beta_test
            }
          }
        }
        yhat_train_mean <- fit$fit_train_mean/control$n.chains
        if(n_test > 0){
          yhat_test_mean <- 
            yhat_test_mean + fit$fit_test_mean/control$n.chains
        }
        if(control$save_samples){
          yhat_train_samples[start_index:end_index,] <- fit$fit_train
          if(n_test > 0){
            yhat_test_samples[start_index:end_index,] <- fit$fit_test
          }
        }
        varcounts_samples[start_index:end_index,,] <- fit$var_count
        if(control$save_trees){
          tree_list <- c(tree_list, fit$trees)
        }
    } # closes if/else checking family
    timing[chain_num] <- tmp_time["elapsed"]
    if(control$verbose){
      cat("Ending chain", chain_num, "at", as.character(round(Sys.time())), "\n")
    }
  }
  ###############################
  # We have to rescale the posterior samples of beta
  # For notational compactness, will keep a copy of the relevant things
  ###############################
  y_mean <- tmp_data$training_info$y_mean
  y_sd <- tmp_data$training_info$y_sd
  z_mean <- tmp_data$training_info$z_mean
  z_sd <- tmp_data$training_info$z_sd
  z_col_id <- tmp_data$training_info$z_col_id
  
  
  if (family == "gaussian" & link == "identity") yhat_train_mean <- y_mean + y_sd * yhat_train_mean
  if (family == "poisson" & link == "log") yhat_train_mean <- yhat_train_mean * exp(offset)
  if(R > 1){
    beta_train_mean <- 
      rescale_beta_mean(raw_beta_train_mean, y_mean, y_sd, z_mean, z_sd, z_col_id)
  }
  
  if(n_test > 0){
    if (family == "gaussian" & link == "identity") yhat_test_mean <- y_mean + y_sd * yhat_test_mean
    if (family == "poisson" & link == "log") yhat_test_mean <- yhat_test_mean * exp(offset)
    if(R > 1){
      beta_test_mean <- 
        rescale_beta_mean(raw_beta_test_mean, y_mean, y_sd, z_mean, z_sd, z_col_id)
    }
  }
  if(control$save_samples){
    if (family == "gaussian" & link == "identity") yhat_train_samples <- y_mean + y_sd * yhat_train_samples
    if (family == "poisson" & link == "log") yhat_train_samples <- yhat_train_samples * exp(offset)
    if(R > 1){
      beta_train_samples <- 
        rescale_beta(raw_beta_train_samples, y_mean, y_sd, z_mean, z_sd, z_col_id)
    }
    if(n_test > 0){
      if (family == "gaussian" & link == "identity") yhat_test_samples <- y_mean + y_sd * yhat_test_samples
      if (family == "poisson" & link == "log") yhat_test_samples <- yhat_test_samples * exp(offset)
      if(R > 1){
        beta_test_samples <- 
          rescale_beta(raw_beta_test_samples, y_mean, y_sd, z_mean, z_sd, z_col_id)
      }
    }
  }
 
  results <- list()
  results[["dinfo"]] <- tmp_data$data_info
  if(control$save_trees) results[["trees"]] <- tree_list
  results[["scaling_info"]] <- 
    list(y_mean = y_mean, y_sd = y_sd,
         z_mean = z_mean, z_sd = z_sd,
         z_col_id = z_col_id, offset = offset)
  results[["M"]] <- hyper$M_vec
  results[["cov_ensm"]] <- cov_ensm
  
  results[["family"]] <- family
  results[["link"]] <- link

  if (family == "binomial") {
    results[["prob.train.mean"]] <- prob_train_mean
    if(R > 1){
      results[["beta.train.mean"]] <- beta_train_mean
      results[["raw_beta.train.mean"]] <- raw_beta_train_mean
    }
    if(control$save_samples){
      results[["prob.train"]] <- prob_train_samples
      if(R > 1){ 
        results[["beta.train"]] <- beta_train_samples
        results[["raw_beta.train"]] <- raw_beta_train_samples
      }
    }
    if(n_test > 0){
      results[["prob.test.mean"]] <- prob_test_mean
      if(R > 1){
        results[["beta.test.mean"]] <- beta_test_mean
        results[["raw_beta.test.mean"]] <- raw_beta_test_mean
      }
      if(control$save_samples){
        results[["prob.test"]] <- prob_test_samples
        if(R > 1){
          results[["beta.test"]] <- beta_test_samples
          results[["raw_beta.test"]] <- raw_beta_test_samples
        }
      }
    }
  } else if (family == "gaussian" & link == "identity") {
      results[["yhat.train.mean"]] <- yhat_train_mean
    if(R > 1){
      results[["beta.train.mean"]] <- beta_train_mean
      results[["raw_beta.train.mean"]] <- raw_beta_train_mean
    }
    if(n_test > 0){
      results[["yhat.test.mean"]] <- yhat_test_mean
      if(R > 1){
        results[["beta.test.mean"]] <- beta_test_mean
        results[["raw_beta.test.mean"]] <- raw_beta_test_mean
      }
    }
    
    if(control$save_samples){
      results[["yhat.train"]] <- yhat_train_samples
      if(R > 1){
        results[["beta.train"]] <- beta_train_samples
        results[["raw_beta.train"]] <- raw_beta_train_samples
      }
      if(n_test > 0){
        results[["yhat.test"]] <- yhat_test_samples
        if(R > 1){
          results[["beta.test"]] <- beta_test_samples
          results[["raw_beta.test"]] <- raw_beta_test_samples
        }
      }
    }

    if (heteroskedastic){
      results[["sigma.train.mean"]] <- sigma_train_mean * y_sd
      if(control$save_samples) results[["sigma.train"]] <- sigma_train_samples * y_sd
      if(n_test > 0){
        results[["sigma.test.mean"]] <- sigma_test_mean * y_sd
        if(control$save_samples) results[["sigma.test"]] <- sigma_test_samples * y_sd
      }
    } else{
      results[["sigma"]] <- sigma_samples * y_sd
    }
  } else if (family == "poisson" & link == "log") {
    results[["yhat.train.mean"]] <- yhat_train_mean
    if(R > 1){
      results[["beta.train.mean"]] <- beta_train_mean
      results[["raw_beta.train.mean"]] <- raw_beta_train_mean
    }
    if(n_test > 0){
      results[["yhat.test.mean"]] <- yhat_test_mean
      if(R > 1){
        results[["beta.test.mean"]] <- beta_test_mean
        results[["raw_beta.test.mean"]] <- raw_beta_test_mean
      }
    }
    
    if(control$save_samples){
      results[["yhat.train"]] <- yhat_train_samples
      if(R > 1){
        results[["beta.train"]] <- beta_train_samples
        results[["raw_beta.train"]] <- raw_beta_train_samples
      }
      if(n_test > 0){
        results[["yhat.test"]] <- yhat_test_samples
        if(R > 1){
          results[["beta.test"]] <- beta_test_samples
          results[["raw_beta.test"]] <- raw_beta_test_samples
        }
      }
    }
  }
  
  results[["varcounts"]] <- varcounts_samples
  results[["timing"]] <- timing
  class(results) <- c(class(results), "flexBART")
  return(results)
  
}
                     
    