rflexBART <- function(train_data,
                      nd,
                      verbose = TRUE,
                      print_every = floor(nd/10),
                      ...)
{
  ###############################
  # Capture additional arguments
  ###############################
  usr_args <- list(...)
  usr_names <- names(usr_args)
  
  #cat(usr_names, "\n")
  
  
  n <- nrow(train_data)
  if("Y_flexBART" %in% colnames(train_data)){
    stop("Y_flexBART is a protected variable name. Re-name this covariate")
  }
  train_data <- cbind(Y_flexBART = rep(0, times = n), train_data)
  frmla <- Y_flexBART~bart(.)
  tmp_form <- parse_formula(frmla, train_data)
  outcome_name <- tmp_form$outcome_name
  cov_ensm <- tmp_form$cov_ensm
  
  tmp_data <- 
    prepare_data(train_data = train_data,
                 outcome_name = outcome_name, 
                 cov_ensm = cov_ensm,...)
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
      warning("[rflexBART]: nesting structure detected but no nest_v argument specified. Defaulting to nest_v=TRUE")
      nest_v <- TRUE
    }
    
    if("nest_v_option" %in% usr_names) nest_v_option <- usr_args[["nest_v_option"]]
    else{
      if(nest_v) warning("[rflexBART]: nesting structure detected but no nest_v_option argument specified. Defaulting to nest_v_option=3")
      nest_v_option <- 3 # may need to change the default
    }

    if("nest_c" %in% usr_names) nest_c <- usr_args[["nest_c"]]
    else{
      warning("[rflexBART]: nesting structure detected but no nest_c argument specified. Defaulting to nest_c=TRUE")
      nest_c <- TRUE
    }
  }
  hyper <- 
    parse_hyper(R = 1,
                y_range = 0,
                nest_v = nest_v, 
                nest_v_option = nest_v_option, 
                nest_c = nest_c,
                tau_vec = c(1), mu0_vec = c(0),...)
  if(verbose){
    if(!is.null(tmp_data$training_info$edge_mat_list)){
      cat("[flexBART]: graph_cut_type = ", hyper$graph_cut_type, "\n")
    }
    if(!is.null(tmp_data$training_info$nest_list)){
      cat("[flexBART]: nest_v = ", hyper$nest_v)
      cat(" nest_v_option = ", hyper$nest_v_option)
      cat(" nest_c = ", hyper$nest_c, "\n")
    }
  }

  tmp <- 
    ._draw_tree(tX_cont = t(tmp_data$training_info$X_cont),
                tX_cat = t(tmp_data$training_info$X_cat),
                cov_ensm = cov_ensm,
                cutpoints_list = tmp_data$training_info$cutpoints,
                cat_levels_list = tmp_data$training_info$cat_levels_list,
                edge_mat_list = tmp_data$training_info$edge_mat_list,
                nest_list = tmp_data$training_info$nest_list,
                graph_cut_type = hyper$graph_cut_type,
                nest_v = hyper$nest_v,
                nest_v_option = hyper$nest_v_option,
                nest_c = hyper$nest_c,
                alpha = hyper$alpha_vec[1], 
                beta = hyper$beta_vec[1],
                nd = nd, 
                verbose = verbose,
                print_every = print_every)
  colnames(tmp$varcounts) <- 
    c(tmp_data$data_info$cont_names, tmp_data$data_info$cat_names)
  return(tmp)
}
