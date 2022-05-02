




flexBART <- function(Y_train,
                     X_cont_train,
                     X_cat_train = matrix(1L, nrow = 1, ncol = 1),
                     unif_cuts = NULL,
                     cutpoints = NULL,
                     graph_split = NULL,
                     adjacency = NULL,
                     
                     sparse = TRUE, 
                     M = NULL, nd = 1000, burn = 1000, thin = 1,
                     save_trees = FALSE, verbose = TRUE, print_every = floor(M/10))
{
  
  if(!is.matrix(X_cont_train)) stop("X_cont_train must be a matrix!")
  if(!is.matrix(X_cat_train)) stop("X_cat_train must be a matrix!")
  
  p_cont <- 0
  p_cat <- 0
  
  if(nrow(X_cont_train) > 1){
    # we are actually passing continuous predictors
    # make sure that everything is scaled to 0/1. 
    # if not, then we need to linearly rescale everything
    p_cont <- ncol(X_cont_train)
    if(is.null(unif_cuts)){
      print("  continuous predictors supplied but unif_cuts argument not provided!")
      print("    falling back to default: all cutpoints drawn uniformly and all continuous predictors re-scaled to [-1,1]")
    } else{
      if(!is.logical(unif_cuts)) stop("unif_cuts must be a logical vector!")
      
      if(length(unif_cuts) != p_cont){
        print(paste("p_cont = ", p_cont, "continuous predictors supplied."))
        print(paste("length of supplied unif_cuts = ", length(unif_cuts)))
        stop("unif_cuts must have length equal to ncol(X_cont_train)")
      }
      
      if(! all(unif_cuts)){
        # we want to have cutpoints
        if(is.null(cutpoints)){
          print("You indicated that for some predictors you did not want to draw cutpoints from a continuous distribution")
          print("However, you did not supply a valid vector of cutpoints for those predictors")
          stop("invalid cutpoints argument")
        }
      }
      
      # if we have reached this point, we have a valid set of unif_cuts and cutpoints
      # check the scaling
      for(j in 1:p_cont){
        if(unif_cuts[j]){
          if(!all(abs(X_cont[,j]) <= 1)){
            print(paste("No cutpoints provided for continuous predictor", j))
            print("Expecting the predictor to be rescaled to interval [-1,1]!")
            stop("Continuous predictor without pre-defined cutpoints takes values outside [-1,1]")
          }
        }
      }
    }
    
    
  }
  
  if(nrow(X_cat_train) > 1){
    # we are actually passing categorical predictors.
    # need to check that (1) everything is an integer, (2) that cat_levels_list is not null,
    # (3) graph_cuts is a vector of booleans of the right size (will default to false), 
    # (4) adjacency is list of right size (will default to list of size p_cat with single entry 0 in each element)
    
  }
  
  
  y_mean <- mean(Y_train)
  y_sd <- sd(Y_train)
  std_Y_train <- (Y_train - y_mean)/y_sd # standardize the output
  tau <- (max(std_Y_train) - min(std_Y_train))/(2 * 2 * sqrt(M)) # CGM10 prior sd on all leaf parameters
  
  # sigma^2 is given an Inv. Gamma(nu/2, nu * lambda/2) prior
  # we're using the standard BART choices for nu & lambda
  nu <- 3
  lambda <- qchisq(0.1, df = nu)/nu
  
  
  fit  <- .flexBARTfast_fit(Y_train = std_Y_train,
                            tX_cont_train = t(X_cont_train), # default value to signal no continuous preds
                            tX_cat_train = t(X_cat_train),
                            tX_cont_test =  t(X_cont_test),
                            tX_cat_test = t(X_cat_test),
                            unif_cuts = unif_cuts, 
                            cutpoints_list = cutpoints_list,
                            cat_levels_list = cat_levels_list,
                            graph_split = graph_split, 
                            graph_cut_type = 0,
                            adj_support_list = adj_support_list,
                            rc_split = FALSE, prob_rc = 0.0, a_rc = 1, b_rc = 1,
                            sparse = sparse, a_u = 0.5, b_u = 1, 
                            mu0 = 0, tau = tau,
                            lambda = lambda, nu = nu,
                            M = M, 
                            nd = nd, burn = burn, thin = thin,
                            save_trees = save_trees, verbose = verbose, print_every = print_every)
  
  
  
  
}