# helper function to map levels of a categorical predictor to integers
get_categorical_mapping <- function(x, name){
  if(is.character(x)){
    message(paste("[get_categorical_mapping]:",name, "was a string and will be converted to a factor"))
    unik_vals <- sort(unique(x))
    raw_x <- factor(x, levels = unik_vals, labels = unik_vals)
  } else{
    raw_x <- x
  }
  
  if(!is.factor(raw_x)){
    stop(paste(name, "must be a factor"))
  }
  
  if(sum(is.na(raw_x)) > 0){
    # detected missing values
    missing_ix <- which(is.na(raw_x))
    cat("NA detected in", name, ". Creating a new level for missing values.\n")
    cat("If you wish to impute missing values, you should do so before calling flexBART.\n")
    if(any(levels(raw_x) == "NA_flexbart")){
      message(paste("Found value of 'NA_flexbart' in variable", name))
      # similar stop to 'b' in stan4bart
      stop("flexBART does not allow the value 'NA_flexbart' for categorical covariates")
    } else{
      levels(raw_x) <- c(levels(raw_x), "NA_flexbart")
      raw_x <- relevel(raw_x, ref = "NA_flexbart")
      raw_x[missing_ix] <- "NA_flexbart"
    } 
  }
  # if we don't have every level represented, this errors out
  #mapping <- data.frame(integer_coding = 0:(max(as.integer(raw_x)-1)), value = levels(raw_x))
  mapping <- 
    data.frame(integer_coding = as.integer(levels(raw_x))-1,
               value = levels(raw_x))
  mapping[which(mapping$value == "NA_flexbart"), "value"] <- NA
  return(mapping)
}

# helper function to convert categorical to integer value, given a mapping
convert_categorical <- function(x, name, mapping = NULL)
{
  if(is.null(mapping)){
    warning("[convert_categorical]: Parsing", name, " but not mapping provided. Mapping levels to integers now")
    mapping <- get_categorical_mapping(x,name)
  } else{
    if(!identical(colnames(mapping), c("integer_coding", "value")) | !is.data.frame(mapping)){
      stop("[convert_categorical]: mapping needs to be a 2-column data_frame with columns named 'value' and 'integer_coding'")
    }
    # check that all values of x are included in first column of mapping
    int_x <- rep(NA_integer_, times = length(x))
    
    if(!all(x %in% mapping$value)){
      message(paste("[convert_categorical]: Parsing categorical variable", name, ": detected new level"))
      stop("Found an unexpected categorical level. Consider re-leveling the factor in the training data")
    } else{
      for(k in 1:nrow(mapping)){
        if(is.na(mapping[k,"value"])){
          #missing is a categorical level
          index <- which(is.na(x))
        } else{
          index <- which(x == mapping[k,"value"])
        }
        if(length(index) > 0) int_x[index] <- mapping[k,"integer_coding"]
      } # closes loop over the levels found in mapping
    } # closes if/else checking whether there are new levels in x that aren't in mapping
    return(int_x)
  }
}
################################
# Check for nesting structure
################################
parse_nesting <- function(X_cat, dinfo)
{
  
  # we need to count the number of levels for each categorical variable
  K <- sapply(dinfo$cat_mapping_list, FUN = nrow)
  tmp <- .detect_nesting(X_cat, K, dinfo$p_cont)
  results <- list()
  if(!is.null(tmp$nesting)){
    results[["nest_list"]] <- tmp$nesting
  } else{
    results["nest_list"] <- list(NULL)
  }
  return(results)

  
}



###############################################
# Parse adjacency structure
###############################################
parse_adjacency <- function(adjacency_list, dinfo)
{
  edge_mat_list <- list()
  for(j in dinfo$cat_names){
    if(j %in% names(adjacency_list)){
      tmp_A <- adjacency_list[[j]]
      if(!identical(rownames(tmp_A), colnames(tmp_A))){
        message(paste("[parse_adjacency]: Dimension names for variable", j, " don't match"))
        stop("Row and column names for elements in adjacency_list must be identical")
      }
      if(any(is.na(dinfo$cat_mapping_list[[j]][,"value"]))){
        message(paste("[parse_adjacency]: Variable", j, "had missing values"))
        stop("Cannot have missing level in network-structured categorical variable")
      }
      if(!identical(sort(dinfo$cat_mapping_list[[j]][,"value"]), sort(colnames(tmp_A)))){
        message(paste("[parse_adjacency]: mismatch b/w detected values and dimnames of adjacency matrix for", j))
        stop("All possible values of the variable must appear in row and column names of adjacency matrix")
      }
      # re-order rows and columns of tmp_A to be the same as in mapping_list
      tmp_A <- tmp_A[dinfo$cat_mapping_list[[j]][,"value"], dinfo$cat_mapping_list[[j]][,"value"]]
      tmp_indices <- which(tmp_A !=0, arr.ind = TRUE)
      edge_mat <- tmp_indices[tmp_indices[,1] < tmp_indices[,2],] # get only the lower triangle of A
      colnames(edge_mat) <- c("from", "to")
      edge_mat_list[[j]] <- edge_mat-1 # remember C++ is 0-indexed
    } else{
      # no adjacency information for this variable
      # so put in a null value
      edge_mat_list[j] <- list(NULL)
    }
  } # close loop over the categorical variables
  return(edge_mat_list)
}

# Determines the max and min ranges of each truly continuous variable 
# This allows us to scale them appropriately
get_continuous_info <- function(x, name, pad = 0.1, n_unik_diffs = 5)
{
  n <- length(x)
  unik_x <- sort(unique(x))
  n_unik <- length(unik_x)
  consecutive_diffs <- unik_x[-1] - unik_x[-n_unik]
  
  if(length(unique(consecutive_diffs)) < n_unik_diffs){
    cat("[get_continuous info]:", name, " suspected to be discrete. Defining a grid of cutpoints for", name, "\n")
    cat("Using the unique values of x as splitting points.\n")
    cat("To use a different grid, manually set the `cutpoints_list` argument of flexBART.\n")
    x_sd <- NA
    x_min <- min(x)
    x_max <- max(x)
  } else{
    x_sd <- sd(x)
    x_min = min(x) - pad*x_sd
    x_max = max(x) + pad*x_sd
  }
  return(list(x_sd = x_sd, x_min = x_min, x_max = x_max))
}

# rescale truly continuous functions and leave discrete, ordinal ones alone
convert_continuous <- function(x, x_min, x_max, discrete = FALSE)
{
  if(discrete){
    std_x <- x
    cutpoints <- sort(unique(x))
  } else{
    std_x <- 2 * (x - x_min)/(x_max - x_min) - 1
    cutpoints <- NULL
  }
  return(list(std_x = std_x, cutpoints = cutpoints))
}





# Function used to parse data frame, figure out the continuous and categorical preds
# and then determine information used to re-scale and re-code these values during
# training and testing
get_covariate_info <- function(cov_data, pad = NULL, n_unik_diffs = NULL)
{
  stopifnot(is.data.frame(cov_data))
  is_cat <- 
    sapply(cov_data, 
           FUN = function(x){return( (is.factor(x) | is.character(x)))})
  # determine number of continuous and categorical covariates
  dinfo <- new_flexBART_data_info()
  dinfo$p_cont <- sum(1-is_cat)
  dinfo$p_cat <- sum(is_cat)
  dinfo$p <- dinfo$p_cont + dinfo$p_cat
  
  if(dinfo$p_cont == 0 & dinfo$p_cat == 0){
    stop("[new_flexBART_data]: no covariates detected!")
  }
  if(ncol(cov_data) != dinfo$p){
    # should *never* be thrown but just in case
    message(paste("[get_covariate_info]: detected", dinfo$p, "covariates but cov_data has", ncol(cov_data), "columns"))
    stop("Incorrect number of columns in cov_data!")
  }
  
  if(dinfo$p_cont > 0){
    # Save continuous variable names so they can be used in later functions
    dinfo$cont_names <- colnames(cov_data)[!is_cat]
    if(is.null(pad)){
      cat("[get_covariate_info]:  No padding for continuous variable ranges provided.\n")
      cat("By default, flexBART adds 0.1 SDs to min and max values of each continuous variable.\n")
      pad <- rep(0.1, times = dinfo$p_cont)
      names(pad) <- dinfo$cont_names
    } else{
      if(length(pad) != dinfo$p_cont){
        stop("[get_covariate_info]: pad must be a vector with length equal to the number of continuous predictors")
      }
      if(!identical(names(pad), dinfo$cont_names)){
        stop("[get_covariate_info]: pad must be a named vector. 
             names should match be the column names of cov_data corresponding 
             to continuous predictors")
      }
    }
    if(is.null(n_unik_diffs)){
      n_unik_diffs <- 5
    }

    
    dinfo$x_min <- rep(NA, times = dinfo$p_cont)
    dinfo$x_max <- rep(NA, times = dinfo$p_cont)
    dinfo$x_sd <- rep(NA, times = dinfo$p_cont)
    
    names(dinfo$x_min) <- dinfo$cont_names
    names(dinfo$x_max) <- dinfo$cont_names
    names(dinfo$x_sd) <- dinfo$cont_names
    
    for(j in dinfo$cont_names){
      tmp <- 
        get_continuous_info(x = cov_data[,j], name = j, 
                            pad = pad[j], n_unik_diffs = n_unik_diffs)
      dinfo$x_min[j] <- tmp$x_min
      dinfo$x_max[j] <- tmp$x_max
      dinfo$x_sd[j] <- tmp$x_sd
    } # closes loop getting information about continuous covariates
  }
  
  if(dinfo$p_cat > 0){
    dinfo$cat_names <- colnames(cov_data)[is_cat]
    # Build mapping of categorical variables
    dinfo$cat_mapping_list <- list()
    for(j in dinfo$cat_names){
      dinfo$cat_mapping_list[[j]] <- get_categorical_mapping(cov_data[,j], j)
    }
  }
  
  return(dinfo)
}