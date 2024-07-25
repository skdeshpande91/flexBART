preprocess <- function(raw_data, 
                       cont_names = c(), 
                       cat_names = c(), out_name = c(), 
                       use_cut_names = c(),
                       log_y = FALSE)
{
  n_all <- nrow(raw_data)
  p_cont <- length(cont_names)
  p_cat <- length(cat_names)
  
  unif_cuts <- rep(TRUE, times = p_cont)
  names(unif_cuts) <- cont_names
  cutpoints_list <- list()
  cat_levels_list <- list()
  
  bart_df <- data.frame(tmp = rep(NA, times = n_all))
  for(j in cont_names){
    if(j %in% use_cut_names){
      # we specify cutpoints for this variable
      # and we do not rescale X
      tmp <- raw_data[,j]
      if(any(is.character(tmp))) stop("found a character")
      cutpoints_list[[j]] <- sort(unique(tmp))
      unif_cuts[j] <- FALSE
      bart_df[,j] <- tmp
    } else{
      tmp <- raw_data[,j]
      tmp_scaled <- scales::rescale(tmp, to = c(-1,1))
      cutpoints_list[[j]] <- c(0)
      bart_df[,j] <- tmp_scaled
    }
  }
  
  for(j in cat_names){
    tmp <- raw_data[,j]
    n_levels <- length(unique(tmp))
    cat_levels_list[[j]] <- 0:(n_levels-1)
    tmp_factor <- factor(tmp, labels = 0:(n_levels-1))
    bart_df[,j] <- tmp_factor
  }
  
  bart_df <- bart_df[,-1] # drop the placeholder
  
  bart_df <- bart_df[which(rowSums(is.na(bart_df)) == 0),]
  n_all <- nrow(bart_df)
  if(p_cont > 0){
    X_cont_all <- as.matrix(bart_df[,cont_names])
  } else{
    X_cont_all <- matrix(0, nrow = 1, ncol = 1)
  }
  
  if(p_cat > 0){
    X_cat_all <- matrix(nrow = n_all, ncol = p_cat, dimnames = list(c(), cat_names))
    for(j in cat_names){
      X_cat_all[,j] <- as.integer(bart_df[,j]) - 1
      cat_levels_list[[j]] <- sort(unique(as.integer(bart_df[,j]) - 1))
    }
  } else{
    X_cat_all <- matrix(0, nrow = 1, ncol  = 1)
    cat_levels_list <- list(c(0))
  }

  X_bart_all <- dbarts::makeModelMatrixFromDataFrame(x = bart_df, drop = FALSE)
  Y_all <- raw_data[,out_name[1]]
  if(log_y){
    if(any(Y_all == 0) && all(Y_all >= 0)) Y_all <- 1 + Y_all
    Y_all <- as.vector(log(Y_all))
  } 
  if(!is.vector(Y_all)) stop("something is weird")
  results <- list(X_cont_all = X_cont_all,
                  X_cat_all = X_cat_all,
                  bart_df_all = bart_df,
                  bart_mm_all = X_bart_all,
                  Y_all = Y_all,
                  n_all = n_all,
                  p_cont = p_cont,
                  p_cat = p_cat,
                  unif_cuts = unif_cuts,
                  cutpoints_list = cutpoints_list,
                  cat_levels_list = cat_levels_list)
  return(results)
}
