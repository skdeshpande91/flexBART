parse_formula <- function(frmla, train_data){
  
  vars <- colnames(train_data)
  ###############################
  # Get the outcome name & check that it's valid
  ###############################
  outcome_name <- all.vars(frmla)[attr(terms(frmla), "response")]
  if(!outcome_name %in% vars){
    stop(paste("[parse_formula]: supplied LHS in formula", outcome_name, "not found in training data"))
  }

  ###############################
  # R formula parser:
  # Pull out the terms and then get the number of ensembles
  # Indicate intercepts
  ###############################

  # extract formula terms 
  trms <- terms(frmla)

  # identify location of 'bart' calls in formula
  location <- gregexpr("bart\\s*\\(([^)]*)\\)", as.character(trms)[3], perl = TRUE)
  # extract individual bart calls
  ensm_ix <- regmatches(as.character(trms)[3], location)[[1]]
  # 'hard count' occurrences of bart calls 
  R <- length(ensm_ix)

  # identify location of 'sigma' calls in formula
  var_loc <- gregexpr("sigma\\s*\\(([^)]*)\\)", as.character(trms)[3], perl = TRUE)
  # extract individual sigma calls
  var_ensm_ix <- regmatches(as.character(trms)[3], var_loc)[[1]]
  # 'hard count' occurrences of sigma calls 
  if (length(var_ensm_ix) > 1){
    stop("[parse_formula]: multiple 'sigma' terms detected. Only one sigma term allowed!")
  } else if(length(var_ensm_ix) == 1){
    heteroskedastic <- TRUE
  } else{
    heteroskedastic <- FALSE
  }

  # extract variables within bart calls
  ensm_terms <-  sub("bart\\s*\\(([^)]*)\\)", "\\1", ensm_ix)
  # remove + between variables store each term as a vector
  ensm_terms <- strsplit(ensm_terms, "\\+")
  # clear any existing white space (as a precaution)
  ensm_terms <- lapply(ensm_terms, trimws)

  # extract variables within sigma calls
  if (heteroskedastic){
    var_terms <-  sub("sigma\\s*\\(([^)]*)\\)", "\\1", var_ensm_ix)
    # remove + between variables store each term as a vector
    var_terms <- strsplit(var_terms, "\\+")
    # clear any existing white space (as a precaution)
    var_terms <- lapply(var_terms, trimws)
  }

  # Extract additional intercepts

  # find *bart() or bart()* terms
  z_names_all_matches <- gregexpr("(?:([[:alnum:]_]+)\\*)?bart\\(([^)]*)\\)(?:\\*([[:alnum:]_]+))?",
                                  gsub(" ", "", as.character(trms)[3], fixed = TRUE), 
                                  perl = TRUE)

  # put in string fomrat
  z_names_matched_strings <- regmatches(gsub(" ", "", as.character(trms)[3], fixed = TRUE), z_names_all_matches)[[1]]

  # now extract relevant strings 
  z_names_capture_groups <- lapply(z_names_matched_strings, function(matched) {
    regmatches(matched, regexec("(?:([[:alnum:]_]+)\\*)?bart\\(([^)]*)\\)(?:\\*([[:alnum:]_]+))?", 
                                matched, 
                                perl = TRUE))[[1]]
  })

  # clean to only include z_name
  z_names <- do.call(c, lapply(z_names_capture_groups, function(x) {
    ifelse(x[[2]] == '', x[4], x[2])
  }))

  z_names[z_names == ''] <- NA_character_

  # check that there are no modifiers on the sigma term
  if (heteroskedastic){
    # find *sigma or sigma()* terms
    tmp_names_all_matches <- gregexpr("(?:([[:alnum:]_]+)\\*)?sigma\\(([^)]*)\\)(?:\\*([[:alnum:]_]+))?",
                                        gsub(" ", "", as.character(trms)[3], fixed = TRUE), 
                                        perl = TRUE)
    
    # put in string fomrat
    tmp_matched_strings <- regmatches(gsub(" ", "", as.character(trms)[3], fixed = TRUE), tmp_names_all_matches)[[1]]
    
    # now extract relevant strings 
    tmp_names_capture_groups <- lapply(tmp_matched_strings, function(matched) {
      regmatches(matched, regexec("(?:([[:alnum:]_]+)\\*)?sigma\\(([^)]*)\\)(?:\\*([[:alnum:]_]+))?", 
                                  matched, 
                                  perl = TRUE))[[1]]
    })
    
    # clean to only include z_name
    tmp_names <- do.call(c, lapply(tmp_names_capture_groups, function(x) {
      ifelse(x[[2]] == '', x[4], x[2])
    }))
    
    tmp_names[tmp_names == ''] <- NA_character_
    
    if(!is.na(tmp_names)){
      stop(paste("[parse_formula]: modifier", tmp_names[1], "detected on sigma term. Modifers are not allowed on sigma term!"))
    }
  }

  ## handle '.' syntax 
  # first we identify all training terms that are not intercepts or the outcome
  X_names <- setdiff(vars, c(outcome_name, z_names))

  # ensm_terms will be a length R list where elements are predictors included in that bart
  ensm_terms <- lapply(1:length(ensm_terms), function(i){
    # if a bart has a . include all eligible predictors (X_names)
    if(any(grepl('\\.', ensm_terms[[i]]))){
      collapsed_terms <- paste(ensm_terms[[i]], collapse = '+')
      # but remove any terms with a -
      remove_terms <-regmatches(collapsed_terms,gregexpr("-[^\\s,]+", collapsed_terms))[[1]]
      remove_terms <- gsub(' ', '', remove_terms)
      if(length(strsplit(remove_terms, split = '-'))> 0){
        remove_terms <- strsplit(remove_terms, split = '-')[[1]]
      }
      # some necessary book keeping
      ensm_terms[[i]] <- NULL
      ensm_terms[[i]] <- setdiff(X_names, remove_terms)
    }else{
      ensm_terms[[i]] <- ensm_terms[[i]]
    }
  })

  # var_terms will be a length 1 list where elements are predictors included in that bart
  if (heteroskedastic){
    var_terms <- lapply(1:length(var_terms), function(i){
      # if a bart has a . include all eligible predictors (X_names)
      if(any(grepl('\\.', var_terms[[i]]))){
        collapsed_terms <- paste(var_terms[[i]], collapse = '+')
        # but remove any terms with a -
        remove_terms <-regmatches(collapsed_terms,gregexpr("-[^\\s,]+", collapsed_terms))[[1]]
        remove_terms <- gsub(' ', '', remove_terms)
        if(length(strsplit(remove_terms, split = '-'))> 0){
          remove_terms <- strsplit(remove_terms, split = '-')[[1]]
        }
        # some necessary book keeping
        var_terms[[i]] <- NULL
        var_terms[[i]] <- setdiff(X_names, remove_terms)
      }else{
        var_terms[[i]] <- var_terms[[i]]
      }
    })
  }

  # identify unique predictors (non-intercepts)
  covariate_names <- unique(unlist(list(ensm_terms, var_terms)))
  p <- length(covariate_names)

  ###############################
  # Check that outcome is not a predictor
  ###############################
  if(outcome_name %in% covariate_names){
    stop(paste("[parse_formula]: supplied response variable", outcome_name,  "as a predictor"))
  }

  ###############################
  # Build the covariate ensemble matrix
  # Convention: for intercept terms, use a column name of NA in cov_ensm and Z
  # Otherwise, we need to check that the z variable name is in vars (i.e., it's in the training data)
  ###############################
  cov_ensm <- matrix(0, nrow = p, ncol = R)
  rownames(cov_ensm) <- covariate_names
  colnames(cov_ensm) <- rep(NA_character_, R)
  if(length(z_names) > 0){
    colnames(cov_ensm)[((R-length(z_names))+1):ncol(cov_ensm)] <- z_names
  }

  for(i in 1:R) {
    cov_ensm[, i] <- rownames(cov_ensm) %in% ensm_terms[[i]]
  }

  ###############################
  # Build the variance ensemble matrix
  # Convention: for intercept terms, use a column name of NA in cov_ensm and Z
  # Otherwise, we need to check that the z variable name is in vars (i.e., it's in the training data)
  ###############################
  if (heteroskedastic){
    cov_var <- matrix(0, nrow = p, ncol = 1)
    rownames(cov_var) <- covariate_names
    colnames(cov_var) <- rep(NA_character_, 1)
    cov_var[, 1] <- rownames(cov_var) %in% var_terms[[1]]
  } else{
    cov_var <- NA
  }

  ###############################
  # Critical that continuous variables precede categorical variables
  ###############################

  if(p == 1){
    tmp_df <- data.frame(train_data[,covariate_names[1]])
    colnames(tmp_df) <- covariate_names[1]
  } else{
    tmp_df <- train_data[,covariate_names]
  }

  is_cat <- 
    sapply(tmp_df, 
          FUN = function(x){return( (is.factor(x) | is.character(x)))})
  p_cont <- sum(1-is_cat)
  cont_names <- NULL
  if(p_cont > 0) cont_names <- covariate_names[!is_cat]

  p_cat <- sum(is_cat)
  cat_names <- NULL
  if(p_cat > 0) cat_names <- covariate_names[is_cat]

  if(R == 1){
    cov_ensm <- matrix(cov_ensm[c(cont_names, cat_names),], 
                      nrow = p, ncol = 1,
                      dimnames = list(c(cont_names, cat_names), colnames(cov_ensm)))
  } else{
    cov_ensm <- cov_ensm[c(cont_names, cat_names),]
  }

  if(heteroskedastic){
    cov_var <- matrix(cov_var[c(cont_names, cat_names),1], 
                      nrow = p, ncol = 1,
                      dimnames = list(c(cont_names, cat_names), colnames(cov_var)))
  }

  out <- list()
  out[["outcome_name"]] <- outcome_name
  out[["cov_ensm"]] <- cov_ensm
  out[["heteroskedastic"]] <- heteroskedastic
  out[["cov_var"]] <- cov_var
}
