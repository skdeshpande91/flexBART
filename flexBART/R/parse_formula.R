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

  # extract variables within bart calls
  ensm_terms <-  sub("bart\\s*\\(([^)]*)\\)", "\\1", ensm_ix)
  # remove + between variables store each term as a vector
  ensm_terms <- strsplit(ensm_terms, "\\+")
  # clear any existing white space (as a precaution)
  ensm_terms <- lapply(ensm_terms, trimws) 

  # Extract additional intercepts
  z_names <- attr(terms(frmla), "term.labels")[!grepl('bart', attr(terms(frmla), "term.labels"))]
  
  ## handle '.' syntax 
  # first we identify all training terms that are not intercepts or the outcome
  X_names <- setdiff(vars, c(outcome_name, z_names))
  
  # ensm_terms will be a length R list where elements are predictors included in that bart
  ensm_terms <- lapply(1:length(ensm_terms), function(i){
    # if a bart has a . include all eligable predictors (X_names)
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
  
  # identify unique predictors (non-intercepts)
  covariate_names <- unique(unlist(ensm_terms))
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
  colnames(cov_ensm) <- rep(NA, R)
  if(length(z_names) > 0){
    colnames(cov_ensm)[((R-length(z_names))+1):ncol(cov_ensm)] <- z_names
  }
  
  for (i in 1:R) {
    cov_ensm[, i] <- rownames(cov_ensm) %in% ensm_terms[[i]]
  }
  
  ###############################
  # Critical that continuous variables precede categorical variables
  ###############################
  is_cat <- 
    sapply(train_data[,covariate_names], 
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
                       dimnames = list(c(cont_names, cat_names), c()))
  } else{
    cov_ensm <- cov_ensm[c(cont_names, cat_names),]
  }
  
  out <- list()
  out[["outcome_name"]] <- outcome_name
  out[["cov_ensm"]] <- cov_ensm
  return(out)
}
