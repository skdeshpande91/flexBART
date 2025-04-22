get_sigma <- function(trinfo){
 
  if(length(trinfo$X_cat) == 1){
    # only continuous covariates supplied
    l1_df <- data.frame(trinfo$X_cont)
  } else{
    if(length(trinfo$X_cont) > 1){
      # both categorical and continuous covariates supplied
      l1_df <- data.frame(trinfo$X_cont, trinfo$X_cat)
    } else{
      # only categorical covariates supplied
      l1_df <- data.frame(trinfo$X_cat)
    }
    for(j in dinfo$cat_names){
      l1_df[,j] <- factor(l1_df[,j], levels = dinfo$cat_mapping_list[[j]][,"integer_coding"])
    }
  }
  l1_X <- stats::model.matrix(~.-1, data = l1_df)
  
  l1_fit <- 
    glmnet::cv.glmnet(x = l1_X, y = trinfo$std_Y)
  fitted <- predict(object = l1_fit, newx = l1_X, s = "lambda.1se")
  return(sqrt(mean( (trinfo$std_Y - fitted)^2 )))

}