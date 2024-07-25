targetBART <- function(Y_train, df_train, df_test,
                        p_cont, p_cat,
                        cat_levels_list,
                        nd = 1000, burn = 1000, sparse = FALSE){
  n_train <- nrow(df_train)
  n_test <- nrow(df_test)
  
  tmp_df_train <- data.frame(df_train[,1:p_cont])
  tmp_df_test <- data.frame(df_test[,1:p_cont])
  
  for(j in 1:p_cat){
    x_train <- rep(NA, times = n_train)
    x_test <- rep(NA, times = n_test)
    for(l in cat_levels_list[[j]]){
      ix_train <- which(df_train[,p_cont+j] == l)
      ix_test <- which(df_test[,p_cont+j] == l)
      if(length(ix_train) > 0){
        ybar <- mean(Y_train[ix_train])
        x_train[ix_train] <- ybar
        if(length(ix_test) > 0){
          x_test[ix_test] <- ybar
        }
      } else{
        if(length(ix_test)> 0){
          # no training data for this category
          x_test[ix_test] <- mean(Y_train)
        }
      }
    }
    tmp_df_train[[colnames(df_train)[p_cont+j]]] <- x_train
    tmp_df_test[[colnames(df_test)[p_cont+j]]] <- x_test
  }
  
  timing <-
    system.time(
      fit <- BART::pbart(x.train = tmp_df_train,
                         y.train = Y_train,
                         x.test = tmp_df_test,
                         sparse = sparse,
                         ndpost = nd, nskip = burn,
                         printevery = nd + burn + 1))["elapsed"]
  
  results <-
    list(prob.train.mean = fit$prob.train.mean,
         prob.test.mean = fit$prob.test.mean,
         timing = timing)
  return(results)
}
