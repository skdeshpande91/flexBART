bart_target <- function(Y_train, df_train, df_test,
                        nd = 1000, burn = 1000, sparse = TRUE){
  n_train <- nrow(df_train)
  n_test <- nrow(df_test)
  tmp_df_train <- df_train[,1:10]
  tmp_df_train[,"X11_new"] <- rep(NA, times = n_train)
  
  tmp_df_test <- df_test[,1:10]
  tmp_df_test[,"X11_new"] <- rep(NA, times = n_test)
  
  for(l in 0:9){
    ix_train <- which(df_train[,"X11"] == l)
    ix_test <- which(df_test[,"X11"] == l)
    if(length(ix_train) > 0){
      tmp_df_train[ix_train,"X11_new"] <- mean(Y_train[ix_train])
      if(length(ix_test) > 0){
        tmp_df_test[ix_test, "X11_new"] <- mean(Y_train[ix_train])
      }
    } else{
      if(length(ix_test > 0)){
        # when there is no data in this category, resort to overall mean
        tmp_df_test[ix_test,"X11_new"] <- mean(Y_train)
      }
    }
  }
  
  timing <-
    system.time(
      fit <- BART::wbart(x.train = tmp_df_train,
                         y.train = Y_train,
                         x.test = tmp_df_test,
                         sparse = sparse,
                         ndpost = nd, nskip = burn,
                         printevery = nd + burn + 1))["elapsed"]
  
  results <-
    list(yhat.train.mean = fit$yhat.train.mean,
         yhat.train = fit$yhat.train,
         yhat.test.mean = fit$yhat.test.mean,
         yhat.test = fit$yhat.test,
         sigma = fit$sigma,
         timing = timing)
  return(results)
}
