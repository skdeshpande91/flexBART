bart_oracle <- function(Y_train, df_train, df_test, ix_list,
                        nd = 1000, burn = 1000, sparse = TRUE)
{
  n_fits <- length(ix_list)
  n_train <- nrow(df_train)
  n_test <- nrow(df_test)
  
  yhat.train <- matrix(nrow = nd, ncol = n_train)
  yhat.test <- matrix(nrow = nd, ncol = n_test)
  
  for(fit_ix in 1:n_fits){
    train_ix <- ix_list[[fit_ix]][["train"]]
    test_ix <- ix_list[[fit_ix]][["test"]]
    
    if(length(train_ix) == 11){
      fit <- BART::wbart(x.train = df_train[train_ix,1:10],
                         y.train = Y_train[train_ix],
                         x.test = df_test[test_ix,1:10],
                         sigest = 1,
                         ndpost = nd, nskip = burn, sparse = sparse,
                         printevery = 1 + nd + burn)
    } else{
      fit <- BART::wbart(x.train = df_train[train_ix,1:10],
                         y.train = Y_train[train_ix],
                         x.test = df_test[test_ix,1:10],
                         ndpost = nd, nskip = burn, sparse = sparse,
                         printevery = 1 + nd + burn)
    }
    yhat.train[,train_ix] <- fit$yhat.train
    yhat.test[,test_ix] <- fit$yhat.test
  }
  results <-
    list(yhat.train.mean = colMeans(yhat.train),
         yhat.train = yhat.train,
         yhat.test.mean = colMeans(yhat.test),
         yhat.test = yhat.test)
  return(results)
}
