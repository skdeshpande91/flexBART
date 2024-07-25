if(cat_unif){
  cat_probs <- rep(1/10, times = 10)
} else{
  # some relatively rare probabilities 
  cat_probs <- c(0.01, 0.1, 0.02, 0.2, 0.15, 0.03, 0.05, 0.15, 0.25, 0.04)
}
X_cont_train <- matrix(runif(n_all*10, min = -1, max = 1), 
                     nrow = n_all, ncol = 10,
                     dimnames = list(c(), paste0("X", 1:10)))
X_cat_train <- matrix(sample(0:9, size = n_all, replace = TRUE, prob = cat_probs), 
                    nrow = n_all, ncol = 1,
                    dimnames = list(c(), "X11"))

bart_df_train <- data.frame(X_cont_train,
                          X_cat_train)
bart_df_train[,"X11"] <- factor(X_cat_train[,1], levels = 0:9)

# For testing, create a set of 500 X's and use the same levels
X_cont_test <- matrix(nrow = 500 * 10, ncol = 10,
                      dimnames = list(c(), paste0("X", 1:10)))
X_cat_test <- matrix(rep(0:9, each = 500), ncol = 1,
                     dimnames = list(c(), "X11"))
tmp_X <- matrix(runif(500*10, min = -1, max = 1),
                nrow = 500, ncol = 10)
for(l in 0:9){
  X_cont_test[(1 + l*500):(500*(l+1)),] <- tmp_X
}
bart_df_test <- data.frame(X_cont_test, X_cat_test)
bart_df_test[,"X11"] <- factor(X_cat_test[,1], levels = 0:9)

# Generate mu

oracle_ix_list <- list()
if(dgp_ix == 1){
  mu_train <- dgp1(X_cont_train, X_cat_train,
                   offset_fried = 0.75,
                   offset_test = 15, scale = 2)
  
  mu_test <- dgp1(X_cont_test, X_cat_test,
                  offset_fried = 0.75,
                  offset_test = 15, scale = 2)
  
  oracle_ix_list[[1]] <- 
    list(train = sort(which(X_cat_train[,1] %in% c(0, 2, 4, 8))),
         test = sort(which(X_cat_test[,1] %in% c(0,2,4,8))))
  oracle_ix_list[[2]] <-
    list(train = sort(which(X_cat_train[,1] %in% c(1,3,5,6,7,9))),
         test = sort(which(X_cat_test[,1] %in% c(1,3,5,6,7,9))))
} else if(dgp_ix == 2){
  mu_train <- dgp2(X_cont_train, X_cat_train,
                   offset_fried = 0.75,
                   offset_test = 15, scale = 2)
  mu_test <- dgp2(X_cont_test, X_cat_test,
                  offset_fried = 0.75,
                  offset_test = 15, scale = 2)
  
  oracle_ix_list[[1]] <-
    list(train = sort(which(X_cat_train[,1] == 0)),
         test = sort(which(X_cat_test[,1] == 0)))
  oracle_ix_list[[2]] <-
    list(train = sort(which(X_cat_train[,1] != 0)),
         test = sort(which(X_cat_test[,1] != 0)))
} else if(dgp_ix == 3){
  mu_train <- dgp3(X_cont_train, X_cat_train,
                   offset_test = 15, scale_test = 2)
  mu_test <- dgp3(X_cont_test, X_cat_test, 
                  offset_test = 15, scale_test = 2)
  
  oracle_ix_list <- list()
  for(l in 0:6){
    oracle_ix_list[[l+1]] <-
      list(train = sort(which(X_cat_train[,1] == l)),
           test = sort(which(X_cat_test[,1] == l)))
  }
  oracle_ix_list[[8]] <-
    list(train = sort(which(X_cat_train[,1] %in% c(7,8,9))),
         test = sort(which(X_cat_test[,1] %in% c(7,8,9))))
} else if(dgp_ix == 4){
  mu_train <- dgp4(X_cont_train, X_cat_train,
                   offset_fried = 0.75,
                   offset_test = 15, scale = 2)
  mu_test <- dgp4(X_cont_test, X_cat_test,
                  offset_fried = 0.75,
                  offset_test = 15, scale = 2)
  oracle_ix_list <- list()
  for(l in 0:9){
    oracle_ix_list[[l+1]] <-
      list(train = sort(which(X_cat_train[,1] == l)),
           test = sort(which(X_cat_test[,1] == l)))
  }
}
Y_train <- mu_train + sigma * rnorm(n = n_all, mean = 0, sd = 1)

cat_levels_list <- list(0:9)
unif_cuts <- rep(TRUE, times = 10)
