# Test flexBART against plain vanilla BART
# in a setting with lots of categorical variables

library(Rcpp)
library(RcppArmadillo)
library(BART)

sourceCpp("flexBART/src/flexBART_fit.cpp")
sourceCpp("flexBART/src/flexBARTfast_fit.cpp")


#############
p_cont <- 10
p_cat <- 10

p <- p_cont + p_cat

mu_true <- function(X_cont, X_cat){
  scaled_X_cont <- (X_cont + 1)/2 # moves it to [0,1]
  tmp1 <- sin(pi * scaled_X_cont[,1] * scaled_X_cont[,2])
  tmp2 <- (scaled_X_cont[,3] - 0.5)^2
  tmp3 <- scaled_X_cont[,4]
  tmp4 <- scaled_X_cont[,5]
  

  #return((10 * (X_cat[,1] %in% 0:4))*tmp1 + (20 * (X_cat[,1] %in% c(0,2,4,6, 8)))*tmp2 + 
  #         (10*X_cat[,1] %in% c(1,3,5, 7,9))* tmp3 + (5 * (X_cat[,1] %in% c(1,2,3, 7,8,9)) * tmp4))
  
  
  return(10 * tmp1 + 
           20 * (X_cat[,1] %in% c(0, 2, 4,6,8)) * tmp2 +
           10 * (X_cat[,2] %in% c(1,2,3,6,7,8)) * tmp3 + 
           5 * (X_cat[,1] %in% c(1,2,3,4,5)) * tmp4)
}

###############
n <- 10000
set.seed(129)
X_cont <- matrix(runif(n * p_cont, -1, 1), nrow = n, ncol = p_cont)
X_cat <- matrix(NA, nrow = n, ncol = p_cat)

#X_cat[,1] <- as.integer(sample(0:50, size = n, replace = TRUE))
cat_levels_list <- list()
for(j in 1:p_cat){
  X_cat[,j] <- as.integer(sample(0:9), size = n, replace = TRUE)
  cat_levels_list[[j]] <- 0:9
}

#sigma <- 0.2
sigma <- 2.5
mu_all <- mu_true(X_cont, X_cat)
set.seed(129)
Y_all <- mu_all + sigma * rnorm(n, 0, 1)


rmse_test <- matrix(nrow = 25, ncol = 4, dimnames = list(c(), c("flexBART", "flexBARTfast", "BART", "DART")))
rmse_train <- rmse_test
timing <- rmse_test

for(r in 1:25){
  print(paste("Starting r =", r))
  test_index <- sort(sample(1:n, size = floor(0.2 * n), replace = FALSE))
  train_index <- (1:n)[-test_index]
  
  X_cont_train <- X_cont[train_index,]
  X_cat_train <- X_cat[train_index,]
  
  X_cont_test <- X_cont[test_index,]
  X_cat_test <- X_cat[test_index,]
  
  Y_train <- Y_all[train_index]
  Y_test <- Y_all[test_index]
  
  mu_train <- mu_all[train_index]
  mu_test <- mu_all[test_index]
  
  # Create some data frames to pass as input to BART::wbart
  df_X_train <- data.frame(X_cont_train, X_cat_train)
  df_X_test <- data.frame(X_cont_test, X_cat_test)
  
  colnames(df_X_train) <- c(paste0("X_cont", 1:p_cont), paste0("X_cat", 1:p_cat))
  colnames(df_X_test) <- c(paste0("X_cont", 1:p_cont), paste0("X_cat", 1:p_cat))
  for(j in 1:p_cat){
    df_X_train[,paste0("X_cat",j)] <- factor(X_cat_train[,j], levels = cat_levels_list[[j]])
    df_X_test[,paste0("X_cat",j)] <- factor(X_cat_test[,j], levels = cat_levels_list[[j]])
  }
  
  
  
  M <- 200
  y_mean <- mean(Y_train)
  y_sd <- sd(Y_train)
  std_Y_train <- (Y_train - y_mean)/y_sd # standardize the output
  tau <- (max(std_Y_train) - min(std_Y_train))/(2 * 2 * sqrt(M)) # CGM10 prior sd on all leaf parameters
  
  # sigma^2 is given an Inv. Gamma(nu/2, nu * lambda/2) prior
  # we're using the standard BART choices for nu & lambda
  nu <- 3
  lambda <- qchisq(0.1, df = nu)/nu
  
  
  cutpoints_list <- NULL
  adj_support_list <- NULL

  set.seed(22022022+r)
  #set.seed(31521+r)
  time1 <- system.time(fit1 <- .flexBART_fit(Y_train = std_Y_train,
                                             tX_cont_train = t(X_cont_train),
                                             tX_cat_train = t(X_cat_train),
                                             tX_cont_test = t(X_cont_test),
                                             tX_cat_test = t(X_cat_test),
                                             cutpoints_list = cutpoints_list,
                                             cat_levels_list = cat_levels_list,
                                             adj_support_list = adj_support_list,
                                             unif_cuts = TRUE,
                                             mst_split = FALSE, mst_reweight = FALSE,
                                             prob_aa = 0.5, prob_rc = 0, # prob_aa = 0 tells flexBART not to try to split on continuous variables
                                             mu0 = 0, tau = tau,
                                             lambda = lambda, nu = nu,
                                             M = M, 
                                             nd = 1000, burn = 1000, thin = 1,
                                             save_trees = FALSE, verbose = FALSE, print_every = 50))
  
  set.seed(22022022+r)
  time2 <- system.time(
    fit2 <- .flexBARTfast_fit(Y_train = std_Y_train,
                              tX_cont_train = t(X_cont_train),
                              tX_cat_train = t(X_cat_train),
                              tX_cont_test = t(X_cont_test),
                              tX_cat_test = t(X_cat_test),
                              cutpoints_list = cutpoints_list,
                              cat_levels_list = cat_levels_list,
                              adj_support_list = adj_support_list,
                              unif_cuts = TRUE,
                              mst_split = FALSE, mst_reweight = FALSE,
                              prob_aa = 0.5, prob_rc = 0, # prob_aa = 0 tells flexBART not to try to split on continuous variables
                              mu0 = 0, tau = tau,
                              lambda = lambda, nu = nu,
                              M = M, 
                              nd = 1000, burn = 1000, thin = 1,
                              save_trees = FALSE, verbose = FALSE, print_every = 50))
  
  # Using the same starting point & same random seed
  # flexBART and flexBARTfast should have identical results
  #identical(fit1$fit_train, fit2$fit_train)
  #identical(fit1$fit_test, fit2$fit_test)
  
  #range(fit1$fit_train - fit2$fit_train)
  #range(fit1$fit_test - fit2$fit_test)
  
  post_mean_train1 <- y_mean + y_sd * rowMeans(fit1$fit_train)
  post_mean_train2 <- y_mean + y_sd * rowMeans(fit2$fit_train)
  
  post_mean_test1 <- y_mean + y_sd * rowMeans(fit1$fit_test)
  post_mean_test2 <- y_mean + y_sd * rowMeans(fit2$fit_test)
  
  
  set.seed(22022022+r)
  
  bart1_time <- system.time(
    bart_fit1 <- wbart(x.train = df_X_train, y.train = Y_train,
                       x.test = df_X_test, sparse = FALSE, ndpost = 1000, nskip = 1000,
                       printevery = 5000))
  set.seed(22022022+r)
  bart2_time <- system.time(
    bart_fit2 <- wbart(x.train = df_X_train, y.train = Y_train,
                       x.test = df_X_test, sparse = TRUE, ndpost = 1000, nskip = 1000,
                       printevery = 5000))
  
  bart_mean_train1 <- bart_fit1$yhat.train.mean
  bart_mean_test1 <- bart_fit1$yhat.test.mean
  
  bart_mean_train2 <- bart_fit2$yhat.train.mean
  bart_mean_test2 <- bart_fit2$yhat.test.mean
  
  rmse_train[r, "flexBART"] <- sqrt(mean( (mu_train - post_mean_train1)^2 ))
  rmse_train[r, "flexBARTfast"] <- sqrt(mean( (mu_train - post_mean_train2)^2 ))
  rmse_train[r, "BART"] <- sqrt(mean( (mu_train - bart_mean_train1)^2 ))
  rmse_train[r, "DART"] <- sqrt(mean( (mu_train - bart_mean_train2)^2 ))
  
  rmse_test[r, "flexBART"] <- sqrt(mean( (mu_test - post_mean_test1)^2 ))
  rmse_test[r, "flexBARTfast"] <- sqrt(mean( (mu_test - post_mean_test2)^2 ))
  rmse_test[r, "BART"] <- sqrt(mean( (mu_test - bart_mean_test1)^2 ))
  rmse_test[r, "DART"] <- sqrt(mean( (mu_test - bart_mean_test2)^2 ))
  
  timing[r, "flexBART"] <- time1["elapsed"]
  timing[r, "flexBARTfast"] <- time2["elapsed"]
  timing[r, "BART"] <- bart1_time["elapsed"]
  timing[r, "DART"] <- bart2_time["elapsed"]
  
  print(round(colMeans(rmse_test, na.rm = TRUE), digits = 3))
}

save(rmse_train, rmse_test, timing, Y_all, X_cat, X_cont, mu_true, 
     mu_all, sigma, file = "results/timing_comparison.RData")
