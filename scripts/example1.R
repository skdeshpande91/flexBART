library(igraph)
library(scales)
library(Rcpp)
library(RcppArmadillo)
library(BART)
library(dbarts)
# currently it exports to R as a hidden function (with a . prefix)
sourceCpp("flexBART/src/flexBARTfast_fit.cpp") 
source("scripts/make_partition.R")

#n <- 50

#g <- erdos.renyi.game(n = n, p = 0.5, type = "gnp", directed = FALSE)
n <- 25
g <- erdos.renyi.game(n = n, p = 0.25, type = "gnp")
fixed_layout <- layout_nicely(g)
A <- as_adjacency_matrix(g, type = "both", sparse = FALSE)


partition1 <- make_partition(g, K = 2, min_size = 5)


beta0_true <- function(X_cont){
  scaled_X_cont <- (X_cont + 1)/2 # moves it to [0,1]
  z1 <- scaled_X_cont[,1]
  #z2 <- 1*(X_cont[,2] >= 0)
  z2 <- X_cont[,2]
  return(3 * z1 + (2 - 5 * z2) * sin(pi * z1) - 2 * z2)
}
beta1_true <- function(X_cont){
  scaled_X_cont <- (X_cont + 1)/2
  return( (3 - 3*cos(6*pi*scaled_X_cont[,1]) * scaled_X_cont[,1]^2) * (scaled_X_cont[,1] > 0.6) - (10 * sqrt(scaled_X_cont[,1])) * (scaled_X_cont[,1] < 0.25) )
}

mu_true <- function(X_cont, X_cat, partition1){
  # network structure in X_cat[,1]
  N <- nrow(X_cont)
  beta0 <- beta0_true(X_cont)
  beta1 <- beta1_true(X_cont)
  
  # partition1 has 2 pieces
  
  mu <- rep(NA, times = N)
  for(i in 1:N){
    if(partition1[X_cat[i,1]+1] == 1){
      mu[i] <- beta0[i]
    } else{
      mu[i] <- beta1[i]
    }
  }
  return(mu)
}

#T <- 50 # 50 observations per vertex
T <- 3
N <- n * T

X_cont <- matrix(runif(N*2, min = -1, max = 1), nrow = N, ncol = 2)
X_cat <- matrix(rep( (1:n) - 1, each = T), nrow = N, ncol = 1)
mu <- mu_true(X_cont, X_cat, partition1)

#sigma <- 2.5
sigma <- 1
Y_all <- mu + sigma * rnorm(n =N, mean = 0, sd = 1)
df_X <- data.frame(X_cont, X_cat)
df_X[,"X_cat"] <- factor(df_X[,"X_cat"])

X_bart <- makeModelMatrixFromDataFrame(df_X)


A_lower <- A
A_lower[upper.tri(A_lower)] <- 0
adj_support_list <- list(which(A_lower !=0) -1) # C++ is 0-indexed
cat_levels_list <- list(0:(n-1))


N_sim <- 25
rmse_train <- data.frame(adj = rep(NA, times = N_sim),
                         noadj = rep(NA, times = N_sim),
                         mybart = rep(NA, times = N_sim),
                         bart = rep(NA, times = N_sim))
rmse_test <- rmse_train
timing <- rmse_train
cov_train <- rmse_train
cov_test <- rmse_train

for(sim_number in 1:N_sim){
  # create a training & testing split
  set.seed(415+sim_number) # strictly for replicability purposes
  
  test_index <- sort(sample(1:N, size = floor(0.2 * N), replace = FALSE))
  train_index <- (1:N)[-test_index]
  
  Y_train <- Y_all[train_index]
  Y_test <- Y_all[test_index]
  
  X_cont_train <- X_cont[train_index,]
  X_cont_test <- X_cont[test_index,]
  
  X_cat_train <- matrix(X_cat[train_index,], ncol = 1)
  X_cat_test <- matrix(X_cat[test_index,], ncol = 1)
  
  X_bart_train <- X_bart[train_index,]
  X_bart_test <- X_bart[test_index,]
  
  mu_train <- mu[train_index]
  mu_test <- mu[test_index]
  
  M <- 200
  y_mean <- mean(Y_train)
  y_sd <- sd(Y_train)
  std_Y_train <- (Y_train - y_mean)/y_sd # standardize the output
  tau <- (max(std_Y_train) - min(std_Y_train))/(2 * 2 * sqrt(M)) # CGM10 prior sd on all leaf parameters
  
  # sigma^2 is given an Inv. Gamma(nu/2, nu * lambda/2) prior
  # we're using the standard BART choices for nu & lambda
  nu <- 3
  lambda <- qchisq(0.1, df = nu)/nu
  
  unif_cuts <- c(TRUE, TRUE)
  
  # Fit using adjacency structure
  adj_time <- system.time(
    adj_fit <- .flexBARTfast_fit(Y_train = std_Y_train,
                                 tX_cont_train = t(X_cont_train), # default value to signal no continuous preds
                                 tX_cat_train = t(X_cat_train),
                                 tX_cont_test =  t(X_cont_test),
                                 tX_cat_test = t(X_cat_test),
                                 unif_cuts = unif_cuts, 
                                 cutpoints_list = NULL,
                                 cat_levels_list = cat_levels_list,
                                 graph_split = rep(TRUE, times = 1), graph_cut_type = 0,
                                 adj_support_list = adj_support_list,
                                 rc_split = FALSE, prob_rc = 0.0, a_rc = 1, b_rc = 1,
                                 sparse = FALSE, a_u = 0.5, b_u = 1, 
                                 mu0 = 0, tau = tau,
                                 lambda = lambda, nu = nu,
                                 M = M, nd = 1000, burn = 1000, thin = 1,
                                 save_trees = FALSE, verbose = FALSE, print_every = 200))
  
  noadj_time <- system.time(
    noadj_fit <- .flexBARTfast_fit(Y_train = std_Y_train,
                                   tX_cont_train = t(X_cont_train), # default value to signal no continuous preds
                                   tX_cat_train = t(X_cat_train),
                                   tX_cont_test =  t(X_cont_test),
                                   tX_cat_test = t(X_cat_test),
                                   unif_cuts = unif_cuts, 
                                   cutpoints_list = NULL,
                                   cat_levels_list = cat_levels_list,
                                   graph_split = rep(FALSE, times = 1), graph_cut_type = 0,
                                   adj_support_list = adj_support_list,
                                   rc_split = FALSE, prob_rc = 0.0, a_rc = 1, b_rc = 1,
                                   sparse = FALSE, a_u = 0.5, b_u = 1, 
                                   mu0 = 0, tau = tau,
                                   lambda = lambda, nu = nu,
                                   M = M, 
                                   nd = 1000, burn = 1000, thin = 1,
                                   save_trees = FALSE, verbose = FALSE, print_every = 200))
  cutpoints_list <- list()
  for(j in 1:ncol(X_bart)) cutpoints_list[[j]] <- sort(unique(X_bart[,j]))
  
  mybart_time <- system.time(
    mybart_fit <- .flexBARTfast_fit(Y_train = std_Y_train,
                                    tX_cont_train = t(X_bart_train), # default value to signal no continuous preds
                                    tX_cat_train = matrix(0, nrow = 1, ncol = 1),
                                    tX_cont_test =  t(X_bart_test),
                                    tX_cat_test = matrix(0, nrow = 1, ncol = 1),
                                    unif_cuts = rep(FALSE, times = ncol(X_bart)), 
                                    cutpoints_list = cutpoints_list,
                                    cat_levels_list = cat_levels_list,
                                    graph_split = rep(FALSE, times = 1), graph_cut_type = 0,
                                    adj_support_list = adj_support_list,
                                    rc_split = FALSE, prob_rc = 0.0, a_rc = 1, b_rc = 1,
                                    sparse = FALSE, a_u = 0.5, b_u = 1, 
                                    mu0 = 0, tau = tau,
                                    lambda = lambda, nu = nu,
                                    M = M, 
                                    nd = 1000, burn = 1000, thin = 1,
                                    save_trees = FALSE, verbose = TRUE, print_every = 200))
  
  bart_time <- system.time(
    bart_fit <- wbart(x.train = X_bart_train, 
                      y.train = Y_train,
                      x.test = X_bart_test, 
                      ndpost = 1000, 
                      nskip = 1000, 
                      printevery = 2001,
                      nkeeptreedraws = 0))
  
  
  for(m in c("adj", "noadj", "mybart", "bart")){
    fit <- get(paste0(m, "_fit"))
    if(m != "bart"){
      mean_train <- y_mean + y_sd * apply(fit$fit_train, MARGIN = 2, FUN = mean)
      l95_train <- y_mean + y_sd * apply(fit$fit_train, MARGIN = 2, FUN = quantile, probs = 0.025)
      u95_train <- y_mean + y_sd * apply(fit$fit_train, MARGIN = 2, FUN = quantile, probs = 0.975)
      
      mean_test <- y_mean + y_sd * apply(fit$fit_test, MARGIN = 2, FUN = mean)
      l95_test <- y_mean + y_sd * apply(fit$fit_test, MARGIN = 2, FUN = quantile, probs = 0.025)
      u95_test <- y_mean + y_sd * apply(fit$fit_test, MARGIN = 2, FUN = quantile, probs = 0.975)
    } else{
      mean_train <- apply(fit$yhat.train, MARGIN = 2, FUN = mean)
      l95_train <- apply(fit$yhat.train, MARGIN = 2, FUN = quantile, probs = 0.025)
      u95_train <- apply(fit$yhat.train, MARGIN = 2, FUN = quantile, probs = 0.975)
      
      mean_test <- apply(fit$yhat.test, MARGIN = 2, FUN = mean)
      l95_test <- apply(fit$yhat.test, MARGIN = 2, FUN = quantile, probs = 0.025)
      u95_test <- apply(fit$yhat.test, MARGIN = 2, FUN = quantile, probs = 0.975)
    }
    
    
    rmse_train[sim_number, m] <- mean( (mu_train - mean_train)^2)
    cov_train[sim_number,m] <- mean( (mu_train >= l95_train) & (mu_train <= u95_train))
    
    rmse_test[sim_number, m] <- mean( (mu_test - mean_test)^2)
    cov_test[sim_number,m] <- mean( (mu_test >= l95_test) & (mu_test <= u95_test))
    
    timing[sim_number,m] <- get(paste0(m, "_time"))["elapsed"]
  }
  print(paste("Finished sim_number = ", sim_number))
  print(round(colMeans(rmse_test, na.rm = TRUE), digits = 3))
}

save(g, A, partition1, mu, X_cont, X_cat, X_bart, Y_all, sigma,
     rmse_train, rmse_test, cov_train, cov_test, timing,
     file = "data/example1.RData")


