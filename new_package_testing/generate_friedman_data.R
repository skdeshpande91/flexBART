################################################################################
# True function: Friedman function with p = 50 covariates
################################################################################
mu_true <- function(X_cont){
  # Recenter to [0,1]
  tmp_X <- (X_cont+1)/2
  return(10*sin(pi*tmp_X[,1] * tmp_X[,2]) + 
           20 * (tmp_X[,8] - 0.5)^2 + 
           10 * tmp_X[,17] + 
           5 * tmp_X[,20])
}
################################################################################
# Generate training & testing data
################################################################################
n_train <- 20000
n_test <- 1000
#p_cont <- 50
p_cont <- 25
p_cat <- 25
sigma <- 1

friedman_train <- data.frame(Y = rep(NA, times = n_train))
for(j in 1:p_cont) friedman_train[[paste0("X",j)]] <- runif(n_train, min = -1, max = 1)
for(j in 1:p_cat) friedman_train[[paste0("X",j+p_cont)]] <- 
  factor(sample(0:9, size = n_train, replace = TRUE), levels = 0:9)

mu_train <- mu_true(friedman_train[,paste0("X",1:p_cont)])
friedman_train[,"Y"] <- mu_train + sigma * rnorm(n = n_train, mean = 0, sd = 1)

friedman_test <- data.frame(Y = rep(NA, times = n_test))
for(j in 1:p_cont) friedman_test[[paste0("X",j)]] <- runif(n_test, min = -1, max = 1)
for(j in 1:p_cat) friedman_test[[paste0("X",j+p_cont)]] <- 
  factor(sample(0:9, size = n_test, replace = TRUE), levels = 0:9)

mu_test <- mu_true(friedman_test[,paste0("X",1:p_cont)])