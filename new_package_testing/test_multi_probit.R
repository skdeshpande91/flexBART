beta0_true <- function(df){
  tmp_X1 <- (1+df[,"X1"])/2
  tmp_X2 <- 1*(df[,"X2"]  == 1)
  return(0.5 * (3 * tmp_X1 + (2 - 5 * (tmp_X2)) * sin(pi * tmp_X1) - 2 * (tmp_X2)))
}
set.seed(517)
D <- 5000
omega0 <- rnorm(n = D, mean = 0, sd = 1.5)
b0 <- 2 * pi * runif(n = D, min = 0, max = 1)
weights <- rnorm(n = D, mean = 0, sd = 2*1/sqrt(D))
beta1_true <- function(df){
  # expect x to be in [-1,1]
  # but let's expand the interval to -5,5
  # z <- 10*x-5
  u <- 5*df[,"X1"]
  
  if(length(u) == 1){
    phi <- sqrt(2) * cos(b0 + omega0*u)
    out <- sum(weights * phi)
  } else{
    phi_mat <- matrix(nrow = length(u), ncol = D)
    for(d in 1:D){
      phi_mat[,d] <- sqrt(2) * cos(b0[d] + omega0[d] * u)
    }
    out <- phi_mat %*% weights
  }
  return(0.3 * out)
}

n_train <- 10000
n_test <- 1000
train_data <-
  data.frame(X1 = runif(n_train, min = -1, max = 1),
             X2 = rbinom(n_train, size = 1, prob = 0.5),
             Z = rbinom(n_train, size = 1, prob = 0.3))
test_data <-
  data.frame(X1 = runif(n = n_test, min = -1, max = 1),
             X2 = rbinom(n_test, size = 1, prob = 0.5),
             Z = rbinom(n_test, size = 1, prob = 0.3))

beta0_train <- beta0_true(train_data)
beta0_test <- beta0_true(test_data)

beta1_train <- beta1_true(train_data)
beta1_test <- beta1_true(test_data)

mu_train <-
  beta0_train + beta1_train * train_data$Z
mu_test <-
  beta0_test + beta1_test * test_data$Z

prob_train <- pnorm(mu_train, mean = 0, sd = 1)
prob_test <- pnorm(mu_test, mean = 0, sd = 1)

train_data$Y <- rbinom(n = n_train, size = 1, prob = prob_train)


library(Rcpp)
library(RcppArmadillo)
sourceCpp("../flexBART/src/detect_nesting.cpp")
source("../flexBART/R/get_z_col_id.R")
source("../flexBART/R/flexBART_structures.R")
source("../flexBART/R/get_covariate_info.R")
source("../flexBART/R/parse_formula.R")
source("../flexBART/R/prepare_data.R")
source("../flexBART/R/parse_controls.R")
source("../flexBART/R/parse_hyper.R")
source("../flexBART/R/get_sigma.R")

sourceCpp("../flexBART/src/single_ensm_probit.cpp")
sourceCpp("../flexBART/src/multi_ensm_probit.cpp")
sourceCpp("../flexBART/src/predict_flexBART.cpp")
sourceCpp("../flexBART/src/rescale_beta.cpp")
source("../flexBART/R/probit_flexBART.R")
source("../flexBART/R/predict_flexBART.R")

flex_fit <-
  probit_flexBART(formula = Y~bart(.)+Z*bart(.),
                  train_data = train_data, 
                  test_data = test_data,
                  n.chains = 4, M = 50)

plot(prob_train, flex_fit$prob.train.mean, pch = 16, cex = 0.5)
abline(a = 0, b = 1, col = 'red')
plot(prob_test, flex_fit$prob.test.mean, pch = 16, cex = 0.5)
abline(a = 0, b = 1, col = 'red')


plot(beta0_train, flex_fit$beta.train.mean[,1],pch = 16, cex = 0.5)
abline(a = 0, b = 1, col = 'red')

plot(beta1_train, flex_fit$beta.train.mean[,2], pch = 16, cex = 0.5)
abline(a = 0, b = 1, col = 'red')

plot(beta0_test, flex_fit$beta.test.mean[,1],pch = 16, cex = 0.5)
abline(a = 0, b = 1, col = 'red')
plot(beta1_test, flex_fit$beta.test.mean[,2], pch = 16, cex = 0.5)
abline(a = 0, b = 1, col = 'red')

