################################################################################
# Testing VCBART with redundant ensembles
# Y = (beta0(x) + beta1(x)) + (beta2(x) + beta3(x)) * z1 + (beta4(x) + beta5(x)) * z2
################################################################################
# True functions
################################################################################
beta0_true <- function(df){
  tmp_X1 <- (1+df[,"X1"])/2
  tmp_X11 <- 1*(df[,"X11"] %in% c(0, 1, 4,5,6,8,9))
  return(3 * tmp_X1 + (2 - 5 * (tmp_X11)) * sin(pi * tmp_X1) - 2 * (tmp_X11))
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
  return(out)
}

beta2_true <- function(df){
  tmp_X <- (1+df[,"X1"])/2
  (3 - 3*cos(6*pi*tmp_X) * tmp_X^2) * (tmp_X > 0.6) - (10 * sqrt(tmp_X)) * (tmp_X < 0.25)
}

beta3_true <- function(df){
  tmp_X <- (1+df[,paste0("X", 1:p_cont)])/2
  return(exp(sin((0.9 * (tmp_X[,1] + 0.48))^10)) + tmp_X[,2] * tmp_X[,3] + tmp_X[,4])
}
beta4_true <- function(df){
  tmp_X <- (1+df[,paste0("X", 1:p_cont)])/2
  return(exp(3*tmp_X[,5])*cos(7*pi*tmp_X[,5]/2))
}
beta5_true <- function(df){
  tmp_X <- (1+df[,1:p_cont])/2
  return(10*sin(pi*tmp_X[,6] * tmp_X[,7]) + 
           20 * (tmp_X[,8] - 0.5)^2 + 
           10 * tmp_X[,9] + 
           5 * tmp_X[,10])
}

################################################################################
# Set problem dimensions
################################################################################
n_train <- 10000
n_test <- 1000
p_cont <- 10
p_cat <- 10
sigma <- 1

################################################################################
# Generate X & Z
################################################################################
train_data <- data.frame(Y = rep(NA, times = n_train))
test_data <- data.frame(Y = rep(NA, times = n_test))
for(j in 1:p_cont){
  train_data[,paste0("X",j)] <- 
    runif(n = n_train, min = -1, max = 1)
  test_data[,paste0("X",j)] <- 
    runif(n = n_test, min = -1, max = 1)
}
for(j in 1:p_cat){
  train_data[,paste0("X",j+p_cont)] <- 
    factor(sample(0:9, size = n_train, replace = TRUE), levels = 0:9)
  test_data[,paste0("X",j+p_cont)] <-
    factor(sample(0:9, size = n_test, replace = TRUE), levels = 0:9)
}

train_data[,"Z1"] <-
  rnorm(n = n_train, mean = 0, sd = 1)
test_data[,"Z1"] <-
  rnorm(n = n_test, mean = 0, sd = 1)

train_data[,"Z2"] <-
  rnorm(n = n_train, mean = 0, sd = 1)
test_data[,"Z2"] <-
  rnorm(n = n_test, mean = 0, sd = 1)


beta0_train <- beta0_true(train_data)
beta0_test <- beta0_true(test_data)

beta1_train <- beta1_true(train_data)
beta1_test <- beta1_true(test_data)

beta2_train <- beta2_true(train_data)
beta2_test <- beta2_true(test_data)

beta3_train <- beta3_true(train_data)
beta3_test <- beta3_true(test_data)

beta4_train <- beta4_true(train_data)
beta4_test <- beta4_true(test_data)

beta5_train <- beta5_true(train_data)
beta5_test <- beta5_true(test_data)



mu_train <-
  (beta0_train + beta1_train) +
  train_data[,"Z1"] * (beta2_train + beta3_train) + 
  train_data[,"Z2"] * (beta4_train + beta5_train)

mu_test <-
  (beta0_test + beta1_test) +
  test_data[,"Z1"] * (beta2_test + beta3_test) + 
  test_data[,"Z2"] * (beta4_test + beta5_test)

train_data[,"Y"] <-
  mu_train + sigma * rnorm(n = n_train, mean = 0, sd = 1)

test_data <- test_data[,colnames(test_data) != "Y"]

################################################################################
# Load all the functions for the package functions
################################################################################


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

sourceCpp("../flexBART/src/multi_ensm_fit.cpp")
sourceCpp("../flexBART/src/rescale_beta.cpp")

source("../flexBART/R/flexBART.R")

pkgfit <-
  flexBART(formula = Y ~ bart(.) + Z1 * bart(.) + Z1 * bart(.) + 
             Z2 * bart(.) + Z2 * bart(.),
           train_data = train_data,
           test_data = test_data,
           M = 50, n.chains = 4)

plot(beta0_test + beta1_test, pkgfit$beta.test.mean[,1],
     xlab = "Actual", ylab = "Predicted",
     pch = 16, cex = 0.5)

plot(beta2_test + beta3_test, pkgfit$beta.test.mean[,2],
     xlab = "Actual", ylab = "Predicted",
     pch = 16, cex = 0.5)

plot(beta4_test + beta5_test, pkgfit$beta.test.mean[,3],
     xlab = "Actual", ylab = "Predicted",
     pch = 16, cex = 0.5)


var_prob <- 
  apply(pkgfit$varcounts >= 1, MARGIN = c(2,3), FUN = mean)