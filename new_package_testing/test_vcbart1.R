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

mu_train <-
  beta0_train + train_data[,"Z1"] * beta1_train + train_data[,"Z2"] * beta2_train

mu_test <- 
  beta0_test + test_data[,"Z1"] * beta1_test + test_data[,"Z2"] * beta2_test

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


fit <-
  flexBART(formula = Y ~ bart(.) + Z1 * bart(.) + Z2 * bart(.),
           train_data = train_data,
           test_data = test_data,
           M_vec = c(50, 50, 50),
           inform_sigma = TRUE,
           sparse = TRUE,
           n.chains = 4)

plot(beta0_test, fit$beta.test.mean[,1], 
     pch = 16, cex = 0.5,
     xlab = "Actual", ylab = "Predicted")
plot(beta1_test, fit$beta.test.mean[,2], 
     pch = 16, cex = 0.5,
     xlab = "Actual", ylab = "Predicted")
plot(beta2_test, fit$beta.test.mean[,3], 
     pch = 16, cex = 0.5,
     xlab = "Actual", ylab = "Predicted")


var_probs <- 
  apply(fit$varcounts >= 1, MAR = c(2,3), FUN = mean)


frmla <- Y ~ bart(.) + Z1 * bart(.) + Z2 * bart(.)
tmp_form <- parse_formula(frmla, train_data)
outcome_name <- tmp_form$outcome_name
cov_ensm <- tmp_form$cov_ensm

tmp_data <- prepare_data(train_data = train_data,
                         outcome_name = outcome_name,
                         cov_ensm = cov_ensm,
                         test_data = test_data)

R <- ncol(cov_ensm)
y_range <- 
  max(tmp_data$training_info$std_Y) - min(tmp_data$training_info$std_Y)
sigest <- 
  get_sigma(tmp_data$training_info,
            tmp_data$data_info)
hyper <- 
  parse_hyper(R = R,
              y_range = y_range,
              sigest = sigest, M = 50, sparse = TRUE)


control <- parse_controls()

y_mean <- tmp_data$training_info$y_mean
y_sd <- tmp_data$training_info$y_sd
z_mean <- tmp_data$training_info$z_mean
z_sd <- tmp_data$training_info$z_sd
z_col_id <- tmp_data$training_info$z_col_id

multi_time <-
  system.time(
    multi_fit <- 
      ._multi_fit(Y_train = tmp_data$training_info$std_Y,
                  cov_ensm = cov_ensm,
                  tZ_train = t(tmp_data$training_info$Z),
                  tX_cont_train = t(tmp_data$training_info$X_cont),
                  tX_cat_train = t(tmp_data$training_info$X_cat),
                  tZ_test = t(tmp_data$testing_info$Z),
                  tX_cont_test = t(tmp_data$testing_info$X_cont),
                  tX_cat_test = t(tmp_data$testing_info$X_cat),
                  cutpoints_list = tmp_data$training_info$cutpoints,
                  cat_levels_list = tmp_data$training_info$cat_levels_list,
                  edge_mat_list = tmp_data$training_info$edge_mat_list,
                  nest_list = tmp_data$training_info$nest_list,
                  graph_cut_type = hyper$graph_cut_type,
                  sparse = hyper$sparse, 
                  a_u = hyper$a_u, b_u = hyper$b_u,
                  nest_v = hyper$nest_v,
                  nest_v_option = hyper$nest_v_option,
                  nest_c = hyper$nest_c,
                  M_vec = hyper$M_vec,
                  alpha_vec = hyper$alpha_vec, beta_vec = hyper$beta_vec,
                  mu0_vec = hyper$mu0_vec, tau_vec = hyper$tau_vec,
                  sigest = hyper$sigest,
                  nu = hyper$nu,lambda = hyper$lambda, 
                  nd = control$nd, burn = control$burn, thin = control$thin,
                  save_samples = control$save_samples, save_trees = control$save_trees,
                  verbose = control$verbose, print_every = control$print_every))


muhat_train <- y_mean + y_sd * multi_fit$fit_train_mean
muhat_test <- y_mean + y_sd * multi_fit$fit_test_mean
plot(mu_test, muhat_multi_test, pch = 16, cex = 0.5)


betahat_train <-
  rescale_beta_mean(multi_fit$beta_train_mean,
                    y_mean = y_mean, y_sd = y_sd,
                    x_mean = z_mean, x_sd = z_sd, 
                    z_col_id = z_col_id)

betahat_test <-
  rescale_beta_mean(multi_fit$beta_test_mean,
                    y_mean = y_mean, y_sd = y_sd,
                    x_mean = z_mean, x_sd = z_sd, 
                    z_col_id = z_col_id)

plot(beta0_test, betahat_test[,1], pch = 16, cex = 0.5)
plot(beta1_test, betahat_test[,2], pch = 16, cex = 0.5)
plot(beta2_test, betahat_test[,3], pch = 16, cex = 0.5)
