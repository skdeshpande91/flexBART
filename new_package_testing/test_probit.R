mu_true <- function(df){
  tmp_X1 <- (1+df[,"X1"])/2
  tmp_X2 <- 1*(df[,"X2"]  == 1)
  return(0.5 * (3 * tmp_X1 + (2 - 5 * (tmp_X2)) * sin(pi * tmp_X1) - 2 * (tmp_X2)))
}

n_train <- 10000
n_test <- 1000
train_data <-
  data.frame(X1 = runif(n_train, min = -1, max = 1),
             X2 = rbinom(n_train, size = 1, prob = 0.5))
test_data <-
  data.frame(X1 = runif(n = n_test, min = -1, max = 1),
             X2 = rbinom(n_test, size = 1, prob = 0.5))

mu_train <- mu_true(train_data)
prob_train <- pnorm(mu_train)

mu_test <- mu_true(test_data)
prob_test <- pnorm(mu_test)

train_data[,"Y"] <- rbinom(n = n_train, size = 1, prob = prob_train)


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
source("../flexBART/R/probit_flexBART.R")

flex_fit <-
  probit_flexBART(formula = Y~bart(.),
                  train_data = train_data, 
                  test_data = test_data,
                  n.chains = 4, M = 50)

plot(prob_test, flex_fit$prob.test.mean, pch = 16, cex = 0.5,
     xlab = "Actual", ylab = "Predicted")
abline(a = 0, b = 1, col = 'red')

tmp_bart_train <- rep(0, times = n_train)
tmp_bart_test <- rep(0, times = n_test)
tmp_bart_time <- 0
for(c_ix in 1:4){
  bart_time <-
    system.time(
      bart_fit <- 
        BART::pbart(x.train = train_data[,colnames(train_data) != "Y"],
                    y.train = train_data[,"Y"],
                    x.test = test_data,
                    sparse = TRUE, ntree = 50,
                    ndpost = 1000, nskip = 1000))
  tmp_bart_train <- 
    tmp_bart_train + bart_fit$prob.train.mean
  tmp_bart_test <-
    tmp_bart_test <-
    tmp_bart_test + bart_fit$prob.test.mean
  tmp_bart_time <- tmp_bart_time + bart_time["elapsed"]
}
tmp_bart_train <- tmp_bart_train/4
tmp_bart_test <- tmp_bart_test/4


plot(prob_test, tmp_bart_test, pch = 16, cex = 0.5)

sqrt(mean( (prob_train - flex_fit$prob.train.mean)^2 ))
sqrt( mean( (prob_train - tmp_bart_train)^2 ))


sqrt(mean( (prob_test - flex_fit$prob.test.mean)^2 ))
sqrt( mean( (prob_test - tmp_bart_test)^2 ))
