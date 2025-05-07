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

sourceCpp("../flexBART/src/single_ensm_fit.cpp")
sourceCpp("../flexBART/src/multi_ensm_fit.cpp")
sourceCpp("../flexBART/src/rescale_beta.cpp")

source("../flexBART/R/flexBART.R")


source("generate_friedman_data.R")

rmse_train <- c(flex = NA,  dbart = NA, bart = NA)
rmse_test <- c(flex = NA, dbart = NA, bart = NA)
timing <- c(flex = NA, dbart = NA, bart = NA)

flex_fit <-
  flexBART(formula = Y~bart(.),
           train_data = friedman_train,
           test_data = friedman_test,
           inform_sigma = TRUE, sparse = TRUE, 
           M = 200, 
           n.chains = 4)

rmse_train["flex"] <- sqrt( mean( (mu_train - flex_fit$yhat.train.mean)^2 ))
rmse_test["flex"] <- sqrt( mean( (mu_test - flex_fit$yhat.test.mean)^2 ))
timing["flex"] <- sum(flex_fit$timing)


###
# Run 4 chains of wbart

tmp_bart_train <- rep(0, times = n_train)
tmp_bart_test <- rep(0, times = n_test)
tmp_bart_time <- 0
for(c_ix in 1:4){
  bart_time <-
    system.time(
      bart_fit <- 
        BART::wbart(x.train = friedman_train[,colnames(friedman_train) != "Y"],
                    y.train = friedman_train[,"Y"],
                    x.test = friedman_test[,colnames(friedman_test) != "Y"],
                    sparse = TRUE,
                    ndpost = 1000, nskip = 1000))
  tmp_bart_train <- 
    tmp_bart_train + bart_fit$yhat.train.mean
  tmp_bart_test <-
    tmp_bart_test <-
    tmp_bart_test + bart_fit$yhat.test.mean
  tmp_bart_time <- tmp_bart_time + bart_time["elapsed"]
}
tmp_bart_train <- tmp_bart_train/4
tmp_bart_test <- tmp_bart_test/4

rmse_train["bart"] <- sqrt( mean( (mu_train - tmp_bart_train)^2 ))
rmse_test["bart"] <- sqrt( mean( (mu_test - tmp_bart_test)^2 ))
timing["bart"] <- tmp_bart_time


###
# Run 4 chains of dbart
dbart_time <-
  system.time(
    dbart_fit <- 
      dbarts::bart(x.train = friedman_train[,colnames(friedman_train) != "Y"],
                   y.train = friedman_train[,"Y"],
                   x.test = friedman_test[,colnames(friedman_test) != "Y"],
                   nchain = 4,
                   ndpost = 1000, nskip = 1000, keeptrees = TRUE))

rmse_train["dbart"] <- sqrt( mean( (mu_train - dbart_fit$yhat.train.mean)^2 ))
rmse_test["dbart"] <- sqrt( mean( (mu_test - dbart_fit$yhat.test.mean)^2 ))
timing["dbart"] <- dbart_time["elapsed"]


par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0))
plot(mu_test, flex_fit$yhat.test.mean,
     pch = 16, cex = 0.5, 
     xlab = "Actual", ylab = "Predicted")

var_probs <- 
  apply(flex_fit$varcounts >= 1, MAR = c(2,3), FUN = mean)
which(var_probs[,1] >= 0.5)