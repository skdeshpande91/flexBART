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
sourceCpp("../flexBART/src/vcbart_fit.cpp")
sourceCpp("../flexBART/src/rescale_beta.cpp")
source("../flexBART/R/flexBART.R")


source("generate_friedman_data.R")
# Old
cov_ensm <- 
  matrix(1, nrow = p_cont, ncol = 1,
         dimnames = list(paste0("X",1:p_cont), c(NA_character_)))
tmp <-
  prepare_data(train_data = friedman_train,
               outcome_name = "Y",
               test_data = friedman_test,
               cov_ensm = cov_ensm)
flex_time <-
  system.time(
    flex_fit <- flexBART::flexBART(Y_train = friedman_train[,1],
                                   X_cont_train = tmp$training_info$X_cont,
                                   X_cont_test = tmp$testing_info$X_cont,
                                   sparse = TRUE,
                                   M = 200, nd = 1000, burn = 1000))


test0 <- 
  flexBART(formula = Y ~ bart(.),
           train_data = friedman_train,
           test_data = friedman_test,
           M = 200, inform_sigma = FALSE)


wbart_time <-
  system.time(
    wbart_fit <- BART::wbart(x.train = friedman_train[,-1],
                             y.train = friedman_train[,1],
                             x.test = friedman_test[,-1],
                             sparse = TRUE, ntree = 200, sigest = sd(friedman_train[,1]),
                             ndpost = 1000, nskip = 1000))
