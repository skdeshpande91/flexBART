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
sourceCpp("../flexBART/src/predict_flexBART.cpp")

source("../flexBART/R/flexBART.R")
source("../flexBART/R/predict_flexBART.R")


source("generate_friedman_data.R")

flex_fit <-
  flexBART(formula = Y~bart(.),
           train_data = friedman_train,
           test_data = friedman_test,
           inform_sigma = TRUE, sparse = TRUE, 
           M = 200, 
           n.chains = 4)

yhat_train <-
  predict.flexBART(object = flex_fit,
                   newdata = friedman_train,
                   verbose = TRUE, print_every = 400)

range(yhat_train - flex_fit$yhat.train)

yhat_test <-
  predict.flexBART(object = flex_fit,
                   newdata = friedman_test,
                   verbose = TRUE, print_every = 400)
range(yhat_test - flex_fit$yhat.test)
