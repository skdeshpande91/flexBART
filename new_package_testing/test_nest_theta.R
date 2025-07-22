library(Rcpp)
library(RcppArmadillo)

sourceCpp("../flexBART/src/detect_nesting.cpp")
source("../flexBART/R/flexBART_structures.R")
source("../flexBART/R/parse_formula.R")
source("../flexBART/R/get_covariate_info.R")
source("../flexBART/R/get_z_col_id.R")
source("../flexBART/R/prepare_data.R")

sourceCpp("test_nest_theta.cpp")

######################
# Prepare dataset
######################
source("generate_nested_example.R")

frmla <- Y~bart(.)
tmp_form <- parse_formula(frmla, train_data)
outcome_name <- tmp_form$outcome_name
cov_ensm <- tmp_form$cov_ensm
tmp_data <- 
  prepare_data(train_data = train_data,
               outcome_name = outcome_name, 
               cov_ensm = cov_ensm, 
               test_data = NULL)


test_nest_theta(cov_ensm = cov_ensm,
                tX_cont_train = t(tmp_data$training_info$X_cont),
                tX_cat_train = t(tmp_data$training_info$X_cat),
                cutpoints_list = tmp_data$training_info$cutpoints,
                cat_levels_list = tmp_data$training_info$cat_levels_list,
                edge_mat_list = tmp_data$training_info$edge_mat_list,
                nest_list = tmp_data$training_info$nest_list)
############
frmla <- Y~bart(.)
tmp_form <- parse_formula(frmla, train_data2)
outcome_name <- tmp_form$outcome_name
cov_ensm <- tmp_form$cov_ensm
tmp_data2 <- 
  prepare_data(train_data = train_data2,
               outcome_name = outcome_name, 
               cov_ensm = cov_ensm, 
               test_data = NULL)


test_nest_theta(cov_ensm = cov_ensm,
                tX_cont_train = t(tmp_data2$training_info$X_cont),
                tX_cat_train = t(tmp_data2$training_info$X_cat),
                cutpoints_list = tmp_data2$training_info$cutpoints,
                cat_levels_list = tmp_data2$training_info$cat_levels_list,
                edge_mat_list = tmp_data2$training_info$edge_mat_list,
                nest_list = tmp_data2$training_info$nest_list)
