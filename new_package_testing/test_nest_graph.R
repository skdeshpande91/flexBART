library(Rcpp)
library(RcppArmadillo)

sourceCpp("../flexBART/src/detect_nesting.cpp")
source("../flexBART/R/get_z_col_id.R")
source("../flexBART/R/flexBART_structures.R")
source("../flexBART/R/get_covariate_info.R")
source("../flexBART/R/prepare_data.R")

source("generate_nested_example.R")

tmp_data <-
  prepare_data(train_data, outcome_name = "Y", cov_ensm = cov_ensm)

sourceCpp("test_nest_graph.cpp")

test_nest_graph(cat_levels_list = tmp_data$training_info$cat_levels_list,
                nest_list = tmp_data$training_info$nest_list,
                cov_ensm = cov_ensm,
                p_cont = tmp_data$data_info$p_cont,
                p_cat = tmp_data$data_info$p_cat)
