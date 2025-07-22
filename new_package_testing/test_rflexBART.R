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
sourceCpp("../flexBART/src/rflexBART.cpp")

source("../flexBART/R/rflexBART.R")

source("generate_nested_example.R")

tmp_data <- unique(train_data[,c("Classroom", "School", "District")])

draw_nonest <-
  rflexBART(train_data = tmp_data,
            nd = 1000,
            verbose = TRUE,
            nest_v = FALSE, nest_c = FALSE)

draw_nest0 <- 
  rflexBART(train_data = tmp_data,
            nd = 10000,
            verbose = TRUE,
            nest_v = TRUE, nest_v_option = 0, nest_c = TRUE)

draw_nest1 <- 
  rflexBART(train_data = tmp_data,
            nd = 10000,
            verbose = TRUE,
            nest_v = TRUE, nest_v_option = 1, nest_c = TRUE)

draw_nest2 <- 
  rflexBART(train_data = tmp_data,
            nd = 10000,
            verbose = TRUE,
            nest_v = TRUE, nest_v_option = 2, nest_c = TRUE)

draw_nest3 <- 
  rflexBART(train_data = tmp_data,
            nd = 10000,
            verbose = TRUE,
            nest_v = TRUE, nest_v_option = 3, nest_c = TRUE)

range(draw_nonest$kernel - draw_nest0$kernel)
range(draw_nest0$kernel - draw_nest1$kernel)
range(draw_nest0$kernel - draw_nest2$kernel)
range(draw_nest0$kernel - draw_nest3$kernel)
range(draw_nest1$kernel - draw_nest2$kernel)



