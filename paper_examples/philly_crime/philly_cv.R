library(flexBART)
library(BART)
library(dbarts)

load("philly_crime_data.RData")
load("sim_settings.RData")
source("helpers.R")

method <- "networkBART" # "BART-default", "BART-alt", "flexBART", or "networkBART"
sim_number <- 1 # should be between 1 and 100

set.seed(129+sim_number)
test_index <- sort(sample(1:N, size = floor(0.1*N), replace = FALSE))
train_index <- (1:N)[-test_index]

#test_index <- which(vertex_id_all == vertex)
#train_index <- (1:N)[-test_index]

X_cont_train <- matrix(X_cont_all[train_index,], ncol = 1)
X_cont_test <- matrix(X_cont_all[test_index,], ncol = 1)

X_cat_train <- matrix(X_cat_all[train_index,], ncol = 1)
X_cat_test <- matrix(X_cat_all[test_index,], ncol = 1)

vertex_id_train <- vertex_id_all[train_index]
vertex_id_test <- vertex_id_all[test_index]
Y_train <- Y_all[train_index]
Y_test <- Y_all[test_index]

tmp_mm <- dbarts::makeModelMatrixFromDataFrame(bart_df, drop = FALSE)
mm_train <- tmp_mm[train_index,]
mm_test <- tmp_mm[test_index,]

set.seed(724+sim_number)

if(method == "networkBART"){
  fit_time <- system.time(
    fit <- 
      flexBART::network_BART(Y_train = Y_train,
                             vertex_id_train = vertex_id_train,
                             X_cont_train = X_cont_train,
                             vertex_id_test = vertex_id_test,
                             X_cont_test = X_cont_test,
                             unif_cuts = unif_cuts,
                             cutpoints_list = cutpoints_list,
                             A = A_tract, 
                             graph_cut_type = 1,
                             save_trees = FALSE,
                             verbose = TRUE))
} else if(method == "flexBART"){
  fit_time <- system.time(
    fit <- 
      flexBART::flexBART(Y_train = Y_train,
                         X_cont_train = X_cont_train,
                         X_cat_train = X_cat_train,
                         X_cont_test = X_cont_test,
                         X_cat_test = X_cat_test,
                         unif_cuts = unif_cuts,
                         cutpoints_list = cutpoints_list,
                         cat_levels_list = cat_levels_list,
                         save_trees = FALSE,
                         verbose = TRUE))
} else if(method == "BART-alt"){
  fit_time <- system.time(
    fit <- 
      BART::wbart(x.train = mm_train,
                  y.train = Y_train, x.test = mm_test, 
                  sparse = FALSE, rm.const = FALSE,
                  ndpost = 1000, nskip = 1000))
} else if(method == "BART-default"){
  fit_time <- system.time(
    fit <- 
      BART::wbart(x.train = mm_train,
                  y.train = Y_train, x.test = mm_test, 
                  sparse = FALSE, rm.const = TRUE,
                  ndpost = 1000, nskip = 1000))
} else{
  stop("Method invalid")
}


sum_train <- summarize_post_pred(fit$yhat.train, fit$sigma[-(1:1000)])
sum_test <- summarize_post_pred(fit$yhat.test, fit$sigma[-(1:1000)])

assign(paste0(method, "_cv_", sim_number),
       list(method = method, sim_number = sim_number,
            rmse_train = sqrt(mean( (Y_train - sum_train[,"MEAN"])^2)),
            cov_train = coverage(Y_train, sum_train[,"L95"], sum_train[,"U95"]),
            int_train = interval_score(Y_train, sum_train[,"L95"], sum_train[,"U95"]),
            rmse_test = sqrt(mean( (Y_test - sum_test[,"MEAN"])^2)),
            cov_test = coverage(Y_test, sum_test[,"L95"], sum_test[,"U95"]),
            int_test = interval_score(Y_test, sum_test[,"L95"], sum_test[,"U95"]),
            timing = fit_time))
save(list = paste0(method, "_cv_", sim_number),
     file  = paste0(method, "_cv_", sim_number, ".RData"))


