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
source("../flexBART/R/flexBART.R")

source("generate_grid_data.R")
adjacency_list <- list(X1 = A)

rmse_train <- 
  c(noadj = NA, gs1 = NA, gs2 = NA, gs3 = NA, gs4 = NA, dbart = NA, bart = NA)
rmse_test <- 
  c(noadj = NA, gs1 = NA, gs2 = NA, gs3 = NA, gs4 = NA, dbart = NA, bart = NA)
timing <- 
  c(noadj = NA, gs1 = NA, gs2 = NA, gs3 = NA, gs4 = NA, dbart = NA, bart = NA)




fit_gs1 <-
  flexBART(formula = Y ~ bart(X1),
           train_data = train_data,
           test_data = test_data,
           inform_sigma = TRUE,
           adjacency_list = adjacency_list,
           graph_cut_type = 1,
           n.chains = 4)

rmse_train["gs1"] <-
  sqrt(mean( (fit_gs1$yhat.test.mean[train_vertices] - mu[train_vertices])^2 ))
rmse_test["gs1"] <-
  sqrt(mean( (fit_gs1$yhat.test.mean[test_vertices] - mu[test_vertices])^2 ))
timing["gs1"] <- sum(fit_gs1$timing)

fit_gs2 <-
  flexBART(formula = Y ~ bart(X1),
           train_data = train_data,
           test_data = test_data,
           inform_sigma = TRUE,
           adjacency_list = adjacency_list,
           graph_cut_type = 2,
           n.chains = 4)
rmse_train["gs2"] <-
  sqrt(mean( (fit_gs2$yhat.test.mean[train_vertices] - mu[train_vertices])^2 ))
rmse_test["gs2"] <-
  sqrt(mean( (fit_gs2$yhat.test.mean[test_vertices] - mu[test_vertices])^2 ))
timing["gs2"] <- sum(fit_gs2$timing)



fit_gs3 <-
  flexBART(formula = Y ~ bart(X1),
           train_data = train_data,
           test_data = test_data,
           inform_sigma = TRUE,
           adjacency_list = adjacency_list,
           graph_cut_type = 3,
           n.chains = 4)
rmse_train["gs3"] <-
  sqrt(mean( (fit_gs3$yhat.test.mean[train_vertices] - mu[train_vertices])^2 ))
rmse_test["gs3"] <-
  sqrt(mean( (fit_gs3$yhat.test.mean[test_vertices] - mu[test_vertices])^2 ))
timing["gs3"] <- sum(fit_gs3$timing)

fit_gs4 <-
  flexBART(formula = Y ~ bart(X1),
           train_data = train_data,
           test_data = test_data,
           inform_sigma = TRUE,
           adjacency_list = adjacency_list,
           graph_cut_type = 4,
           n.chains = 4)
rmse_train["gs4"] <-
  sqrt(mean( (fit_gs4$yhat.test.mean[train_vertices] - mu[train_vertices])^2 ))
rmse_test["gs4"] <-
  sqrt(mean( (fit_gs4$yhat.test.mean[test_vertices] - mu[test_vertices])^2 ))
timing["gs4"] <- sum(fit_gs4$timing)


fit_noadj <-
  flexBART(formula = Y ~ bart(X1),
           train_data = train_data,
           test_data = test_data,
           inform_sigma = TRUE,
           n.chains = 4)
rmse_train["noadj"] <-
  sqrt(mean( (fit_noadj$yhat.test.mean[train_vertices] - mu[train_vertices])^2 ))
rmse_test["noadj"] <-
  sqrt(mean( (fit_noadj$yhat.test.mean[test_vertices] - mu[test_vertices])^2 ))
timing["noadj"] <- sum(fit_noadj$timing)


##########
tmp_bart <- matrix(nrow = n, ncol = 4)
tmp_time <- 0
for(cix in 1:4){
  bart_time <-
    system.time(
      bart_fit <- BART::wbart(x.train = train_data[,colnames(train_data) != "Y"],
                              y.train = train_data[,"Y"],
                              x.test = test_data[,colnames(test_data) != "Y"],
                              sparse = TRUE,
                              ndpost = 1000, nskip = 1000))
  tmp_bart[,cix] <- bart_fit$yhat.test.mean
  tmp_time <- tmp_time + bart_time["elapsed"]
}
tmp_bart <- rowMeans(tmp_bart)
rmse_train["bart"] <- 
  sqrt(mean( (tmp_bart[train_vertices] - mu[train_vertices])^2 ))
rmse_test["bart"] <- 
  sqrt(mean( (tmp_bart[test_vertices] - mu[test_vertices])^2 ))
timing["bart"] <- tmp_time

tmp_dbart <- matrix(nrow = n, ncol = 4)
tmp_time <- 0
for(cix in 1:4){
  dbart_time <-
    system.time(
      dbart_fit <- dbarts::bart(x.train = train_model_mat,
                                y.train = train_data[,"Y"],
                                x.test = test_model_mat,
                                ndpost = 1000, nskip = 1000, keeptrees = TRUE))
  
  tmp_dbart[,cix] <- dbart_fit$yhat.test.mean
  tmp_time <- tmp_time + dbart_time["elapsed"]
}
tmp_dbart <- rowMeans(tmp_dbart)
rmse_train["dbart"] <- 
  sqrt(mean( (tmp_dbart[train_vertices] - mu[train_vertices])^2 ))
rmse_test["dbart"] <- 
  sqrt(mean( (tmp_dbart[test_vertices] - mu[test_vertices])^2 ))
timing["dbart"] <- tmp_time

#########################
mu_lim <- 
  c(-1,1) * 
  max(abs(c(mu, fit_noadj$yhat.test.mean, 
            fit_gs1$yhat.test.mean, fit_gs2$yhat.test.mean,
            fit_gs3$yhat.test.mean, fit_gs4$yhat.test.mean,
            tmp_bart, tmp_dbart)))

g_noadj <- g
g_gs1 <- g
g_gs2 <- g
g_gs3 <- g
g_gs4 <- g
g_bart <- g
g_dbart <- g

scaled_mu <- scales::rescale(mu, to = c(0,1), from = mu_lim)
scaled_noadj <-
  scales::rescale(fit_noadj$yhat.test.mean, to = c(0,1), from = mu_lim)
scaled_gs1 <- 
  scales::rescale(fit_gs1$yhat.test.mean, to = c(0,1), from = mu_lim)
scaled_gs2 <- 
  scales::rescale(fit_gs2$yhat.test.mean, to = c(0,1), from = mu_lim)
scaled_gs3 <- 
  scales::rescale(fit_gs3$yhat.test.mean, to = c(0,1), from = mu_lim)
scaled_gs4 <- 
  scales::rescale(fit_gs4$yhat.test.mean, to = c(0,1), from = mu_lim)

scaled_bart <- scales::rescale(tmp_bart, to = c(0,1), from = mu_lim)
scaled_dbart <- scales::rescale(tmp_dbart, to = c(0,1), from = mu_lim)

V(g)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_mu)/255)
V(g_noadj)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_noadj)/255)
V(g_gs1)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_gs1)/255)
V(g_gs2)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_gs2)/255)
V(g_gs3)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_gs3)/255)
V(g_gs4)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_gs4)/255)
V(g_bart)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_bart)/255)
V(g_dbart)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_dbart)/255)

par(mar = c(1,1,1,1), mgp = c(1.8, 0.5, 0), mfrow = c(3,3))
plot(g, layout = layout_on_grid, vertex.label = NA, main = "Truth")
plot(g_noadj, layout = layout_on_grid, vertex.label = NA,main = "noadj")
plot(g_gs1, layout = layout_on_grid, vertex.label = NA,main = "gs1")
plot(g_gs2, layout = layout_on_grid, vertex.label = NA,main = "gs2")
plot(g_gs3, layout = layout_on_grid, vertex.label = NA,main = "gs3")
plot(g_gs4, layout = layout_on_grid, vertex.label = NA,main = "gs4")
plot(g_bart, layout = layout_on_grid, vertex.label = NA,main = "bart")
plot(g_dbart, layout = layout_on_grid, vertex.label = NA,main = "dbart")






