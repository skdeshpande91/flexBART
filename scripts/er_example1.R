library(igraph)
library(RColorBrewer)
library(scales)
library(Rcpp)
library(RcppArmadillo)
library(BART)
library(dbarts)
# currently it exports to R as a hidden function (with a . prefix)
sourceCpp("flexBART/src/flexBART_fit.cpp") 

###################

if(!file.exists("data/er_example1.RData")){
  n <- 25 # number of verticies
  g <- erdos.renyi.game(n = n, p = 0.2, type = "gnp", directed = FALSE)
  fixed_layout <- layout_nicely(g)
  plot(g, layout = fixed_layout)
  A <- as_adjacency_matrix(g, type = "both", sparse = FALSE)
  D <- diag(rowSums(A))
  L <- D - A # unnormalized laplacian
  L_eig <- eigen(L)
  
  mu <-  L_eig$vectors[,(n-1):(n-5)] %*% runif(n = 5, min = 5, max = 10)
  save(n, g, A, mu, fixed_layout, file = "data/er_example2.RData")
} else{
  load("data/er_example1.RData")
}

col_list <- rev(brewer.pal(n = 8, name = "RdYlBu")) # color scale
g_true <- g
scaled_mu <- rescale(mu, to = c(0,1), from = max(abs(mu)) * c(-1.1, 1.1))
V(g_true)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_mu)/255)
plot(g_true, layout = fixed_layout, main = "True underlying function", vertex.label = "")


# We will observe T = 5 observations at each vertex

#sigma <- 2
#T <- 5

sigma <- 1
#T <- 1
T <- 5

N <- T*n
#N <- 10 * n

X_cat <- matrix(rep(0:(n-1), each = T), nrow = N, ncol = 1)
Y <- rep(mu, each = T) + sigma * rnorm(N, mean = 0, sd = 1)

df_X <- data.frame(X = rep(0:(n-1), each = T))
df_X[,"X"] <- as.factor(df_X[,"X"])
X_dummy <- makeModelMatrixFromDataFrame(df_X)

# We have to tell flexBART about the adjacency structure
# we do this by numbering the edges in the adjacency graph
# and passing a list of these edge numbers
A_lower <- A
A_lower[upper.tri(A_lower)] <- 0
adj_support_list <- list(which(A_lower !=0) -1) # C++ is 0-indexed
cat_levels_list <- list(0:(n-1))

M <- 200
y_mean <- mean(Y)
y_sd <- sd(Y)
std_Y_train <- (Y - y_mean)/y_sd # standardize the output
tau <- (max(std_Y_train) - min(std_Y_train))/(2 * 2 * sqrt(M)) # CGM10 prior sd on all leaf parameters

# sigma^2 is given an Inv. Gamma(nu/2, nu * lambda/2) prior
# we're using the standard BART choices for nu & lambda
nu <- 3
lambda <- qchisq(0.1, df = nu)/nu

# Fit 0: predict Y using only the adjacency structure
# graph_split = TRUE tells flexBART to partition levels by cutting an MST
# graph_cut_type = 0: delete an edge from the MST by uniformly deleting an edge
# unif_cuts = TRUE: ignored in this fit
adj_fit <- .flexBART_fit(Y_train = std_Y_train,
                         tX_cont_train = matrix(0, nrow = 1, ncol = 1), # default value to signal no continuous preds
                         tX_cat_train = t(X_cat),
                         tX_cont_test =  matrix(0, nrow = 1, ncol = 1),
                         tX_cat_test = matrix(0, nrow = 1, ncol = 1),
                         unif_cuts = rep(FALSE, times = 1), 
                         cutpoints_list = NULL,
                         cat_levels_list = cat_levels_list,
                         graph_split = rep(TRUE, times = 1), graph_cut_type = 0,
                         adj_support_list = adj_support_list,
                         rc_split = FALSE, prob_rc = 0.0, a_rc = 1, b_rc = 1,
                         sparse = FALSE, a_u = 0.5, b_u = 1, 
                         mu0 = 0, tau = tau,
                         lambda = lambda, nu = nu,
                         M = M, 
                         nd = 1000, burn = 1000, thin = 1,
                         save_trees = FALSE, verbose = TRUE, print_every = 50)

noadj_fit <- .flexBART_fit(Y_train = std_Y_train,
                           tX_cont_train = matrix(0, nrow = 1, ncol = 1), # default value to signal no continuous preds
                           tX_cat_train = t(X_cat),
                           tX_cont_test =  matrix(0, nrow = 1, ncol = 1),
                           tX_cat_test = matrix(0, nrow = 1, ncol = 1),
                           unif_cuts = rep(FALSE, times = 1), 
                           cutpoints_list = NULL,
                           cat_levels_list = cat_levels_list,
                           graph_split = rep(FALSE, times = 1), graph_cut_type = 0, # don't use adjacency
                           adj_support_list = adj_support_list,
                           rc_split = FALSE, prob_rc = 0.0, a_rc = 1, b_rc = 1,
                           sparse = FALSE, a_u = 0.5, b_u = 1, 
                           mu0 = 0, tau = tau,
                           lambda = lambda, nu = nu,
                           M = M, 
                           nd = 1000, burn = 1000, thin = 1,
                           save_trees = FALSE, verbose = TRUE, print_every = 50)
cutpoints_list <- list()
for(i in 1:n) cutpoints_list[[i]] <- c(0,1)
my_bart_fit <- .flexBART_fit(Y_train = std_Y_train,
                             tX_cont_train = t(X_dummy), 
                             tX_cat_train = matrix(0, nrow = 1, ncol = 1), # default value to signal no continuous preds
                             tX_cont_test =  matrix(0, nrow = 1, ncol = 1),
                             tX_cat_test = matrix(0, nrow = 1, ncol = 1),
                             unif_cuts = rep(FALSE, times = n), 
                             cutpoints_list = cutpoints_list,
                             cat_levels_list = NULL,
                             graph_split = rep(FALSE, times = 1), graph_cut_type = 1, # graph_cut ignored
                             adj_support_list = NULL,
                             rc_split = FALSE, prob_rc = 0.0, a_rc = 1, b_rc = 1,
                             sparse = FALSE, a_u = 0.5, b_u = 1, 
                             mu0 = 0, tau = tau,
                             lambda = lambda, nu = nu,
                             M = M, 
                             nd = 1000, burn = 1000, thin = 1,
                             save_trees = FALSE, verbose = TRUE, print_every = 50)


bart_fit <- wbart(x.train = X_dummy, y.train = Y,
                  ntree = M, ndpost = 1000, nskip = 1000, printevery = 100)

# There's a lot of duplicate entries

index <- 1 + (1:n - 1) * T
adj_mean <- y_mean + y_sd * colMeans(adj_fit$fit_train[,index])

noadj_mean <- y_mean + y_sd * colMeans(noadj_fit$fit_train[,index])

my_bart_mean <- y_mean + y_sd * colMeans(my_bart_fit$fit_train[,index])

bart_mean <- bart_fit$yhat.train.mean[index]



mean( (mu - adj_mean)^2 )
mean( (mu - noadj_mean)^2 )
mean( (mu - bart_mean)^2 )
mean( (mu - my_bart_mean)^2 )


#################
muhat_lim <- max(abs(c(mu,adj_mean,noadj_mean, my_bart_mean,bart_mean)))

g_adj <- g
scaled_mu <- rescale(adj_mean, to = c(0,1), from = muhat_lim * c(-1.01, 1.01))
V(g_adj)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_mu)/255)
#plot(g_adj, layout = fixed_layout, main = "flexBART (w/ adjacency)", vertex.label = "", vertex.size = 15)

g_noadj <- g
scaled_mu <- rescale(noadj_mean, to = c(0,1), from = muhat_lim * c(-1.01, 1.01))
V(g_noadj)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_mu)/255)
#plot(g_noadj, layout = fixed_layout, main = "flexBART (w/o adjacency)", vertex.label = "", vertex.size = 15)

g_bart <- g
scaled_mu <- rescale(bart_mean, to = c(0,1), from = muhat_lim * c(-1.01, 1.01))
V(g_bart)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_mu)/255)
#plot(g_bart, layout = fixed_layout, main = "BART", vertex.label = "", vertex.size = 15)


g_mybart <- g
scaled_mu <- rescale(my_bart_mean, to = c(0,1), from = muhat_lim * c(-1.01, 1.01))
V(g_mybart)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_mu)/255)
#plot(g_mybart, layout = fixed_layout, main = "myBART", vertex.label = "", vertex.size = 15)

scaled_mu <- rescale(mu, to = c(0,1), from = muhat_lim * c(-1.01, 1.01))
V(g_true)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_mu)/255)


par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(2,3))

plot(g_true, layout = fixed_layout, main = "True function", vertex.label = "", vertex.size = 15)

plot(g_bart, layout = fixed_layout, main = "BART", vertex.label = "", vertex.size = 15)

plot(g_mybart, layout = fixed_layout, main = "BART (no node min)", vertex.label = "", vertex.size = 15)

plot(g_noadj, layout = fixed_layout, main = "flexBART (w/o adjacency)", vertex.label = "", vertex.size = 15)

plot(g_adj, layout = fixed_layout, main = "flexBART (w/ adjacency)", vertex.label = "", vertex.size = 15)
