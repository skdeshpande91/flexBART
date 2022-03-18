library(igraph)
library(RColorBrewer)
library(scales)
library(Rcpp)
library(RcppArmadillo)

# currently it exports to R as a hidden function (with a . prefix)
sourceCpp("flexBART/src/flexBART_fit.cpp") 

######################
# Generate data on a 10 x 10 grid
n_side <- 10
n <- 100

g <- make_lattice(length = n_side, dim = 2)
A <- as_adjacency_matrix(g, type = "both", sparse = FALSE)

# make up 5 clusters
cluster1 <- 81:100
cluster1 <- cluster1[!cluster1 %in% c(89,90, 100)]
cluster2 <- c(89, 90, 100, 79, 80, 69, 70, 59, 60, 50, 40)
cluster3 <- rep(1:4, times = 4) + 10 * rep(0:3, each = 4) # lower 4x4 grid
cluster4 <- c(52, 53,63,54, 64, 55,65,56, 66,46, 
              47, 57, 67, 37, 38, 48, 28, 29, 27, 
              17, 16, 18, 6, 7, 8,30, 20, 10, 9, 19, 39, 49)
cluster5 <- which(!1:100 %in% c(cluster1, cluster2, cluster3, cluster4))

cluster_means <- c(7.5, 3, -2, -5, 0)
mu <- rep(NA, times = n)
mu[cluster1] <- cluster_means[1]
mu[cluster2] <- cluster_means[2]
mu[cluster3] <- cluster_means[3]
mu[cluster4] <- cluster_means[4]
mu[cluster5] <- cluster_means[5]


sigma <- 1
set.seed(129)
Y_all <- mu + sigma * rnorm(n, 0, 1)

# quickly visualize the true underlying function
col_list <- rev(brewer.pal(n = 8, name = "RdYlBu")) # color scale
g_true <- g
scaled_mu <- rescale(mu, to = c(0,1), from = c(-8,8))
V(g_true)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_mu)/255)
plot(g_true, layout = layout_on_grid, main = "True underlying function", vertex.label = "")

# Represent each spatial unit as level of a single categorical variable
X_cat_all <- matrix(as.integer(1:n) - 1, nrow = n, ncol = 1) # note C++ is 0-indexed

# For comparison with/ regular BART we need to also
# represent each spatial unit with a dummy variable
X_cont_all <- matrix(0, nrow = n, ncol = n)
diag(X_cont_all) <- 1


# To facilitate comparison with regular BART, we need to pre-specify cut-points
cutpoints_list <- list()
for(i in 1:n) cutpoints_list[[i]] <- as.numeric(0:1)

# We have to tell flexBART about the adjacency structure
# we do this by numbering the edges in the adjacency graph
# and passing a list of these edge numbers
A_lower <- A
A_lower[upper.tri(A_lower)] <- 0
adj_support_list <- list(which(A_lower !=0) -1) # C++ is 0-indexed

# flexBART needs to label the levels of the categorical variable
# this labelling can be arbitrary so long as it's consistent 
# with the adjacency structure encoded by A
# in this case it's just 0, ..., n-1 (remembering C++ is 0-index)
cat_levels_list <- list(0:(n-1))

###############
# Create a single test/train split:
#set.seed(211)
set.seed(31822)
test_index <- sort(sample(1:n, size = 10, replace = FALSE))
train_index <- (1:n)[-test_index]

X_cont_train <- X_cont_all[train_index,]
X_cat_train <- X_cat_all[train_index,]
X_cont_test <- X_cont_all[test_index,]
X_cat_test <- X_cat_all[test_index,]
Y_train <- Y_all[train_index]
Y_test <- Y_all[test_index]

###
# Some pre-processing
# to set certain prior hyperparameters
M <- 200
y_mean <- mean(Y_train)
y_sd <- sd(Y_train)
std_Y_train <- (Y_train - y_mean)/y_sd # standardize the output
tau <- (max(std_Y_train) - min(std_Y_train))/(2 * 2 * sqrt(M)) # CGM10 prior sd on all leaf parameters

# sigma^2 is given an Inv. Gamma(nu/2, nu * lambda/2) prior
# we're using the standard BART choices for nu & lambda
nu <- 3
lambda <- qchisq(0.1, df = nu)/nu


# Fit 0: predict Y using only the adjacency structure
# graph_split = TRUE tells flexBART to partition levels by cutting an MST
# graph_cut_type = 0: delete an edge from the MST by uniformly deleting an edge
# unif_cuts = TRUE: ignored in this fit
fit0 <- .flexBART_fit(Y_train = std_Y_train,
                      tX_cont_train = matrix(0, nrow = 1, ncol = 1), # default value to signal no continuous preds
                      tX_cat_train = t(X_cat_train),
                      tX_cont_test =  matrix(0, nrow = 1, ncol = 1),
                      tX_cat_test = t(X_cat_test),
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
# Fit 1: predict Y using only the adjacency structure
# graph_split = TRUE tells flexBART to partition levels by cutting an randomly weighted MST
# graph_cut_type = 1: delete an edge from the MST w.p. proportional to the (random) uniform weight
# unif_cuts = TRUE: ignored in this fit
fit1 <- .flexBART_fit(Y_train = std_Y_train,
                      tX_cont_train = matrix(0, nrow = 1, ncol = 1), # default value to signal no continuous preds
                      tX_cat_train = t(X_cat_train),
                      tX_cont_test =  matrix(0, nrow = 1, ncol = 1),
                      tX_cat_test = t(X_cat_test),
                      unif_cuts = rep(FALSE, times = 1), 
                      cutpoints_list = NULL,
                      cat_levels_list = cat_levels_list,
                      graph_split = rep(TRUE, times = 1), graph_cut_type = 1,
                      adj_support_list = adj_support_list,
                      rc_split = FALSE, prob_rc = 0.0, a_rc = 1, b_rc = 1,
                      sparse = FALSE, a_u = 0.5, b_u = 1, 
                      mu0 = 0, tau = tau,
                      lambda = lambda, nu = nu,
                      M = M, 
                      nd = 1000, burn = 1000, thin = 1,
                      save_trees = FALSE, verbose = TRUE, print_every = 50)

# Fit 2: predict Y using only the adjacency structure
# graph_split = TRUE tells flexBART to partition levels by cutting an randomly weighted MST
# graph_cut_type = 2: delete an edge from the MST w.p. proportional to size of smallest cluster formed
# unif_cuts = TRUE: ignored in this fit
fit2 <- .flexBART_fit(Y_train = std_Y_train,
                      tX_cont_train = matrix(0, nrow = 1, ncol = 1), # default value to signal no continuous preds
                      tX_cat_train = t(X_cat_train),
                      tX_cont_test =  matrix(0, nrow = 1, ncol = 1),
                      tX_cat_test = t(X_cat_test),
                      unif_cuts = rep(FALSE, times = 1), 
                      cutpoints_list = NULL,
                      cat_levels_list = cat_levels_list,
                      graph_split = rep(TRUE, times = 1), graph_cut_type =2,
                      adj_support_list = adj_support_list,
                      rc_split = FALSE, prob_rc = 0.0, a_rc = 1, b_rc = 1,
                      sparse = FALSE, a_u = 0.5, b_u = 1, 
                      mu0 = 0, tau = tau,
                      lambda = lambda, nu = nu,
                      M = M, 
                      nd = 1000, burn = 1000, thin = 1,
                      save_trees = FALSE, verbose = TRUE, print_every = 50)

# Fit 3: predict Y using only the adjacency structure
# graph_split = TRUE tells flexBART to partition levels by cutting an randomly weighted MST
# graph_cut_type = 2: delete the edge from the MST with largest (random) weight
# unif_cuts = TRUE: ignored in this fit
fit3 <- .flexBART_fit(Y_train = std_Y_train,
                      tX_cont_train = matrix(0, nrow = 1, ncol = 1), # default value to signal no continuous preds
                      tX_cat_train = t(X_cat_train),
                      tX_cont_test =  matrix(0, nrow = 1, ncol = 1),
                      tX_cat_test = t(X_cat_test),
                      unif_cuts = rep(FALSE, times = 1), 
                      cutpoints_list = NULL,
                      cat_levels_list = cat_levels_list,
                      graph_split = rep(TRUE, times = 1), graph_cut_type = 3,
                      adj_support_list = adj_support_list,
                      rc_split = FALSE, prob_rc = 0.0, a_rc = 1, b_rc = 1,
                      sparse = FALSE, a_u = 0.5, b_u = 1, 
                      mu0 = 0, tau = tau,
                      lambda = lambda, nu = nu,
                      M = M, 
                      nd = 1000, burn = 1000, thin = 1,
                      save_trees = FALSE, verbose = TRUE, print_every = 50)

# Fit 4: predict Y using only the index of the observation (i.e. ignore adjacency)
# instead levels are partitionined uniformly at random
fit4 <- .flexBART_fit(Y_train = std_Y_train,
                      tX_cont_train = matrix(0, nrow = 1, ncol = 1), # default value to signal no continuous preds
                      tX_cat_train = t(X_cat_train),
                      tX_cont_test =  matrix(0, nrow = 1, ncol = 1),
                      tX_cat_test = t(X_cat_test),
                      unif_cuts = rep(FALSE, times = 1), 
                      cutpoints_list = NULL,
                      cat_levels_list = cat_levels_list,
                      graph_split = rep(FALSE, times = 1), graph_cut_type = 1, # graph_cut ignored
                      adj_support_list = adj_support_list,
                      rc_split = FALSE, prob_rc = 0.0, a_rc = 1, b_rc = 1,
                      sparse = FALSE, a_u = 0.5, b_u = 1, 
                      mu0 = 0, tau = tau,
                      lambda = lambda, nu = nu,
                      M = M, 
                      nd = 1000, burn = 1000, thin = 1,
                      save_trees = FALSE, verbose = TRUE, print_every = 50)

# Fit5: vanilla BART
fit5 <- .flexBART_fit(Y_train = std_Y_train,
                      tX_cont_train = t(X_cont_train), 
                      tX_cat_train = matrix(0, nrow = 1, ncol = 1), # default value to signal no continuous preds
                      tX_cont_test =  t(X_cont_test),
                      tX_cat_test = matrix(0, nrow = 1, ncol = 1),
                      unif_cuts = rep(FALSE, times = 1), 
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


##########
# flexBART works with the standardized responses, so we need to convert back 
# to the original scale
# here is the posterior mean on the training data
post_mean_train0 <- y_mean + y_sd * colMeans(fit0$fit_train)
post_mean_train1 <- y_mean + y_sd * colMeans(fit1$fit_train)
post_mean_train2 <- y_mean + y_sd * colMeans(fit2$fit_train)
post_mean_train3 <- y_mean + y_sd * colMeans(fit3$fit_train)
post_mean_train4 <- y_mean + y_sd * colMeans(fit4$fit_train)
post_mean_train5 <- y_mean + y_sd * colMeans(fit5$fit_train)

# and here is the posterior mean on the testing data
post_mean_test0 <- y_mean + y_sd * colMeans(fit0$fit_test)
post_mean_test1 <- y_mean + y_sd * colMeans(fit1$fit_test)
post_mean_test2 <- y_mean + y_sd * colMeans(fit2$fit_test)
post_mean_test3 <- y_mean + y_sd * colMeans(fit3$fit_test)
post_mean_test4 <- y_mean + y_sd * colMeans(fit4$fit_test)
post_mean_test5 <- y_mean + y_sd * colMeans(fit5$fit_test)

#######
# Expectation: fit1 and fit2 should perform quite similarly -- they both
# take the adjacency structure into account. The only difference
# is how they choose which edge in the MST to split; fit1 is biased towards
# creating equally sized clusters while fit2 is slightly biased towards creating singletons
# I suspect fit3 and fit4 will be pretty bad: they will tend to cluster units that are not
# adjacent and we should see some strong shrinkage towards the global average;


# In-sample RMSE:
sqrt(mean( (mu[train_index] - post_mean_train0)^2 ))
sqrt(mean( (mu[train_index] - post_mean_train1)^2 ))
sqrt(mean( (mu[train_index] - post_mean_train2)^2 ))
sqrt(mean( (mu[train_index] - post_mean_train3)^2 ))
sqrt(mean( (mu[train_index] - post_mean_train4)^2 ))
sqrt(mean( (mu[train_index] - post_mean_train5)^2 ))


# Out-of-sample RMSE
sqrt(mean( (mu[test_index] - post_mean_test0)^2 ))
sqrt(mean( (mu[test_index] - post_mean_test1)^2 ))
sqrt(mean( (mu[test_index] - post_mean_test2)^2 ))
sqrt(mean( (mu[test_index] - post_mean_test3)^2 ))
sqrt(mean( (mu[test_index] - post_mean_test4)^2 ))
sqrt(mean( (mu[test_index] - post_mean_test5)^2 ))

# Standardized mean square errors (1 - SMSE is akin to "predictive R^2")
# 1 - SMSE is essentially "how much more variation am I explaining than the mean
# of the training data"
# small values of SMSE are good
mean( (mu[test_index] - post_mean_test0)^2 ) / mean( (mu[test_index] - y_mean)^2 )
mean( (mu[test_index] - post_mean_test1)^2 ) / mean( (mu[test_index] - y_mean)^2 )
mean( (mu[test_index] - post_mean_test2)^2 ) / mean( (mu[test_index] - y_mean)^2 )
mean( (mu[test_index] - post_mean_test3)^2 ) / mean( (mu[test_index] - y_mean)^2 )
mean( (mu[test_index] - post_mean_test4)^2 ) / mean( (mu[test_index] - y_mean)^2 )
mean( (mu[test_index] - post_mean_test5)^2 ) / mean( (mu[test_index] - y_mean)^2 )


###########
# Visualizations

col_list <- rev(brewer.pal(n = 8, name = "RdYlBu"))
max_value <- max(abs(c(Y_all, mu,
                       post_mean_train0, post_mean_test0,
                       post_mean_train1, post_mean_test1, 
                       post_mean_train2, post_mean_test2,
                       post_mean_train3, post_mean_test3,
                       post_mean_train4, post_mean_test4,
                       post_mean_train5, post_mean_test5)))



# Plot the true function mu but this time with the enhanced color scale
scaled_mu <- rescale(mu, to = c(0,1), from = max_value * c(-1,1))
V(g_true)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_mu)/255)
plot(g_true, layout = layout_on_grid, main = "True mu", vertex.label = "")

# Let's plot the original data
g_all <- g
scaled_Y <- scales::rescale(Y_all, to = c(0,1), from = max_value * c(-1,1))
V(g_all)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_Y)/255)
plot(g_all, layout = layout_on_grid, main = "Data", vertex.label = "")

# Let's see testing dataset
g_train <- g
scaled_Y <- scales::rescale(Y_all, to = c(0,1), from = max_value * c(-1,1))
V(g_train)$color <- rgb(colorRamp(col_list, bias = 1)(scaled_Y)/255)
V(g_train)$color[test_index] <- "darkgray" 
plot(g_train, layout = layout_on_grid, main = "Training", vertex.label = "")

for(fit_ix in 0:5){
  fit <- get(paste0("fit",fit_ix))
  tmp_g <- g
  scaled_train <- scales::rescale(get(paste0("post_mean_train",fit_ix)), to = c(0,1), from = max_value * c(-1,1))
  scaled_test <- scales::rescale(get(paste0("post_mean_test",fit_ix)), to = c(0,1), from = max_value * c(-1,1))
  V(tmp_g)$color <- rep(NA, times = n)
  V(tmp_g)$color[train_index] <- rgb(colorRamp(col_list, bias = 1)(scaled_train)/255)
  V(tmp_g)$color[test_index] <- rgb(colorRamp(col_list, bias = 1)(scaled_test)/255)
  assign(paste0("g_fit", fit_ix), tmp_g)
}

png("figures/grid_example.png", width = 8, height = 6, units = "in", res = 600)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(2,3))
plot(g_true, layout = layout_on_grid, main = "True mu", vertex.label = "")
plot(g_all, layout = layout_on_grid, main = "Data", vertex.label = "")
plot(g_train, layout = layout_on_grid, main = "Training", vertex.label = "")
plot(g_fit2, layout = layout_on_grid, main = "flexBART (w/ adjacency)", vertex.label = "")
plot(g_fit4, layout = layout_on_grid, main = "flexBART (w/o adjacency)", vertex.label = "")
plot(g_fit5, layout = layout_on_grid, main = "Regular BART", vertex.label = "")
dev.off()

