library(Rcpp)
library(RcppArmadillo)
library(BART)

sourceCpp("flexBART/src/flexBARTfast_fit.cpp")

p_cont <- 2
p_cat <- 0
n <- 10000
set.seed(129)
X_cont <- matrix(runif(n * p_cont, -1, 1), nrow = n, ncol = p_cont)
X_cont[,2] <- 1*(X_cont[,2] >= 0)
true_mu <- function(X_cont){
  scaled_X_cont <- (X_cont + 1)/2 # moves it to [0,1]
  z1 <- scaled_X_cont[,1]
  #z2 <- 1*(X_cont[,2] >= 0)
  z2 <- X_cont[,2]
  return(3 * z1 + (2 - 5 * z2) * sin(pi * z1) - 2 * z2)
}

mu <- true_mu(X_cont)

sigma <- 2
Y_train <- mu + sigma* rnorm(n = n, mean = 0, sd = 1)


M <- 200
y_mean <- mean(Y_train)
y_sd <- sd(Y_train)
std_Y_train <- (Y_train - y_mean)/y_sd # standardize the output
tau <- (max(std_Y_train) - min(std_Y_train))/(2 * 2 * sqrt(M)) # CGM10 prior sd on all leaf parameters
nu <- 3
lambda <- qchisq(0.1, df = nu)/nu



unif_cuts <- c(FALSE, FALSE)
cutpoints_list <- list()
cutpoints_list[[1]] <- seq(min(X_cont[,1]), max(X_cont[,1]), length = 100)
cutpoints_list[[2]] <- c(0,1)

bart_fit <- wbart(x.train = X_cont, y.train = Y_train)


start_time <- Sys.time()
my_bart_fit <- .flexBARTfast_fit(Y_train = std_Y_train,
                                 tX_cont_train = t(X_cont),
                                 tX_cat_train = matrix(1, nrow = 1, ncol = 1),
                                 tX_cont_test = matrix(1, nrow = 1, ncol = 1),
                                 tX_cat_test = matrix(1, nrow = 1, ncol = 1),
                                 unif_cuts = unif_cuts,
                                 cutpoints_list = cutpoints_list,
                                 cat_levels_list = NULL,
                                 graph_split = c(FALSE,FALSE),
                                 graph_cut_type = 0, 
                                 adj_support_list = NULL,
                                 rc_split = FALSE, prob_rc = 0, a_rc = 1, b_rc = 1,
                                 sparse = FALSE, a_u = 0.5, b_u = 1, 
                                 mu0 = 0, tau = tau,
                                 lambda = lambda, nu = nu,
                                 M = M, nd = 1000, burn = 1000, thin = 1,
                                 save_trees = FALSE, verbose = FALSE, print_every = 50)
end_time <- Sys.time()

bart_mean <- bart_fit$yhat.train.mean
my_mean <- y_mean + y_sd * colMeans(my_bart_fit$fit_train)

plot(my_mean, bart_mean, pch = 16, cex = 0.5, col = 'red')


mean( (mu - my_mean)^2)
mean( (mu - bart_mean)^2)
