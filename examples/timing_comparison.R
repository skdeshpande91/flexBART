# Comparing run times between flexBART and BART::wbart
# For this demonstration, we'll use the Friedman function, which we know BART fits well

n <- 10000
n_test <- 1000
p_cont <- 10
set.seed(99)
X <- matrix(runif(n*p_cont, min = -1, max = 1), nrow = n, ncol = p_cont)
X_test <- matrix(runif(n_test * p_cont, min = -1, max = 1), nrow = n, ncol = p_cont)

f <- function(x_cont){
  if(!all(abs(x_cont) <= 1)){
    stop("all entries in x_cont must be between -1 and 1")
  } else{
    x <- (x_cont+1)/2 # convert to [0,1]
    return(10 * sin(pi*x[,1]*x[,2]) + 20 * (x[,3] - 0.5)^2 + 10*x[,4] + 5 * x[,5])
  }
}

sigma <- 2.5
set.seed(99)
mu <- f(X)
mu_test <- f(X_test)
y <- mu + sigma * rnorm(n, mean = 0, sd = 1)

flex_fit <- flexBART::flexBART(Y.train = y, X_cont_train = X)