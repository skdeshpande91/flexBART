timing_comparison
================
2023-05-12

**flexBART** is often much faster than **BART** because it avoids a lot
of redundant computations. The speedup really manifests when you have
lots of data. To illustrate, let’s look compare the two implementations
on an easy task (or at least, easy for BART): the [Friedman
function](https://www.sfu.ca/~ssurjano/fried.html)!

``` r
f <- function(x_cont){
  if(ncol(x_cont) < 5) stop("x_cont needs to have at least five columns")
  
  if(!all(abs(x_cont) <= 1)){
    stop("all entries in x_cont must be between -1 and 1")
  } else{
    x <- (x_cont+1)/2 # convert to [0,1]
    return(10 * sin(pi*x[,1]*x[,2]) + 20 * (x[,3] - 0.5)^2 + 10*x[,4] + 5 * x[,5])
  }
}
```

Note: the Friedman function is defined over the unit square but
**flexBART** assumes that continuous predictors lie in \[-1,1\]. So in
this example, we had to transform our predictors to the unit interval.

Now let’s generate some data.

``` r
n <- 10000
n_test <- 1000
p_cont <- 10
set.seed(99)
X <- matrix(runif(n*p_cont, min = -1, max = 1), nrow = n, ncol = p_cont)
X_test <- matrix(runif(n_test * p_cont, min = -1, max = 1), nrow = n, ncol = p_cont)

mu <- f(X)
mu_test <- f(X_test)

sigma <- 1
set.seed(99)

y <- mu + sigma * rnorm(n, mean = 0, sd = 1)

# Containers to store the performance results
rmse_train <- c("flexBART" = NA, "BART" = NA)
rmse_test <- c("flexBART" = NA, "BART" = NA)
timing <- c("flexBART" = NA, "BART" = NA)
```

We’re now ready to run both `BART::wbart` and `flexBART::flexBART`. Note
that we have hidden the printed output.

``` r
flex_time <-
  system.time(
    flex_fit <- 
      flexBART::flexBART(Y_train = y, 
                         X_cont_train = X, X_cont_test = X_test))
rmse_train["flexBART"] <- sqrt(mean( (mu - flex_fit$yhat.train.mean)^2 ))
rmse_test["flexBART"] <- sqrt(mean( (mu_test - flex_fit$yhat.test.mean)^2 ))
timing["flexBART"] <- flex_time["elapsed"]
```

``` r
bart_time <-
  system.time(
    bart_fit <-
      BART::wbart(x.train = X, y.train = y, x.test = X_test,
                  ndpost = 1000, nskip = 1000))
rmse_train["BART"] <- sqrt(mean( (mu - bart_fit$yhat.train.mean)^2 ))
rmse_test["BART"] <- sqrt(mean( (mu_test - bart_fit$yhat.test.mean)^2 ))
timing["BART"] <- bart_time["elapsed"]
```

``` r
print("Training RMSE")
print(round(rmse_train, digits = 3))

print("Testing RMSE")
print(round(rmse_test, digits = 3))

print("Timing (seconds):")
print(round(timing, digits = 3))
```

    ## [1] "Training RMSE"
    ## flexBART     BART 
    ##    0.305    0.303 
    ## [1] "Testing RMSE"
    ## flexBART     BART 
    ##    0.343    0.338 
    ## [1] "Timing (seconds):"
    ## flexBART     BART 
    ##   35.251  111.539
