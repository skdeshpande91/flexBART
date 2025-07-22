timing_comparison
================
2025-05-07

**flexBART** is often much faster than **BART** because it avoids a lot
of redundant computations. The speedup really manifests when you have
lots of data. To illustrate, let’s look compare the two implementations
on an easy task (or at least, easy for BART): a slightly modified
version of the [Friedman
function](https://www.sfu.ca/~ssurjano/fried.html)!

``` r
mu_true <- function(df){
  if(!all(abs(df) <= 1)){
    stop("all entries in x_cont must be between -1 and 1")
  } else{
    tmp_X <- (df+1)/2
    return(10*sin(pi*tmp_X[,1] * tmp_X[,2]) + 
           20 * (tmp_X[,8] - 0.5)^2 + 
           10 * tmp_X[,17] + 
           5 * tmp_X[,20])
  }
}
```

Note: the Friedman function is defined over the unit square but
**flexBART** assumes that continuous predictors lie in \[-1,1\]. So in
this example, we had to transform our predictors to the unit interval.

Now let’s generate some data.

``` r
n_train <- 10000
n_test <- 1000
p_cont <- 50
sigma <- 1
train_data <- data.frame(Y = rep(NA, times = n_train))
for(j in 1:p_cont) train_data[[paste0("X",j)]] <- runif(n_train, min = -1, max = 1)
mu_train <- mu_true(train_data[,paste0("X",1:p_cont)])
train_data[,"Y"] <- mu_train + sigma * rnorm(n = n_train, mean = 0, sd = 1)

test_data <- data.frame(Y = rep(NA, times = n_test))
for(j in 1:p_cont) test_data[[paste0("X",j)]] <- runif(n_test, min = -1, max = 1)
mu_test <- mu_true(test_data[,paste0("X",1:p_cont)])

# Containers to store the performance results
rmse_train <- c("flexBART" = NA, "BART" = NA)
rmse_test <- c("flexBART" = NA, "BART" = NA)
timing <- c("flexBART" = NA, "BART" = NA)
```

By default, `flexBART::flexBART` simulates 4 Markov chains for 2000
iterations each.

``` r
flex_fit <-
  flexBART::flexBART(formula = Y~bart(.),
                     train_data = train_data,
                     test_data = test_data,
                     M = 200)
rmse_train["flexBART"] <- sqrt(mean( (mu_train - flex_fit$yhat.train.mean)^2 ))
rmse_test["flexBART"] <- sqrt(mean( (mu_test - flex_fit$yhat.test.mean)^2 ))
timing["flexBART"] <- sum(flex_fit$timing) # total run time over all chains
```

To make the timing comparison fair, we’ll run 4 `BART::wbart` chains for
the same number of iterations.

``` r
bart_time <- rep(NA, times = 4)
bart_train <- rep(0, times = n_train)
bart_test <- rep(0, times = n_test)
for(cix in 1:4){
  tmp_time <-
    system.time(
      bart_fit <- 
        BART::wbart(x.train = train_data[,colnames(train_data) != "Y"], 
                    y.train = train_data[,"Y"], 
                    x.test = test_data[,colnames(test_data) != "Y"],
                  ndpost = 1000, nskip = 1000))
  bart_train <-
    bart_train + bart_fit$yhat.train.mean/4
  bart_test <-
    bart_test + bart_fit$yhat.test.mean/4
  bart_time[cix] <- tmp_time["elapsed"]
}
rmse_train["BART"] <- sqrt(mean( (mu_train - bart_train)^2 ))
rmse_test["BART"] <- sqrt(mean( (mu_test - bart_test)^2 ))
timing["BART"] <- sum(bart_time)
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
    ##    0.239    0.255 
    ## [1] "Testing RMSE"
    ## flexBART     BART 
    ##    0.241    0.258 
    ## [1] "Timing (seconds):"
    ## flexBART     BART 
    ##   87.524  383.035
