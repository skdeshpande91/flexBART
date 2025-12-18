# Friedman function example

We will demonstrate the basic functionality of **flexBART** using a
slightly modified version of the [Friedman
function](https://www.sfu.ca/~ssurjano/fried.html), which is often used
to check and benchmark BART implementations.

We begin by defining the Friedman function.

``` r
friedman_func <- function(df){
  if(!all(abs(df[,1:5]-0.5) <= 1)){
    stop("all entries in the first 5 columns of df must be between 0 and 1")
  } else{
    
    return(10*sin(pi*df[,1] * df[,2]) + 
           20 * (df[,3] - 0.5)^2 + 
           10 * df[,4] + 
           5 * df[,5])
  }
}
```

Although the function depends on only 5 covariates, for this
demonstration, we will create a total of \\p = 50\\ predictors, each
drawn uniformly from the interval \\\[0,1\].\\ We will also add
\\\mathcal{N}(0, 2.5^{2})\\ noise

``` r
set.seed(724)
n_train <- 10000
n_test <- 5000
p_cont <- 50
sigma <- 2.5

train_data <- data.frame(Y = rep(NA, times = n_train))
for(j in 1:p_cont) train_data[[paste0("X",j)]] <- runif(n_train, min = 0, max = 1)
mu_train <- friedman_func(train_data[,paste0("X",1:p_cont)])
train_data[,"Y"] <- mu_train + sigma * rnorm(n = n_train, mean = 0, sd = 1)

test_data <- data.frame(Y = rep(NA, times = n_test))
for(j in 1:p_cont) test_data[[paste0("X",j)]] <- runif(n_test, min = 0, max = 1)
mu_test <- friedman_func(test_data[,paste0("X",1:p_cont)])

# Containers to store the performance results
rmse_train <- c("flexBART" = NA, "BART" = NA)
rmse_test <- c("flexBART" = NA, "BART" = NA)
timing <- c("flexBART" = NA, "BART" = NA)
```

By default,
[`flexBART::flexBART`](https://skdeshpande91.github.io/flexBART/reference/flexBART.md)
simulates 4 Markov chains for 2000 iterations each. It also performs
variable selection using [Linero
(2018)](doi:10.1080/01621459.2016.1264957)’s sparse Dirichlet prior on
splitting probabilities instead of using uniform splitting
probabilities.

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

To make the comparison fair, we’ll run 4
[`BART::wbart`](https://rdrr.io/pkg/BART/man/wbart.html) chains for the
same number of iterations with the argument `sparse = TRUE`.

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
                    sparse = TRUE,
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

Besides the handling of categorical covariates, there are some important
implementation differences between **BART** and **flexBART**. First,
**BART** creates a grid of 100 potential cutpoints (this is controlled
by the `numcut` argument of
[`BART::wbart()`](https://rdrr.io/pkg/BART/man/wbart.html)) for each
continuous predictor **flexBART**, by contrast, draws cutpoints
uniformly from the interval of available values. Second, **BART**
initializes the residual standard deviation \\\sigma\\ based on an
estimated linear model (see
[here](https://github.com/cran/BART/blob/master/R/wbart.R#L95-L101)).
**flexBART** instead uses the root mean square of the residuals from a
cross-validated LASSO fit. Because of these differences, we would not
expect to obtain identical results from **flexxBART** and **BART**.
Nevertheless, at least for these data, the two implementations display
similar in-sample and out-of-sample predictive performance. We see,
additionally, that **flexBART** is much faster than **BART**.

``` r
print("Training RMSE")
#> [1] "Training RMSE"
print(round(rmse_train, digits = 3))
#> flexBART     BART 
#>    0.418    0.414

print("Testing RMSE")
#> [1] "Testing RMSE"
print(round(rmse_test, digits = 3))
#> flexBART     BART 
#>    0.423    0.425

print("Timing (seconds):")
#> [1] "Timing (seconds):"
print(round(timing, digits = 3))
#> flexBART     BART 
#>  142.784  561.390
```
