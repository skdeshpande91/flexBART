\name{flexBART}
\alias{flexBART}

\title{
A more flexible BART
}
\description{
Implements Chipman et al. (2010)'s Bayesian additive regression trees (BART) method for nonparametric regression with continuous outcomes and Deshpande et al. (2024)'s varying coefficient BART (VCBART) model. 
Unknown functions are represented as a sum of binary regression trees. 
flexBART handles categorical outcomes more flexibly than other implementations of BART.
}
\usage{
flexBART(formula, train_data, test_data = NULL, inform_sigma = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class \code{\link[stats:formula]{formula}} (or one that can be coerced to the class): 
  a symbolic description of the model to be fitted. 
  The details of model specification are given under 'Details'}
  \item{train_data}{an object of class \code{data.frame} containing data used to train the model. 
  As usual, rows (resp. columns) correspond to observations (resp. variables)}
  \item{test_data}{an optional object of class \code{data.frame} containing test-set (i.e., out-of-sample) data. 
  Default is \code{NULL}.}
  \item{inform_sigma}{a logical value for whether the residual variance 
  should be initialized using regularized regression (\code{inform_sigm=TRUE}) or not (\code{inform_sigma=FALSE}). Default is \code{TRUE}.}
  \item{\dots}{additional arguments to be passed to the low-level model fitting function (see below)}
}
\details{
Default implementations of Bayesian Additive Regression Trees (BART) represent categorical predictors 
using several binary indicators, one for each level of each categorical predictor. 
Axis-aligned decision rules are well-defined with these indicators; they send one level of a categorical predictor to the left and all other levels to the right (or vice versa). 
Regression trees built with these rules partition the set of all levels of a categorical predictor by recursively removing one level at a time. Unfortunately, most partitions of 
the levels cannot be built with this ``remove one at a time'' strategy, 
meaning that default implementations of BART are extremely limited in their ability to ``borrow strength'' 
across groups of levels.

\code{flexBART} overcomes this limitation using a new prior on regression trees.
Under this new prior, conditional on splitting on a categorical predictor at a particular node in the tree, levels of the predictor are sent to the left and right child uniformly at random. 
In this way, multiple levels of a categorical predictor are able to be clustered together.
}
\value{
An object of \code{class} "flexBART" containing
\item{dinfo}{TO DO}
\item{trees}{A list (or length \code{nd}) of character vectors (of lenght \code{M}) containing textual representations of the regression trees. These strings are parsed by \code{predict_flexBART} to reconstruct the C++ representations of the sampled trees.}
\item{scaling_info}{TO DO: stuff needed for \code{\link{predict.flexBART}}}
\item{M}{Number of trees}
\item{cov_ensm}{TO DO}
\item{is.probit}{a logical value indicating whether we fit a probit model or not}
\item{yhat.train.mean}{Vector containing posterior mean of evaluations of regression function on training data.}
\item{yhat.train}{Matrix with \code{nd} rows and \code{length(Y_train)} columns. Each row corresponds to a posterior sample of the regression function and each column corresponds to a training observation. Only returned if \code{save_samples == TRUE}.}
\item{yhat.test.mean}{Vector containing posterior mean of evaluations of regression function on testing data, if testing data is provided.}
\item{yhat.test}{If testing data was supplied, matrix containing posterior samples of the regression function evaluated on the testing data. Structure is similar to that of \code{yhat_train}. Only returned if testing data is passed and \code{save_samples == TRUE}.}
\item{beta.train.mean}{Matrix containing posterior mean of evaluations of coefficient function on training data.}
\item{beta.train}{Array of posterior samples of slope function evaluations. 
Only returned if \code{save_samples == TRUE}.}
\item{beta.test.mean}{Matrix containing posterior mean of evaluations of coefficient function on test data.}
\item{beta.test}{Array of posterior samples of slope function evaluations on test data. 
Only returned if \code{save_samples == TRUE}.}
\item{sigma}{Vector containing post-burnin samples of the residual standard deviation.}
\item{varcounts}{Array that counts the number of times a variable was used in a decision rule in each MCMC iteration. Structure is similar to that of \code{beta_train}, with rows corresponding to MCMC iteration and columns corresponding to predictors, with continuous predictors listed first followed by categorical predictors}
\item{timing}{Vector of runtimes for each chain}
}

\seealso{
\code{\link{probit_flexBART}} for binary outcomes.
}

\examples{
\dontrun{
## A modified version of Friedman's function (from the MARS paper)
# with 50 predictors in [-1,1]
set.seed(99)
mu_true <- function(df){
  # Recenter to [0,1]
  tmp_X <- (df+1)/2
  return(10*sin(pi*tmp_X[,1] * tmp_X[,2]) + 
           20 * (tmp_X[,8] - 0.5)^2 + 
           10 * tmp_X[,17] + 
           5 * tmp_X[,20])
}
n_train <- 1000
n_test <- 100
p_cont <- 50
sigma <- 1
train_data <- data.frame(Y = rep(NA, times = n_train))
for(j in 1:p_cont) train_data[[paste0("X",j)]] <- runif(n_train, min = -1, max = 1)
mu_train <- mu_true(train_data[,paste0("X",1:p_cont)])
train_data[,"Y"] <- mu_train + sigma * rnorm(n = n_train, mean = 0, sd = 1)

test_data <- data.frame(Y = rep(NA, times = n_test))
for(j in 1:p_cont) test_data[[paste0("X",j)]] <- runif(n_test, min = -1, max = 1)
mu_test <- mu_true(test_data[,paste0("X",1:p_cont)])
fit <-
  flexBART(formula = Y~bart(.),
           train_data = train_data,
           test_data = test_data,
           inform_sigma = TRUE, sparse = TRUE, 
           save_samples = FALSE)
par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(1,2))
plot(mu_train, fit$yhat.train.mean, ,
     pch = 16, cex = 0.5,
     xlab = "Actual", ylab = "Predicted", main = "Training")
abline(a = 0, b = 1, col = 'blue')
plot(mu_test, fit$yhat.test.mean, ,
     pch = 16, cex = 0.5,
     xlab = "Actual", ylab = "Predicted", main = "Testing")
abline(a = 0, b = 1, col = 'blue')
}
}