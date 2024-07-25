\name{flexBART}
\alias{flexBART}

\title{
A more flexible BART
}
\description{
Implements Chipman et al. (2010)'s Bayesian additive regression trees (BART) method for nonparametric regression with continuous outcomes. The regression function is represented as a sum of binary regression trees. flexBART handles categorical outcomes more flexibly than other implementations of BART.
}
\usage{
flexBART(Y_train, 
         X_cont_train = matrix(0, nrow = 1, ncol = 1),
         X_cat_train = matrix(0L, nrow = 1, ncol = 1),
         X_cont_test = matrix(0, nrow = 1, ncol = 1),
         X_cat_test = matrix(0L, nrow = 1, ncol = 1),
         unif_cuts = rep(TRUE, times = ncol(X_cont_train)),
         cutpoints_list = NULL, cat_levels_list = NULL, 
         sparse = FALSE,
         M = 200,nd = 1000, burn = 1000, thin = 1, save_samples = TRUE,
         save_trees = TRUE, verbose = TRUE, print_every = floor( (nd*thin + burn))/10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y_train}{Vector of continuous responses for training data}
  \item{X_cont_train}{Matrix of continuous predictors for training data. Note, predictors must be re-scaled to lie in the interval [-1,1]. Default is a 1x1 matrix, which signals that there are no continuous predictors in the training data.}
  \item{X_cat_train}{Integer matrix of categorical predictors for training data. Note categorical levels should be 0-indexed. That is, if a categorical predictor has 10 levels, the values should run from 0 to 9. Default is a 1x1 matrix, which signals that there are no categorical predictors in the training data.}
  \item{X_cont_test}{Matrix of continuous predictors for testing data. Default is a 1x1 matrix, which signals that testing data is not provided.}
  \item{X_cat_test}{Integer matrix of categorical predictors for testing data. Default is a 1x1 matrix, which signals that testing data is not provided.}
  \item{unif_cuts}{Vector of logical values indicating whether cutpoints for each continuous predictor should be drawn from a continuous uniform distribution (\code{TRUE}) or a discrete set (\code{FALSE}) specified in \code{cutpoints_list}. Default is \code{TRUE} for each variable in \code{X_cont_train}}
  \item{cutpoints_list}{List of length \code{ncol(X_cont_train)} containing a vector of cutpoints for each continuous predictor. By default, this is set to \code{NULL} so that cutpoints are drawn uniformly from a continuous distribution.}
  \item{cat_levels_list}{List of length \code{ncol(X_cat_train)} containing a vector of levels for each categorical predictor. If the j-th categorical predictor contains L levels, \code{cat_levels_list[[j]]} should be the vector \code{0:(L-1)}. Default is \code{NULL}, which corresponds to the case that no categorical predictors are available.}
  \item{sparse}{Logical, indicating whether or not to perform variable selection based on a sparse Dirichlet prior rather than uniform prior; see Linero 2018. Default is \code{FALSE}}
  \item{M}{Number of trees in the ensemble. Default is 200.}
  \item{nd}{Number of posterior draws to return. Default is 1000.}
  \item{burn}{Number of MCMC iterations to be treated as "warmup" or "burn-in". Default is 1000.}
  \item{thin}{Number of post-warmup MCMC iteration by which to thin. Default is 1.}
  \item{save_samples}{Logical, indicating whether to return all posterior samples. Default is \code{TRUE}. If \code{FALSE}, only posterior mean is returned.}
  \item{save_trees}{Logical, indicating whether or not to save a text-based representation of the tree samples. This representation can be passed to \code{predict_flexBART} to make predictions at a later time. Default is \code{FALSE}.}
  \item{verbose}{Logical, inciating whether to print progress to R console. Default is \code{TRUE}.}
  \item{print_every}{As the MCMC runs, a message is printed every \code{print_every} iterations. Default is \code{floor( (nd*thin + burn)/10)} so that only 10 messages are printed.}
}
\details{
Default implementations of Bayesian Additive Regression Trees (BART) represent categorical predictors using several binary indicators, one for each level of each categorical predictor. Axis-aligned decision rules are well-defined with these indicators; they send one level of a categorical predictor to the left and all other levels to the right (or vice versa). Regression trees built with these rules partition the set of all levels of a categorical predictor by recursively removing one level at a time. Unfortunately, most partitions of the levels cannot be built with this ``remove one at a time'' strategy, meaning that default implementations of BART are extremely limited in their ability to ``borrow strength'' across groups of levels.

\code{flexBART} overcomes this limitation using a new prior on regression trees. Under this new prior, conditional on splitting on a categorical predictor at a particular node in the tree, levels of the predictor are sent to the left and right child uniformly at random. In this way, multiple levels of a categorical predictor are able to be clustered together.
}
\value{
A list containing
\item{y_mean}{Mean of the training observations (needed by \code{predict_flexBART})}
\item{y_sd}{Standard deviation of the training observations (needed by \code{predict_flexBART})}
\item{yhat.train.mean}{Vector containing posterior mean of evaluations of regression function on training data.}
\item{yhat.train}{Matrix with \code{nd} rows and \code{length(Y_train)} columns. Each row corresponds to a posterior sample of the regression function and each column corresponds to a training observation. Only returned if \code{save_samples == TRUE}.}
\item{yhat.test.mean}{Vector containing posterior mean of evaluations of regression function on testing data, if testing data is provided.}
\item{yhat.test}{If testing data was supplied, matrix containing posterior samples of the regression function evaluated on the testing data. Structure is similar to that of \code{yhat_train}. Only returned if testing data is passed and \code{save_samples == TRUE}.}
\item{sigma}{Vector containing ALL samples of the residual standard deviation, including burn-in.}
\item{varcounts}{Matrix that counts the number of times a variable was used in a decision rule in each MCMC iteration. Structure is similar to that of \code{yhat_train}, with rows corresponding to MCMC iteration and columns corresponding to predictors, with continuous predictors listed first followed by categorical predictors}
\item{diag}{A list containing diagnostic information about the Metropolis-Hastings step. Separate elements for numbers of tree proposals accepted and the number of axis-aligned and categorical rules proposed and rejected in each MCMC iteration.}
\item{trees}{A list (or length \code{nd}) of character vectors (of lenght \code{M}) containing textual representations of the regression trees. These strings are parsed by \code{predict_flexBART} to reconstruct the C++ representations of the sampled trees.}
}

\seealso{
\code{\link{probit_flexBART}} for binary outcomes and \code{\link{network_BART}} for regression with continuous outcomes on a network.
}

\examples{
## Friedman's function (from the MARS paper)
# Slight modification: flexBART requires all X's to be in [-1,1]
# Friendman's function
f0 <- function(x_cont, x_cat){
  if(!all(abs(x_cont) <= 1)){
    stop("all entries in x_cont must be between -1 and 1")
  } else{
    x <- (x_cont+1)/2 # convert to [0,1]
    return(10 * sin(pi*x[,1]*x[,2]) + 20 * (x[,3] - 0.5)^2 + 10*x[,4] + 5 * x[,5])
  }
}
f1 <- function(x_cont, x_cat){
  if(!all(abs(x_cont) <= 1)){
    stop("all entries in x_cont must be between -1 and 1")
  } else{
    x <- (x_cont+1)/2
    z1 <- x[,1]
    z2 <- 1*(x[,2] > 0.5)
    return(3 * z1 + (2 - 5 * z2) * sin(pi * z1) - 2 * z2)
  }
}
f <- function(x_cont, x_cat){
    mu <- rep(NA, times = nrow(x_cont))
    tmp_f0 <- f0(x_cont, x_cat)
    tmp_f1 <- f1(x_cont, x_cat)
    
    index0 <- which(x_cat[,1] \%in\% c(0,1,3, 5,6, 9))
    index1 <- which(x_cat[,1] \%in\% c(2,4,7,8))
    mu[index0] <- tmp_f0[index0]
    mu[index1] <- tmp_f1[index1]
    return(mu)
}
sigma <- 0.5
n <- 2000
n_test <- 100
p_cont <- 10
p_cat <- 10
set.seed(99)
X_cont <- matrix(runif(n*p_cont, min = -1, max = 1), nrow = n, ncol = p_cont)
X_cat <- matrix(sample(0:9, size = n*p_cat, replace = TRUE), nrow = n, ncol = p_cat)

mu <- f(X_cont, X_cat)
y <- mu + sigma * rnorm(n, mean = 0, sd = 1)

# In this example, all categorical variables take values in \{0, 1, ..., 9\}
cat_levels_list <- list()
for(j in 1:p_cat) cat_levels_list[[j]] <- 0:9

# Token run that ensures installation was successful
set.seed(99)
flexBART_fit <- flexBART(Y_train = y, X_cont_train = X_cont, X_cat_train = X_cat,
                         cat_levels_list = cat_levels_list,
                         nd = 5, burn = 5)
\dontrun{
X_cont_test <- matrix(runif(n_test*p_cont, min = -1, max = 1), nrow = n_test, ncol = p_cont)
X_cat_test <- matrix(sample(0:9, size = n_test*p_cat, replace = TRUE), nrow = n_test, ncol = p_cat)
mu_test <- f(X_cont_test, X_cat_test)

set.seed(99)
flexBART_fit <- flexBART(Y_train = y, X_cont_train = X_cont, X_cat_train = X_cat,
                         X_cont_test = X_cont_test, X_cat_test = X_cat_test,
                         cat_levels_list = cat_levels_list)

par(mar = c(3,3,2,1), mgp = c(1.8, 0.5, 0), mfrow = c(1,2))
plot(mu, flexBART_fit$yhat.train.mean, pch = 16, cex = 0.5, 
     xlab = "Truth", ylab = "Estimate", main = "Training")
abline(a = 0, b = 1, col = 'blue')
plot(mu_test, flexBART_fit$yhat.test.mean, pch = 16, cex = 0.5, 
     xlab = "Truth", ylab = "Estimate", main = "Testing")
abline(a = 0, b = 1, col = 'blue')

}
}