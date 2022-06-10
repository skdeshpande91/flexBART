\name{probit_flexBART}
\alias{probit_flexBART}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probit flexBART for binary outcomes
}
\description{
Implements flexBART with the Albert \& Chip (1993) data augmentation for probit regression.
}
\usage{
probit_flexBART(Y_train,
                X_cont_train = matrix(0, nrow = 1, ncol = 1),
                X_cat_train = matrix(0, nrow = 1, ncol = 1),
                X_cont_test = matrix(0, nrow = 1, ncol = 1),
                X_cat_test = matrix(0L, nrow = 1, ncol = 1),
                unif_cuts = rep(TRUE, times = ncol(X_cont_train)),
                cutpoints_list = NULL,
                cat_levels_list,
                sparse = FALSE,
                M = 200,
                mu0 = stats::qnorm(mean(Y_train)), tau = 1/sqrt(M),
                nd = 1000, burn = 1000, thin = 1,
                save_trees = TRUE, verbose = TRUE, 
                print_every = floor( (nd*thin + burn))/10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y_train}{Integer vector of binary (i.e. 0/1) responses for training data}
  \item{X_cont_train}{Matrix of continuous predictors for training data. Note, predictors must be re-scaled to lie in the interval [-1,1]. Default is a 1x1 matrix, which signals that there are no continuous predictors in the training data.}
  \item{X_cat_train}{Integer matrix of categorical predictors for training data. Note categorical levels should be 0-indexed. That is, if a categorical predictor has 10 levels, the values should run from 0 to 9. Default is a 1x1 matrix, which signals that there are no categorical predictors in the training data.}
  \item{X_cont_test}{Matrix of continuous predictors for testing data. Default is a 1x1 matrix, which signals that testing data is not provided.}
  \item{X_cat_test}{Integer matrix of categorical predictors for testing data. Default is a 1x1 matrix, which signals that testing data is not provided.}
  \item{unif_cuts}{Vector of logical values indicating whether cutpoints for each continuous predictor should be drawn from a continuous uniform distribution (\code{TRUE}) or a discrete set (\code{FALSE}) specified in \code{cutpoints_list}. Default is \code{TRUE} for each variable in \code{X_cont_train}}
  \item{cutpoints_list}{List of length \code{ncol(X_cont_train)} containing a vector of cutpoints for each continuous predictor. By default, this is set to \code{NULL} so that cutpoints are drawn uniformly from a continuous distribution.}
  \item{cat_levels_list}{List of length \code{ncol(X_cat_train)} containing a vector of levels for each categorical predictor. If the j-th categorical predictor contains L levels, \code{cat_levels_list[[j]]} should be the vector \code{0:(L-1)}. Default is \code{NULL}, which corresponds to the case that no categorical predictors are available.}
  \item{sparse}{Logical, indicating whether or not to perform variable selection based on a sparse Dirichlet prior rather than uniform prior; see Linero 2018. Default is \code{FALSE}}
  \item{M}{Number of trees in the ensemble. Default is 200.}
  \item{mu0}{Prior mean on the jumps/leaf parameters. See Details. Default is \code{qnorm(mean(Y_train))} so that the prior on the regression function is shrunk toward the empirical probability that Y = 1.}
  \item{tau}{Prior standard deviation on the jumps/leaf parameters. See Details. Default is \code{1/sqrt(M)}. Smaller values of \code{tau} induce stronger shrinkage towards the empirical probability that Y = 1.}
  \item{nd}{Number of posterior draws to return. Default is 1000.}
  \item{burn}{Number of MCMC iterations to be treated as "warmup" or "burn-in". Default is 1000.}
  \item{thin}{Number of post-warmup MCMC iteration by which to thin. Default is 1.}
  \item{save_trees}{Logical, indicating whether or not to save a text-based representation of the tree samples. This representation can be passed to \code{predict_flexBART} to make predictions at a later time. Default is \code{FALSE}.}
  \item{verbose}{Logical, inciating whether to print progress to R console. Default is \code{TRUE}.}
  \item{print_every}{As the MCMC runs, a message is printed every \code{print_every} iterations. Default is \code{floor( (nd*thin + burn)/10)} so that only 10 messages are printed.}
}
\details{
Implements the Albert & Chib (1993) data augmentation strategy for probit regression and models the regression function with a sum-of-trees.
The marginal prior of any evaluation of the regression function f(x) is a normal distribution centered at \code{mu0} with standard deviation \code{tau * sqrt(M)}. 
As such, for each x, the induced prior for P(Y = 1 | x) places 95\% probability on the interval \code{pnorm(mu0 -2 * tau * sqrt(M)), pnorm(mu0 + 2 * tau * sqrt(M))}. 
By default, we set \code{tau = 1/sqrt(M)} and \code{mu0 = qnorm(mean(Y_train))} to shrink towards the observed mean.
}
\value{
A list containing
\item{yhat_train}{Matrix with \code{nd} rows and \code{length(Y_train)} columns containing posterior samples of P(y = 1 | x) for the training data Each row corresponds to a posterior sample of the regression functionand each column corresponds to a training observation.}
\item{yhat_test}{If testing data was supplied, matrix containing posterior samples of P(y = 1 | x) for the testing data. Structure is similar to that of \code{yhat_train}}
\item{varcounts}{Matrix that counts the number of times a variable was used in a decision rule in each MCMC iteration. Structure is similar to that of \code{yhat_train}, with rows corresponding to MCMC iteration and columns corresponding to predictors}
\item{trees}{A list (or length \code{nd}) of character vectors (of lenght \code{M}) containing textual representations of the regression trees. These strings are parsed by \code{predict_flexBART} to reconstruct the C++ representations of the sampled trees.}

}
\references{
%% ~put references to the literature/web site here ~
Albert, J.H. and Chib, S. (1993), "Bayesian analysis of binary and polychotomous data". \emph{Journal of the American Statistical Association}. 88(422):669--679.
}