\name{probit_flexBART}
\alias{probit_flexBART}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Probit flexBART for binary outcomes
}
\description{
Fit a BART model of a binary responses using the the Albert \& Chib (1993) data augmentation for probit models.
}
\usage{
probit_flexBART(formula, train_data, test_data = NULL,...)
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
  \item{...}{additional arguments to be passed to the low-level model fitting function (see below)}
}
\details{
Implements the Albert & Chib (1993) data augmentation strategy for probit regression and models the regression function with a sum-of-trees.
The marginal prior of any evaluation of the regression function f(x) is a normal distribution centered at \code{mu0*M} with standard deviation \code{tau * sqrt(M)}. 
As such, for each x, the induced prior for P(Y = 1 | x) places 95\% probability on the interval \code{pnorm(mu0*M -2 * tau * sqrt(M)), pnorm(mu0*M + 2 * tau * sqrt(M))}. 
By default, we set \code{tau = 1/sqrt(M)} and \code{mu0 = qnorm(mean(Y_train))/M} to shrink towards the observed mean.
}
\value{
An object of \code{class} "flexBART" containing
\item{dinfo}{TO DO}
\item{trees}{A list (or length \code{nd}) of character vectors (of lenght \code{M}) containing textual representations of the regression trees. These strings are parsed by \code{\link{predict_flexBART}} to reconstruct the C++ representations of the sampled trees.}
\item{M}{number of trees used}
\item{cov_ensm}{TO DO}
\item{is.probit}{Logical with value \code{TRUE}. Used by \code{\link{predict.flexBART}}.}
\item{prob.train.mean}{Vector containing posterior mean of P(y = 1 | x) for the training data.}
\item{prob.train}{Matrix with \code{nd} rows and \code{length(Y_train)} columns containing posterior samples of P(y = 1 | x) for the training data Each row corresponds to a posterior sample of the regression functionand each column corresponds to a training observation. Only returned if \code{save_samples == TRUE}.}
\item{prob.test.mean}{Vector containing posterior mean of P(y = 1 | x) on testing data, if testing data is provided.}
\item{prob.test}{If testing data was supplied, matrix containing posterior samples of the regression function evaluated on the testing data. Structure is similar to that of \code{yhat_train}. Only returned if testing data is passed and \code{save_samples == TRUE}.}
\item{varcounts}{Matrix that counts the number of times a variable was used in a decision rule in each MCMC iteration. Structure is similar to that of \code{prob_train}, with rows corresponding to MCMC iteration and columns corresponding to predictors}
\item{timing}{Vector of runtimes for each chain}
}
\seealso{
\code{\link{flexBART}} for continuous outcomes.
}
\references{
%% ~put references to the literature/web site here ~
Albert, J.H. and Chib, S. (1993), "Bayesian analysis of binary and polychotomous data". \emph{Journal of the American Statistical Association}. 88(422):669--679.
}