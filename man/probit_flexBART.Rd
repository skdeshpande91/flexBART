\name{probit_flexBART}
\alias{probit_flexBART}
\title{
Probit flexBART for binary outcomes
}
\description{
Fit a BART model of a binary responses using the the Albert and Chib (1993) data augmentation for probit models.
}
\usage{
probit_flexBART(formula, train_data, test_data = NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{an object of class \code{\link[stats:formula]{formula}} (or one that can be coerced to the class): 
  a symbolic description of the model to be fitted. 
  The details of model specification are given under 'Details'.}
  \item{train_data}{an object of class \code{data.frame} containing data used to train the model. 
  As usual, rows (resp. columns) correspond to observations (resp. variables)}
  \item{test_data}{an optional object of class \code{data.frame} containing test-set (i.e., out-of-sample) data. 
  Default is \code{NULL}.}
  \item{\dots}{Additional arguments for setting prior hyperparameters (e.g., number of trees, \eqn{\mu_{0}}{mu_0}, \eqn{\tau}{tau}, etc.) and MCMC control parameters (e.g., number of chains, iterations, etc.) 
  See \code{\link{flexBART}} for details about additional arguments.}
}
\details{
For a binary response \eqn{Y}, 
\eqn{p}{p} predictors \eqn{X_{1}, \ldots, X_{p}}{X_1, ..., X_p},
\code{probit_flexBART} models
\eqn{P(Y=1 \vert X = x) = \Phi(f(x)}{P(Y=1|x) = Phi(f(x))},
where \eqn{\Phi}{Phi} is the standard normal cumulative distribution function. 
\code{probit_flexBART} combines the Albert & Chib (1993) data augmentation strategy for probit regression
with the usual Bayesian backfitting used to fit (VC)BART models.

\subsection{The formula argument}{
\strong{Currently, \code{probit_flexBART} only supports fitting single ensemble probit BART models.}
So, the only valid formula will look something like \code{Y~bart(.)} or \code{Y ~ bart(x1+x2)}.
As with \code{\link{flexBART}}, you must include the string \dQuote{bart} on the right-hand side of the formula object.
}

\subsection{Prior specification}{
\code{probit_flexBART} approximates the function \eqn{f(x)}{f(x)} with an ensemble of binary regression trees.
It also specifies independent priors on the trees in the ensemble that are essentially identical to those deployed by \code{\link{flexBART}}.
That is, the tree structure is generated using a branching process in which
the probability that a node at depth \eqn{d}{d} is non-terminal is \eqn{\alpha \times (1 + d)^{-\beta}}{alpha * (1 + d)^(-beta)}.
Then, decision rules are drawn sequentially from the root down to each leaf.
Finally, independent \eqn{N(\mu_0, \tau^2)}{N(mu_0, tau^2)} priors are specified for the outputs in each leaf.

With this specification, the marginal prior of any evaluation of the regression function 
\eqn{f(x)}{f(x)} is \eqn{N(M \times \mu_{0}, \tau^{2} \times M)}{N(M * mu_0, M*tau^2)}, where \eqn{M}{M} is the number of trees in the ensemble.
Thus, for each \eqn{x}{x}, the induced prior for \eqn{P(Y = 1 \vert X = x)}{P(Y = 1|X = x)} places 95\% probability on the interval
\eqn{[\Phi(M \times \mu_{0} - 2 \times \tau \times \sqrt{M}), \Phi(M \times \mu_{0} + 2 \times \tau \times \sqrt{M})]}{[Phi(M*mu_0 - 2*tau*sqrt(M)),Phi(M*mu0 + 2 *tau*sqrt(M))]}.
By default, \code{probit_flexBART} sets \eqn{\tau = 1/\sqrt{M}}{tau=1/sqrt(M)} and \eqn{\Phi^{-1}(\overline{y})/M}{Phi^{-1}(bar{y})/M} 
(i.e., \code{qnorm(mean(Y))/M}).
Use the \code{mu0_vec} and \code{tau_vec} arguments to set other hyperparameter values.
}

}
\value{
An object of \code{class} \dQuote{flexBART} (essentially a list) containing
\item{dinfo}{Essentially a list containing information about the input and output variables. Used by \code{\link{predict.flexBART}}.}
\item{trees}{A list (or length \code{nd}) of character vectors (of lenght \code{M}) containing textual representations of the regression trees. 
These strings are parsed by \code{\link{predict.flexBART}} to reconstruct the C++ representations of the sampled trees.}
\item{scaling_info}{Essentially a list containing information for re-scaling raw MCMC output to the original outcome scale. Used by \code{\link{predict.flexBART}}.}
\item{M}{A copy of the argument \code{M_vec}. Used by \code{\link{predict.flexBART}}.}
\item{cov_ensm}{An \eqn{p \times R}{p x R} binary matrix encoding whose (j,r)-element is 1 if trees in the ensemble for \eqn{\beta_{r}(X)}{beta_r(X)} can split on \eqn{X_{j}}{X_j}.}
\item{prob.train.mean}{Vector containing posterior mean of \eqn{P(Y=1 \vert X = x)}{P(Y=1|X=x)} for the training data.}
\item{prob.train}{Matrix with \code{nd} rows and \code{length(Y_train)} columns containing posterior samples of \eqn{P(Y=1 \vert X = x)}{P(Y=1|X=x)} for the training data.
Each row corresponds to a posterior sample of the regression functionand each column corresponds to a training observation. 
Only returned if \code{save_samples = TRUE}.}
\item{prob.test.mean}{Vector containing posterior mean of \eqn{P(Y=1 \vert X = x)}{P(Y=1|X=x)} on testing data, if testing data is provided.}
\item{prob.test}{If testing data was supplied, matrix containing posterior samples of the regression function evaluated on the testing data. 
Structure is similar to that of \code{prob.train}. 
Only returned if testing data is passed and \code{save_samples = TRUE}.}
\item{varcounts}{Matrix that counts the number of times a variable was used in a decision rule in each MCMC iteration. 
Structure is similar to that of \code{prob.train}, with rows corresponding to MCMC iteration and columns corresponding to predictors}
\item{timing}{Vector of runtimes for each chain}
}
\seealso{
\code{\link{flexBART}} for continuous outcomes.
}
\references{
%% ~put references to the literature/web site here ~
Albert, J.H. and Chib, S. (1993) Bayesian analysis of binary and polychotomous data. 
\emph{Journal of the American Statistical Association}. \bold{88}(422):669--679.
\doi{doi:10.2307/2290350}.
}