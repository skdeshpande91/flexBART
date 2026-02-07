\name{flexBART}
\alias{flexBART}

\title{
A more flexible BART
}
\description{
Implements a more flexible version of Bayesian Additive Regression Trees (BART) and Deshpande et al. (2024)'s varying coefficient BART model (VCBART) that, at a high-level, represents regression functions as the sum of several weighted ensembles of binary regression trees. 
Through the formula interface, users can carefully control the splitting variables used in each ensemble. 
The optional \code{family} argument allows users to fit generalized linear varying coefficient models (e.g., logistic regression models in which the covariate effects vary as functions of effect modifiers).
Conditionally Gaussian heteroskedastic regression models are also supported. 
This implementation includes several improvements on Deshpande (2025)'s priors for decision rules for categorical predictors. 
Trees directly partition the levels of categorical predictors, enabling much more flexible \dQuote{partial pooling} across groups/categories (esp., when they are network-structured or nested) than one-hot encoding.
}
\usage{
flexBART(formula, train_data, 
         test_data = NULL, 
         initialize_sigma = TRUE, ...)
}
\arguments{
  \item{formula}{an object of class \code{\link[stats:formula]{formula}} (or one that can be coerced to the class): a symbolic description of the model to be fitted. The details of model specification are given under \sQuote{Details}.}
  \item{train_data}{an object of class \code{data.frame} containing data used to train the model. 
  As usual, rows (resp. columns) correspond to observations (resp. predictors)}
  \item{test_data}{an optional object of class \code{data.frame} containing test-set (i.e., out-of-sample) data. 
  Default is \code{NULL}.}
  \item{initialize_sigma}{a logical value for whether the residual variance 
  should be initialized using regularized regression (\code{inform_sigma = TRUE}) or not (\code{inform_sigma = FALSE}). 
  Default is \code{TRUE}. When \code{inform_sigma = FALSE}, residual variance initialized to marginal variance of the outcome/response. See \sQuote{Details} below.}
  \item{\dots}{Additional arguments for specifying error distribution (\code{family}) and setting prior hyperparameters and MCMC control parameters (e.g., number of chains, iterations, etc.). 
  See \sQuote{Details} below.}
}
\details{
\code{flexBART} allows users to specify and fit several generalized linear varying coefficient models using ensembles of binary regression trees.
Given variables \eqn{Y}, \eqn{X_{1}, \ldots, X_{p}}{X_1, ..., X_p}, and \eqn{Z_{1}, \ldots, Z_{R}}{Z_1, ..., Z_R} the model asserts that

\eqn{g(\mathbb{E}[Y \vert X,Z]) = \beta_{0}(X) + \beta_{1}(X)Z_{1} + \cdots + \beta_{R}(X)Z_{R}},

for a known link function \eqn{g}.
\code{flexBART} approximates each function \eqn{\beta_{r}(X)}{beta_r(X)} with its own regression tree ensemble and returns posterior samples of all R tree ensembles.

\subsection{The formula argument}{

Models for \code{flexBART} are specified symbolically. 
A simple, single ensemble model has the form \code{response ~ bart(terms)} where \code{response} is the numeric response and \code{terms} specifies the predictors on which the trees in the ensemble
are allowed to split.
A \code{terms} specification of the form \code{.} will include all variables \emph{except} \code{response} that are found in \code{train_data}.
A \code{terms} specification of the form \code{x1 + x2} will include just the variables \code{x1} and \code{x2}.
A \code{terms} specification of the form \code{.-x1-x2} will include all variables found in \code{train_data} \emph{except} \code{response}, \code{x1}, and \code{x2}.
Of course, you can specify more than two variables in \code{terms}.

To fit a varying coefficient model, the formula argument takes the form \code{response ~ z1 * bart(terms1) + z2 * bart(terms2) + \dots} where \code{response} is the numeric response and 
\code{z1}, \code{z2}, etc. are numeric covariates and \code{terms1}, \code{terms2} are like above.
For varying coefficient models, if you specify terms using \code{.}, then the trees are not able to split on the covariates used to weight the different ensembles (i.e., \code{z1},\code{z2}, \dots).
To allow trees to split on these outer weighting covariates, you must manually specify them inside of the terms.

You \strong{must} include the string \dQuote{bart} on the right-hand side of the supplied formula. 
Expressions like \code{Y~Z1+Z2} or \code{Y~Z1 + Z2*bart(.)} will not work.

\subsection{Heteroskedastic Regression}{

\code{flexBART} supports fitting models of the form \eqn{Y \sim N(\mu(X), \sigma^{2}(X))}{Y ~ N(mu(X), sigma^2(X))} in which both \eqn{\mu(X)}{mu(X)} and \eqn{\log \sigma(X)}{log(sigma(X))} are approximated with regression tree ensembles.
Such a model can be specified through the \code{formula} argument as \code{Y~bart(terms1)+sigma(terms2)} where \code{terms1} and \code{terms2} are like above. 
By default, \code{flexBART} fits a regression model with homoskedastic Gaussian errors. 
}


}

\subsection{Supported families}{

\code{flexBART} currently supports models with continuous and binary outcomes.
For continuous outcomes, \code{flexBART} fits regression models with conditionally Gaussian errors and for binary outcomes, \code{flexBART} can fit logistic and probit regression models.
Users can specify the appropriate error distribution using the optional arguments
\itemize{
  \item \code{family}: a family object, or the result of a call to a family function. Currently only \dQuote{gaussian} and \dQuote{binomial} are supported. Default is \dQuote{gaussian}.
}
}

\subsection{Numerical predictors}{
Decision rules based on a numerical predictor \eqn{X}{X} take the form \eqn{\{X < c\}}{X<c}.
The \dQuote{cutpoint} \eqn{c}{c} in such a decision rule at a given node in a regression tree is drawn uniformly from a set of available values of \eqn{X}{X}.
This set of available values is determined by the decision rules at the node's ancestors in the tree.
This set can be a discrete set of continuous interval, depending on the nature of \eqn{X}{X}.

\code{flexBART()} initially determines whether \eqn{X}{X} is discrete (i.e., ordinal) or continuous by first sorting the unique values of \eqn{X}{X} and computing consecutive differences.
If there are fewer than (resp. more than) \code{n_unik_diffs} unique consecutive differences, \eqn{X}{X} is treated as a discrete (resp., continuous) predictor.
For discrete, numerical \eqn{X}{X}, the initial set of available cutpoints is set of unique values of \eqn{X}{X}.
\code{flexBART()} rescales continuous \eqn{X}{X} values to the interval [-1,1], which forms the initial set of available cutpoint values. 
When re-scaling a continuous \eqn{X}{X}, \code{flexBART()} \code{flexBART()} adds (resp. subtracts) \code{pad} standard deviations to the maximal (resp. minimal) value of \eqn{X}{X}.
The maximal and minimal values and the standard deviation of \eqn{X}{X} is determined using training data (i.e., \code{train_data}) and, if provided, the testing data (i.e., \code{test_data}).

The amount of padding and the cut-off for 

\itemize{
  \item \code{n_unik_diffs}: Threshold for the number of unique consecutive differences for a numerical predictor to be considered discrete. Default is 5.
  \item \code{pad}: Number of standard deviations with which to pad range of continuous predictors before internally re-scaling to [-1,1]. Default is 0.2. 
}
}

\subsection{Categorical predictors}{
Many implementations of BART and its extensions represent categorical predictors using several binary indicators, one for each level of each categorical predictor. 
Axis-aligned decision rules are well-defined with these indicators: they send one level of a categorical predictor to the left and all other levels to the right (or vice versa). 
Regression trees built with these rules partition the set of all levels of a categorical predictor by recursively removing one level at a time. Unfortunately, most partitions of 
the levels cannot be built with this \dQuote{remove one at a time} strategy, meaning these implementations are extremely limited in their ability to "borrow strength'' across groups of levels.

\code{flexBART()} overcomes this limitation using a new prior on regression trees.
Under this new prior, conditional on splitting on a categorical predictor at a particular node in the tree, levels of the predictor are sent to the left and right child uniformly at random. 
In this way, multiple levels of a categorical predictor are able to be clustered together.
\code{flexBART()} implements several decision rule priors specifically for nested categorical predictors and network-structured categorical predictors.

\subsection{Decision rules for network-structured predictors}{
To expose the network structure between levels of categorical predictor, users should specify the \code{adjacency_list} argument.
This argument should be a \strong{named} \code{list} with one element per network-structured predictor.
Each element should be a binary or weighted adjacency matrix whose row and column names correspond to the unique values of the corresponding predictor. 
\code{flexBART()} implements four different priors over decision rules for network-structured predictors.
Each prior recursively partitions the network into two pieces.
The argument \code{graph_cut_type} determines which prior is implemented:
\itemize{
  \item \code{graph_cut_type = 1}: A deterministic partition is formed based on the signs of entries in the second smallest eigenvector of the unweighted graph Laplacian (i.e., the Fiedler vector).
  \item \code{graph_cut_type = 2}: The partition is formed by first drawing a uniformly random spanning tree using Wilson's algorithm and then deleting a uniformly random edge from the tree.
  \item \code{graph_cut_type = 3}: Like \code{graph_cut_type=2} but the probability of deleting a spanning tree edge is proportional to the size of the smallest cluster that results from that deletion.
  \item \code{graph_cut_type = 4}: The partition is formed by first drawing a uniformly random spanning tree and computing the Fiedler vector of the spanning tree
}

If unspecified, the default value of \code{graph_cut_type} is 2. 
If no adjacency information is provided, \code{graph_cut_type} is ignored.
}
\subsection{Decision rules for nested predictors}{
\code{flexBART()} can automatically detect potential nesting structure between categorical predictors.
It implements eight different decision rule priors, which are based on the arguments \code{nest_v}, \code{nest_v_option}, and \code{nest_c}.

The Boolean argument \code{nest_v} indicates whether nesting structure is used when selecting the splitting variable (\code{nest_v = TRUE}) or not (\code{nest_v = FALSE}).
When \code{nest_v = TRUE}, the argument \code{nest_v_option} determines how \code{flexBART()} selects a splitting variable.
Say it is trying to draw a decision rule at a tree node whose ancestor splits on a nested predictor \eqn{X_v}{X_v}.
In addition to predictors that are not (i) nested within \eqn{X_v}{X_v} and (ii) nest \eqn{X_v}{X_v}, \code{flexBART()} places positive probability on splitting on
\itemize{
  \item \code{nest_v_option = 0}: \eqn{X_v}{X_v} but not variables that are nested within or that nest \eqn{X_v}{X_v}.
  \item \code{nest_v_option = 1}: \eqn{X_v}{X_v} and variables that are nested within \eqn{X_v}{X_v} but not variables that nest \eqn{X_v}{X_v}.
  \item \code{nest_v_option = 2}: \eqn{X_v}{X_v} and variables that nest \eqn{X_v}{X_v} but not variables nested within \eqn{X_v}{X_v}.
  \item \code{nest_v_option = 3}: \eqn{X_v}{X_v} and variables that nest or are nested within \eqn{X_v}{X_v}.
}
The Boolean argument \code{nest_c} controls whether nesting structure is used when selecting the set of levels assigned to the left branch of 
a decision node (\code{nest_c = TRUE}) or not (\code{nest_c = FALSE}).
Default is \code{nest_c = TRUE}.
This argument is ignored when nested structure is not detected.
}
}

\subsection{Prior specification and standardization}{
Internally, \code{flexBART()} re-centers and re-scales continuous outcomes Y to have mean 0 and standard deviation 1.
Except when \eqn{Z_{r}}{Z_r} is all ones (i.e., for intercept terms), \code{flexBART()} also standardizes the covariates to have mean 0 and standard deviation 1.
It then places independent BART priors on the coefficients for the varying coefficient model on the standardized scale. 

No transformations are applied when \code{family} is \dQuote{binomial}.
}

\subsection{Regression tree priors}{
\code{flexBART()} specifies independent priors on each tree in the ensemble approximating \eqn{\beta_{r}(x)}{beta_r(x)}.
Under this prior, the tree structure is generated using a branching process in which
the probability that a node at depth \eqn{d}{d} is non-terminal is \eqn{\alpha \times (1 + d)^{-\beta}}{alpha * (1 + d)^(-beta)}.
Then, decision rules are drawn sequentially from the root down to each leaf.
Finally, independent \eqn{N(\mu_0, \tau^2)}{N(mu_0, tau^2)} priors are specified for the outputs in each leaf.

Users can specify prior hyperparameters for each ensemble using the following optional arguments (passed through \dots):
\itemize{
  \item \code{M_vec}: Vector of number of trees used in each ensemble. Default is 50 trees for each ensemble.
  \item \code{alpha_vec}: Vector of base parameter \eqn{\alpha}{alpha} in branching process tree prior for each ensemble. Default is 0.95 for all ensembles.
  \item \code{beta_vec}: Vector of power parameter \eqn{\beta}{beta} in the branching process tree prior for each ensemble. Default is 2 for all ensembles.
  \item \code{mu0_vec}: Vector of prior means \eqn{\mu_0}{mu_0} at each leaf in each ensemble. Default is 0 for all ensembles.
  \item \code{tau_vec}: Vector of prior standard deviations in each leaf \eqn{\tau}{tau} in each ensemble. 
  When there is only one ensemble, default value is \code{y_range/(2 * 2 * sqrt(M_vec[1])} where \code{y_range} is the range of the standardized response. 
  When there are multiple ensembles, default is \code{1/sqrt(M_vec)} for all ensembles.
}
}
\subsection{Prior for the residual variance}{
\code{flexBART()} specifies an \eqn{Inv.\ Gamma(\nu/2, \nu*\lambda/2)}{Inv. Gamma(nu/2, nu*lambda/2)} prior on the residual variance \emph{on the standardized outcome scale}.
The prior is calibrated so that this prior places a user-specified amount of prior probability on the event that this residual variance 
is less than some initial over-estimate \eqn{\hat{\sigma}^2_0}{hat(sigma)^2_0}.
Users can control this prior with the optional arguments
\itemize{
  \item \code{nu}: Degrees of freedom \eqn{\nu}{nu} for inverse gamma prior on the residual variance (on the standardized scale). Default is 3.
  \item \code{sigest}: Initial over-estimate of the residual variance (\eqn{\hat{\sigma}^2_0}{hat(sigma)^2_0}).
  If not provided, it is set based on a fitted l1-regularized regression model (if \code{initialize_sigma = TRUE}) or 
  to the variance of the outcome variable (if \code{initialize_sigma = FALSE}).
  \item \code{sigquant}: Amount of prior probability on the event that residual variance is less than initial over-estimate 
  (i.e., \eqn{\sigma < \hat{\sigma}_0}{sigma < hat(sigma)_0}). Default is 0.9. 
}

}

\subsection{Saving sampled trees and function evaluations}{

The arguments \code{save_samples} and \code{save_trees} respectively control the
amount of output returned by \code{flexBART()}.

\itemize{
  \item \code{save_samples}: Logical, indicating whether to return all posterior samples. 
  Default is \code{TRUE}. If \code{FALSE}, only posterior mean is returned. 
  \item \code{save_trees}: Logical, indicating whether or not to save a text-based representation of the tree samples. 
  This representation can be passed to \code{\link{predict.flexBART}} to make predictions on other data and/or in another R session. 
  Default is \code{TRUE}.
}

When \code{save_samples = TRUE}, \code{flexBART()} internally creates, populates, and 
outputs an array containing all post-"burn-in" posterior draws of each function 
specified in the \code{formula} evaluated at every training and testing observation.
Storing this array can require considerable memory if the number of training or testing samples is large.
For this reason, when there are more than 10,000 observations across the training and testing datasets,
it is \strong{highly} recommended to set \code{save_samples = FALSE} and \code{save_trees = TRUE}.
With these settings, \code{flexBART()} just returns the posterior mean of the function evaluations 
for each observation and a text-based representation of the regression trees.
To access posterior samples for individual observations, pass the fitted object 
returned by \code{flexBART} to \code{\link{predict.flexBART}} along with a data 
frame containing the values of the predictors.
See the package vignettes for examples of this workflow.
}

\subsection{Additional arguments}{
The following arguments, which are passed to internal pre-processing and model fitting functions, 
can be supplied using the \code{\dots}:

\itemize{
  \item \code{sparse}: Whether to perform variable selection based on the sparse Dirichlet 
  prior rather than uniform (see Linero (2018)). Default is \code{TRUE} but is ignored if \code{nest_v = TRUE}. 
  \item \code{a_u}: Hyper-parameter for the \eqn{Beta(a_u,b_u)}{Beta(a_u,b_u)} hyper-prior 
  used to specify Linero (2018)'s sparse Dirichlet prior. Default is 0.5. Ignored if \code{nest_v=TRUE}.
  \item \code{b_u}: Hyper-parameter for the \eqn{Beta(a_u,b_u)}{Beta(a_u,b_u)} hyper-prior
  used to specify Linero (2018)'s sparse Dirichlet prior. Default is 1. Ignored if \code{nest_v=TRUE}.
  \item \code{n.chains}: Number of MCMC chains to run. Default is 4.
  \item \code{nd}: Number of posterior draws to return per chain. Default is 1000.
  \item \code{burn}: Number of "burn-in" or "warm-up" iterations per chain. Default is 1000.
  \item \code{thin}: Number of post-warmup MCMC iterations by which to thin. Default is 1.
  \item \code{verbose}: Logical, indicating whether to print progress to the \R console. Default is \code{TRUE}.
  \item \code{print_every}: As the MCMC runs, a message is printed every \code{print_every} iterations. 
  Default is \code{floor( (nd*thin + burn)/10)} so that only 10 messages are printed.
}
}
}
\value{
An object of \code{class} \dQuote{flexBART} (essentially a list) containing
\item{dinfo}{Essentially a list containing information about the input and output variables. Used by \code{\link{predict.flexBART}}.}
\item{trees}{A list (or length \code{nd}) of character vectors (of length \code{M}) containing textual representations of the regression trees. These strings are parsed by \code{\link{predict.flexBART}} to reconstruct the C++ representations of the sampled trees.}
\item{scaling_info}{Essentially a list containing information for re-scaling raw MCMC output to the original outcome scale. Used by \code{\link{predict.flexBART}}.}
\item{M}{A copy of the argument \code{M_vec}. Used by \code{\link{predict.flexBART}}.}
\item{family}{Records the \code{family} argument passed when the model was fit.}
\item{link}{Records the link function specified in the \code{family} object passed when the model was fit.}
\item{heteroskedastic}{A Boolean identifying whether or not there was a \code{sigma()} ensemble in the model.}
\item{cov_ensm}{An \eqn{p \times R}{p x R} binary matrix encoding whose (j,r)-element is 1 if trees in the ensemble for \eqn{\beta_{r}(X)}{beta_r(X)} can split on \eqn{X_{j}}{X_j}.}
\item{yhat.train.mean}{Vector containing posterior mean of evaluations of regression function for each observation in the training data.}
\item{yhat.train}{Matrix with \code{nd} rows and \eqn{n}{n} columns. Each row corresponds to a posterior sample of the regression function and each column corresponds to a training observation. Only returned if \code{save_samples = TRUE}.}
\item{yhat.test.mean}{Vector containing posterior mean of evaluations of regression function on testing data, if testing data is provided.}
\item{yhat.test}{If testing data was supplied, matrix containing posterior samples of the regression function evaluated on the testing data. Structure is similar to that of \code{yhat.train}. Only returned if testing data is passed and \code{save_samples = TRUE}.}
\item{beta.train.mean}{Matrix containing posterior mean of evaluations of coefficient function on training data.}
\item{beta.train}{Array of posterior samples of slope function evaluations. 
Only returned if \code{save_samples = TRUE}.}
\item{beta.test.mean}{Matrix containing posterior mean of evaluations of coefficient function on test data.}
\item{beta.test}{Array of posterior samples of slope function evaluations on test data. 
Only returned if \code{save_samples = TRUE}.}
\item{sigma}{Vector containing post-burnin samples of the residual standard deviation.}
\item{varcounts}{Array that counts the number of times a variable was used in a decision rule in each MCMC iteration. Structure is similar to that of \code{beta_train}, with rows corresponding to MCMC iteration and columns corresponding to predictors, with continuous predictors listed first followed by categorical predictors}
\item{timing}{Vector of runtimes for each chain.}
}


\examples{

## A modified version of Friedman's function with 50 predictors in [-1,1]
set.seed(99)
mu_true <- function(df){
  # Recenter to [0,1]
  tmp_X <- (df+1)/2
  return(10*sin(pi*tmp_X[,1] * tmp_X[,2]) + 
           20 * (tmp_X[,8] - 0.5)^2 + 
           10 * tmp_X[,17] + 
           5 * tmp_X[,20])
}
## Set problem dimensions
n_train <- 500
p_cont <- 50

## Set residual error variance
sigma <- 1

## Generate training data
train_data <- data.frame(Y = rep(NA, times = n_train))
for(j in 1:p_cont) train_data[[paste0("X",j)]] <- runif(n_train, min = -1, max = 1)
mu_train <- mu_true(train_data[,paste0("X",1:p_cont)])
train_data[,"Y"] <- mu_train + sigma * rnorm(n = n_train, mean = 0, sd = 1)


## Fit flexBART model
fit <-
  flexBART(formula = Y~bart(.),
           train_data = train_data)
           
\donttest{
  ## Plot the posterior mean regression function evaluations (i.e., fitted values)
  ## against the actual values. The points should cluster around the 45-degree diagonal
  ## line y=x
  plot(mu_train, fit$yhat.train.mean, ,
       pch = 16, cex = 0.5,
       xlab = "Actual", ylab = "Predicted", main = "Training")
  abline(a = 0, b = 1, col = 'blue')
}
}
\references{
Chipman, H, George, E.I., and McCulloch, R. (2008) BART: Bayesian Additive Regression Trees. 
\emph{Annals of Applied Statistics}. \bold{4}(1):266--298. \doi{10.1214/09-AOAS285}.

Deshpande, S.K. (2025) \bold{flexBART:} Flexible Bayesian regression trees with categorical predictors. 
\emph{Journal of Computational and Graphical Statistics.} \bold{34}(3):1117--1126. \doi{doi:10.1080/10618600.2024.2431072}.

Deshpande, S.K., Bai, R., Balocchi, C., Starling, J.E., and Weiss, J. (2024). VCBART:Bayesian trees for varying coefficients. 
\emph{Bayesian Analysis}. \doi{doi:10.1214/24-BA1470}.

Linero, A.R. (2018) Bayesian regression trees for high-dimensional prediction and variable selection. 
\emph{Journal of the American Statistical Association.} \bold{113}(522):626--636. 
\doi{doi:10.1080/01621459.2016.1264957}.

Pratola, M.T., Chipman, H.A., George, E.I., and McCulloch R.E. (2020). Heteroskedastic BART via multiplicative regression trees.
\emph{Journal of Computational and Graphical Statistics.} \bold{29}(2):405--417. 
\doi{doi:10.1080/10618600.2019.1677243}.
}