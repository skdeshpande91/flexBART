\name{probit_networkBART}
\alias{probit_networkBART}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a BART model of a binary response with continuous predictors and network-indexed observations
}
\description{
Fit a BART model of a binary response with continuous predictors and network-indexed observations using the the Albert \& Chib (1993) data augmentation for probit models.
}
\usage{
probit_networkBART(Y_train,
                   vertex_id_train, # observation i belongs to node vertex_id_train[i].
                   X_cont_train = matrix(0, nrow = 1, ncol = 1),
                   vertex_id_test = NULL,
                   X_cont_test = matrix(0, nrow = 1, ncol = 1),
                   unif_cuts = rep(TRUE, times = ncol(X_cont_train)),
                   cutpoints_list = NULL,
                   A = matrix(0, nrow = 1, ncol = 1),
                   graph_split = TRUE,
                   graph_cut_type = 1,
                   sparse = FALSE,
                   M = 200,
                   mu0 = stats::qnorm(mean(Y_train)), tau = 1/sqrt(M),
                   nd = 1000, burn = 1000, thin = 1,
                   save_samples = TRUE,
                   save_trees = TRUE, verbose = TRUE, 
                   print_every = floor( (nd*thin + burn))/10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y_train}{Integer vector of binary (i.e. 0/1) responses for training data}
      \item{vertex_id_train}{An integer vector of the same length as \code{Y_train} that records the vertex label for each observation. Assumed to take values in the set {1, 2, ..., n_vertex} where n_vertex is the number of vertices in the network}
  \item{X_cont_train}{Matrix of continuous predictors for training data. Note, predictors must be re-scaled to lie in the interval [-1,1]. Default is a 1x1 matrix, which signals that there are no continuous predictors in the training data.}
  \item{vertex_id_test}{An integer vector recording the vertex labels for each observation in the testing data. Default is \code{NULL}.}
  \item{X_cont_test}{Matrix of continuous predictors for testing data. Default is a 1x1 matrix, which signals that testing data is not provided.}
  \item{unif_cuts}{Vector of logical values indicating whether cutpoints for each continuous predictor should be drawn from a continuous uniform distribution (\code{TRUE}) or a discrete set (\code{FALSE}) specified in \code{cutpoints_list}. Default is \code{TRUE} for each variable in \code{X_cont_train}}
  \item{cutpoints_list}{List of length \code{ncol(X_cont_train)} containing a vector of cutpoints for each continuous predictor. By default, this is set to \code{NULL} so that cutpoints are drawn uniformly from a continuous distribution.}
  \item{A}{Binary adjacency matrix of the network. Assumed to be symmetric.}
  \item{graph_split}{Logical, indicating whether to constrain the regression tree prior to respect adjacency when splitting on vertex labels. Default is \code{TRUE}.}
  \item{graph_cut_type}{An integer (0, 1, 2, or 3) that determines how partitions of the network are constructed. Default is 1. See Details of \code{\link{network_BART}}.}
  \item{sparse}{Logical, indicating whether or not to perform variable selection based on a sparse Dirichlet prior rather than uniform prior; see Linero 2018. Default is \code{FALSE}}
  \item{M}{Number of trees in the ensemble. Default is 200.}
  \item{mu0}{Prior mean on the jumps/leaf parameters. See Details. Default is \code{qnorm(mean(Y_train))} so that the prior on the regression function is shrunk toward the empirical probability that Y = 1.}
  \item{tau}{Prior standard deviation on the jumps/leaf parameters. See Details. Default is \code{1/sqrt(M)}. Smaller values of \code{tau} induce stronger shrinkage towards the empirical probability that Y = 1.}
  \item{nd}{Number of posterior draws to return. Default is 1000.}
  \item{burn}{Number of MCMC iterations to be treated as "warmup" or "burn-in". Default is 1000.}
  \item{thin}{Number of post-warmup MCMC iteration by which to thin. Default is 1.}
  \item{save_samples}{Logical, indicating whether to return all posterior samples. Default is \code{TRUE}. If \code{FALSE}, only posterior mean is returned.}
  \item{save_trees}{Logical, indicating whether or not to save a text-based representation of the tree samples. This representation can be passed to \code{predict_flexBART} to make predictions at a later time. Default is \code{FALSE}.}
  \item{verbose}{Logical, inciating whether to print progress to R console. Default is \code{TRUE}.}
  \item{print_every}{As the MCMC runs, a message is printed every \code{print_every} iterations. Default is \code{floor( (nd*thin + burn)/10)} so that only 10 messages are printed.}
}
\details{
Implements the Albert & Chib (1993) data augmentation strategy for probit regression and models the regression function with a sum-of-trees for network-indexed data.
The marginal prior of any evaluation of the regression function f(x) is a normal distribution centered at \code{mu0} with standard deviation \code{tau * sqrt(M)}. 
As such, for each x, the induced prior for P(Y = 1 | x) places 95\% probability on the interval \code{pnorm(mu0 -2 * tau * sqrt(M)), pnorm(mu0 + 2 * tau * sqrt(M))}. 
By default, we set \code{tau = 1/sqrt(M)} and \code{mu0 = qnorm(mean(Y_train))} to shrink towards the observed mean.

The flexBART prior is modified to handle network-structured categorical predictors as follows.
To split on the vertex label, we randomly draw a partition of the network into two (possibly disconnected) subgraphs.
See the documentation for \code{network_BART} for further details about how the partition is drawn according to the value of \code{graph_cut_type}.
}
\value{
A list containing
\item{prob.train.mean}{Vector containing posterior mean of P(y = 1 | x) for the training data.}
\item{prob.train}{Matrix with \code{nd} rows and \code{length(Y_train)} columns containing posterior samples of P(y = 1 | x) for the training data Each row corresponds to a posterior sample of the regression functionand each column corresponds to a training observation. Only returned if \code{save_samples == TRUE}.}
\item{prob.test.mean}{Vector containing posterior mean of P(y = 1 | x) on testing data, if testing data is provided.}
\item{prob.test}{If testing data was supplied, matrix containing posterior samples of the regression function evaluated on the testing data. Structure is similar to that of \code{yhat_train}. Only returned if testing data is passed and \code{save_samples == TRUE}.}
\item{varcounts}{Matrix that counts the number of times a variable was used in a decision rule in each MCMC iteration. Structure is similar to that of \code{prob_train}, with rows corresponding to MCMC iteration and columns corresponding to predictors}
\item{trees}{A list (or length \code{nd}) of character vectors (of lenght \code{M}) containing textual representations of the regression trees. These strings are parsed by \code{predict_flexBART} to reconstruct the C++ representations of the sampled trees.}
\item{diag}{A list containing diagnostic information about the Metropolis-Hastings step. Separate elements for numbers of tree proposals accepted and the number of axis-aligned and categorical rules proposed and rejected in each MCMC iteration.}
\item{is.probit}{Logical with value \code{TRUE}. Used by \code{predict_flexBART}.}
}
\references{
%% ~put references to the literature/web site here ~
Albert, J.H. and Chib, S. (1993), "Bayesian analysis of binary and polychotomous data". \emph{Journal of the American Statistical Association}. 88(422):669--679.
}