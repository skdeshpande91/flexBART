\name{network_BART}
\alias{network_BART}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a BART model of a continuous response with continuous predictors and network-indexed observations
}
\description{
In some applications, the predictor-outcome pairs (x,y) are network-indexed. That is, the data may be observed at the vertices of a fixed network. This function extends BART to treat location in the network as an additional predictor.
The vertex labels are passed as categorical predictors and the regression tree prior creates contiguous splits on that variable that respect the underlying adjacency. See Details for more. 
}
\usage{
network_BART(Y_train,
             vertex_id_train,
             X_cont_train = matrix(0, nrow = 1, ncol = 1),
             vertex_id_test = NULL,
             X_cont_test = matrix(0, nrow = 1, ncol = 1),
             unif_cuts = rep(TRUE, times = ncol(X_cont_train)),
             cutpoints_list = NULL,
             A = matrix(0, nrow = 1, ncol = 1),
             graph_split = TRUE,
             graph_cut_type = 0,
             sparse = FALSE,
             M = 200,nd = 1000, burn = 1000, thin = 1,
             save_samples = TRUE,
             save_trees = TRUE, verbose = TRUE, print_every = floor( (nd*thin + burn))/10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y_train}{Vector of continous responses for training data}
  \item{vertex_id_train}{An integer vector of the same length as \code{Y_train} that records the vertex label for each observation. Assumed to take values in the set {1, 2, ..., n_vertex} where n_vertex is the number of vertices in the network}
  \item{X_cont_train}{Matrix of continuous predictors for training data. Note, predictors must be re-scaled to lie in the interval [-1,1]. Default is a 1x1 matrix, which signals that there are no continuous predictors in the training data.}
  \item{vertex_id_test}{An integer vector recording the vertex labels for each observation in the testing data. Default is \code{NULL}.}
  \item{X_cont_test}{Matrix of continuous predictors for testing data. Default is a 1x1 matrix, which signals that testing data is not provided.}
  \item{unif_cuts}{Vector of logical values indicating whether cutpoints for each continuous predictor should be drawn from a continuous uniform distribution (\code{TRUE}) or a discrete set (\code{FALSE}) specified in \code{cutpoints_list}. Default is \code{TRUE} for each variable in \code{X_cont_train}}
  \item{cutpoints_list}{List of length \code{ncol(X_cont_train)} containing a vector of cutpoints for each continuous predictor. By default, this is set to \code{NULL} so that cutpoints are drawn uniformly from a continuous distribution.}
  \item{A}{Binary adjacency matrix of the network. Assumed to be symmetric.}
  \item{graph_split}{Logical, indicating whether to constrain the regression tree prior to respect adjacency when splitting on vertex labels. Default is \code{TRUE}.}
  \item{graph_cut_type}{An integer (0 or 1) that determines how partitions of the network are constructed. See Details.}
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
The flexBART prior is modified to handle network-structured categorical predictors as follows.
To split on the vertex label, we randomly draw a partition of the network into two connected subgraphs.
When \code{graph_cut_type == 0}, the partition is drawn by first drawing a uniform random spanning tree using Wilson's algorithm and then deleting a uniformly selected edge from the tree. When \code{graph_cut_type == 1}, the partition is drawn by first drawing uniform weights for each edge, computing the minimum spanning tree with Boruvka's algorithm, and deleting an edge uniformly at random from the MST. In this way, each regression tree can partition the network in connected subgraphs. 
}
\value{
A list containing
\item{y_mean}{Mean of the training observations (needed by \code{predict_flexBART})}
\item{y_sd}{Standard deviation of the training observations (needed by \code{predict_flexBART})}
item{yhat.train.mean}{Vector containing posterior mean of evaluations of regression function on training data.}
\item{yhat.train}{Matrix with \code{nd} rows and \code{length(Y_train)} columns. Each row corresponds to a posterior sample of the regression function and each column corresponds to a training observation. Only returned if \code{save_samples == TRUE}.}
\item{yhat.test.mean}{Vector containing posterior mean of evaluations of regression function on testing data, if testing data is provided.}
\item{yhat.test}{If testing data was supplied, matrix containing posterior samples of the regression function evaluated on the testing data. Structure is similar to that of \code{yhat_train}. Only returned if testing data is passed and \code{save_samples == TRUE}.}
\item{sigma}{Vector containing ALL samples of the residual standard deviation, including burn-in.}
\item{varcounts}{Matrix that counts the number of times a variable was used in a decision rule in each MCMC iteration. Structure is similar to that of \code{yhat_train}, with rows corresponding to MCMC iteration and columns corresponding to predictors}
\item{trees}{A list (or length \code{nd}) of character vectors (of lenght \code{M}) containing textual representations of the regression trees. These strings are parsed by \code{predict_flexBART} to reconstruct the C++ representations of the sampled trees.}
}
