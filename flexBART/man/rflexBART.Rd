\name{rflexBART}
\alias{rflexBART}
\title{
Draw a single sample from the ensemble of tree prior
}
\description{
Draws M binary regression trees from the flexBART tree prior, resulting in a single prior draw from the sum-of-trees prior.
}
\usage{
rflexBART(M, 
          X_cont = matrix(0, nrow = 1, ncol = 1),
          X_cat = matrix(0L, nrow = 1, ncol = 1),
          unif_cuts = rep(TRUE, times = ncol(X_cont)),
          cutpoints_list = NULL,
          cat_levels_list = NULL,
          alpha = 0.95, beta = 2, mu0 = 0, tau = 1/sqrt(M),
          verbose = TRUE, print_every = floor(M/10))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{Number of trees to draw}
  \item{X_cont}{Matrix of continuous predictors. Note all predictors must be re-scaled to lie in the interval [-1,1]. Default is a 1x1 matrix, which signals that no continuos predictors are available.}
  \item{X_cat}{Integer matrix of categorical predictors for training data. Note categorical levels should be 0-indexed. That is, if a categorical predictor has 10 levels, the values should run from 0 to 9. Default is a 1x1 matrix, which signals that no categorical predictors.}
  \item{unif_cuts}{Vector of logical values indicating whether cutpoints for each continuous predictor should be drawn from a continuous uniform distribution (\code{TRUE}) or a discrete set (\code{FALSE}) specified in \code{cutpoints_list}. Default is \code{TRUE} for each variable in \code{X_cont_train}}
  \item{cutpoints_list}{List of length \code{ncol(X_cont_train)} containing a vector of cutpoints for each continuous predictor. By default, this is set to \code{NULL} so that cutpoints are drawn uniformly from a continuous distribution.}
  \item{cat_levels_list}{List of length \code{ncol(X_cat_train)} containing a vector of levels for each categorical predictor. If the j-th categorical predictor contains L levels, \code{cat_levels_list[[j]]} should be the vector \code{0:(L-1)}. Default is \code{NULL}, which corresponds to the case that no categorical predictors are available.}
  \item{alpha}{In the decision tree prior, the probability that a node at depth d is non-terminal is \code{alpha * (1 + d)^(-beta)}. Default is 0.95}
  \item{beta}{In the decision tree prior, the probability that node at depth d is non-terminal is \code{alpha * ( 1 + d)^(-beta). Default is 2}}
  \item{mu0}{Prior mean for the jumps (i.e. the leaf parameters) in each tree. Default is 0.}
  \item{tau}{Prior standard deviation for the jumps (i.e. the leaf parameters) in each tree. Default is \code{1/sqrt(M)}$.}
  \item{verbose}{Logical, indicating whether a message should be printed to the R console after every \code{print_every} trees have been sampled. Default is \code{FALSE}.}
  \item{print_every}{If \code{verbose == TRUE}, then a message is printed to the R console after every \code{print_every} trees have been sampled. Default is \code{floor(M/10)}.}
}
\details{
This function is useful for drawing a M samples from the regression tree prior underpinning flexBART. Together, these M sampled trees form a single ensemble of trees.
The main utility of this function is to study how often certain observations are clustered together in individual trees. 
This is key to understanding how flexBART ``borrows strenght'' across observations. 
}
\value{
A list containing the following elements
\item{fit}{A vector containing the value of sampled regression function evaluated at all of the supplied inputs}
\item{trees}{A character vector (of length \code{M}) containing text representations of the sampled trees. Currently this vector cannot be passed immediately to \code{predict_flexBART}}
\item{tree_fits}{A matrix with \code{M} columns containing the evaluations of each sampled regression trees. Note that \code{rowSums(tree_fits)} will be identical to \code{fit}.}
\item{num_clusters}{An integer vector of length \code{M} recording the number of leaf nodes in each tree of the ensemble.}
\item{num_singletons}{An integer vector of length \code{M} recording the number of leaf nodes in each tree that contain only one observation.}
\item{num_empty}{An integer vector of length \code{M} recording the number of leaf nodes in each tree that contain no observations.}
\item{max_leaf_size}{An integer vector of length \code{M} recording the maximum number of observation contained in a leaf in each tree of the ensemble.}
\item{min_leaf_size}{An integer vector of length \code{M} recording the minimum number of observation contained in a leaf in each tree of the ensemble.}
\item{kernel}{An n x n matrix whose (i,j) entry is the proportion of tree draws in which observations i and j land in the same leaf of the tree.}
}