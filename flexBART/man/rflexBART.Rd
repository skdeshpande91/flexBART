\name{rflexBART}
\alias{rflexBART}
\title{
Draw samples from the regression tree prior
}
\description{
Draws M binary regression trees from the flexBART tree prior, resulting in a single prior draw from the sum-of-trees prior.
}
\usage{
rflexBART(train_data,nd,verbose = TRUE,
          print_every = floor(nd/10),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{train_data}{an object of class \code{data.frame} containing data used to train the model. 
  As usual, rows (resp. columns) correspond to observations (resp. variables)}
  \item{nd}{Number of trees to draw}
  \item{verbose}{Logical, indicating whether a message should be printed to the R console after every \code{print_every} trees have been sampled. Default is \code{FALSE}.}
  \item{print_every}{If \code{verbose == TRUE}, then a message is printed to the R console after every \code{print_every} trees have been sampled. Default is \code{floor(nd/10)}.}
    \item{\dots}{additional arguments to be passed to the low-level tree sampling functions (see below)}

}
\details{
This function is useful for drawing a M samples from the regression tree prior underpinning flexBART. Together, these M sampled trees form a single ensemble of trees.
The main utility of this function is to study how often certain observations are clustered together in individual trees. 
This is key to understanding how flexBART ``borrows strenght'' across observations. 
}
\value{
A list containing the following elements
\item{num_leafs}{An integer vector of length \code{M} recording the number of leaf nodes in each tree of the ensemble.}
\item{num_singletons}{An integer vector of length \code{M} recording the number of leaf nodes in each tree that contain only one observation.}
\item{num_empty}{An integer vector of length \code{M} recording the number of leaf nodes in each tree that contain no observations.}
\item{max_leaf_size}{An integer vector of length \code{M} recording the maximum number of observation contained in a leaf in each tree of the ensemble.}
\item{min_leaf_size}{An integer vector of length \code{M} recording the minimum number of observation contained in a leaf in each tree of the ensemble.}
\item{kernel}{An n x n matrix whose (i,j) entry is the proportion of tree draws in which observations i and j land in the same leaf of the tree.}
}