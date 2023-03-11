\name{rnetwork_BART}
\alias{rnetwork_BART}
\title{
Draw a single sample from the ensemble of trees prior over a network
}
\description{
Draws M binary regression trees from the network_BART tree prior, resulting in a single prior draw from the sum-of-trees prior.
}
\usage{
rnetwork_BART(M, vertex_id,
              X_cont = matrix(0, nrow = 1, ncol = 1),
              unif_cuts = rep(TRUE, times = ncol(X_cont)),
              cutpoints_list = NULL,
              A = matrix(0, nrow = 1, ncol = 1),
              graph_split = TRUE,
              graph_cut_type = 0,
              alpha = 0.95, beta = 2, mu0 = 0, tau = 1/sqrt(M),
              verbose = TRUE, print_every = floor(M/10))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{Number of trees to draw}
    \item{vertex_id}{An integer vector of the same length as \code{Y_train} that records the vertex label for each observation. Assumed to take values in the set {1, 2, ..., n_vertex} where n_vertex is the number of vertices in the network}
  \item{X_cont}{Matrix of continuous predictors. Note all predictors must be re-scaled to lie in the interval [-1,1]. Default is a 1x1 matrix, which signals that no continuos predictors are available.}
  \item{unif_cuts}{Vector of logical values indicating whether cutpoints for each continuous predictor should be drawn from a continuous uniform distribution (\code{TRUE}) or a discrete set (\code{FALSE}) specified in \code{cutpoints_list}. Default is \code{TRUE} for each variable in \code{X_cont_train}}
  \item{cutpoints_list}{List of length \code{ncol(X_cont_train)} containing a vector of cutpoints for each continuous predictor. By default, this is set to \code{NULL} so that cutpoints are drawn uniformly from a continuous distribution.}
    \item{A}{Binary adjacency matrix of the network. Assumed to be symmetric.}
  \item{graph_split}{Logical, indicating whether to constrain the regression tree prior to respect adjacency when splitting on vertex labels. Default is \code{TRUE}.}
  \item{graph_cut_type}{An integer (0,1,2,3,4) that determines how partitions of the network are constructed. See Details.}
  \item{alpha}{In the decision tree prior, the probability that a node at depth d is non-terminal is \code{alpha * (1 + d)^(-beta)}. Default is 0.95}
  \item{beta}{In the decision tree prior, the probability that node at depth d is non-terminal is \code{alpha * ( 1 + d)^(-beta). Default is 2}}
  \item{mu0}{Prior mean for the jumps (i.e. the leaf parameters) in each tree. Default is 0.}
  \item{tau}{Prior standard deviation for the jumps (i.e. the leaf parameters) in each tree. Default is \code{1/sqrt(M)}$.}
  \item{verbose}{Logical, indicating whether a message should be printed to the R console after every \code{print_every} trees have been sampled. Default is \code{FALSE}.}
  \item{print_every}{If \code{verbose == TRUE}, then a message is printed to the R console after every \code{print_every} trees have been sampled. Default is \code{floor(M/10)}.}
}
\details{
This function is useful for drawing a M samples from the regression tree prior underpinning network_BART. Together, these M sampled trees form a single ensemble of trees.
The main utility of this function is to study how often certain observations are clustered together in individual trees. 
This is key to understanding how network_BART ``borrows strength'' across the vertices of a network.

To split on the vertex label, we randomly draw a partition of the network into two (possibly disconnected) subgraphs.
When \code{graph_cut_type == 1}, the partition is formed by first drawing a uniform random spanning tree using Wilson's algorithm and then deleting a uniformly selected edge from the tree.
When \code{graph_cut_type == 2}, the partition is formed by drawing a uniform random spanning tree using Wilson's algorithm, computing the eigenvectors of the unweighted graph Laplacian, and identifying the positive entries in the eigenvector corresponding to the second smallest eigenvalue.
When \code{graph_cut_type == 3}, the partition is formed by randomly selecting a single vertex and taking all vertices within a randomly drawn graph distance of the selected vertex. This procedure can result in disconnected subgraphs.
When \code{graph_cut_type == 4}, the partition is deterministic. It is formed by looking at the signs of entries in the second smallest eigenvector.
Finally, when \code{graph_cut_type == 0}, we uniformly select one of the above four partitioning procedures. 

}
\value{
A list containing the following elements
\item{fit}{A vector containing the value of sampled regression function evaluated at all of the supplied inputs}
\item{trees}{A character vector (of length \code{M}) containing text representations of the sampled trees. Currently this vector cannot be passed immediately to \code{predict_flexBART}}
\item{tree_fits}{A matrix with \code{M} columns containing the evaluations of each sampled regression trees. Note that \code{rowSums(tree_fits)} will be identical to \code{fit}.}
\item{leaf}{An integer matrix with recording the ID number of the leaf to which each observation is assigned. Nodes in a regression tree can be numbered as follows: the root is labelled with 1 and the left and right children of a node labelled x are labelled 2x and 2x+1. The rows of this matrix correspond to the rows of \code{X_cont} and \code{X_cat} and columns correspond to the individual tree.}
\item{kernel}{An n x n matrix whose (i,j) entry is the proportion of tree draws in which observations i and j land in the same leaf of the tree}
}