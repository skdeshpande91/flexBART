# Draw samples from the regression tree prior

Draws M binary regression trees from the flexBART tree prior, resulting
in a single prior draw from the sum-of-trees prior.

## Usage

``` r
rflexBART(train_data,nd,verbose = TRUE,
          print_every = floor(nd/10),...)
```

## Arguments

- train_data:

  an object of class `data.frame` containing data used to train the
  model. As usual, rows (resp. columns) correspond to observations
  (resp. variables)

- nd:

  Number of trees to draw

- verbose:

  Logical, indicating whether a message should be printed to the R
  console after every `print_every` trees have been sampled. Default is
  `FALSE`.

- print_every:

  If `verbose == TRUE`, then a message is printed to the R console after
  every `print_every` trees have been sampled. Default is
  `floor(nd/10)`.

- ...:

  Additional arguments for setting prior hyperparameters (e.g., number
  of trees, \\\mu\_{0}\\, \\\tau\\, etc.). See
  [`flexBART()`](https://skdeshpande91.github.io/flexBART/reference/flexBART.md)
  for details about additional arguments.

## Details

This function is useful for drawing samples from the regression tree
prior underpinning flexBART. Together, these sampled trees form a single
ensemble of trees. The main utility of this function is to study how
often certain observations are clustered together in individual trees.
This is key to understanding how flexBART “borrows strength” across
observations.

## Value

A list containing

- num_leafs:

  An integer vector of length `M` recording the number of leaf nodes in
  each tree of the ensemble.

- num_singletons:

  An integer vector of length `M` recording the number of leaf nodes in
  each tree that contain only one observation.

- num_empty:

  An integer vector of length `M` recording the number of leaf nodes in
  each tree that contain no observations.

- max_leaf_size:

  An integer vector of length `M` recording the maximum number of
  observation contained in a leaf in each tree of the ensemble.

- min_leaf_size:

  An integer vector of length `M` recording the minimum number of
  observation contained in a leaf in each tree of the ensemble.

- kernel:

  An n x n matrix whose (i,j) entry is the proportion of tree draws in
  which observations i and j land in the same leaf of the tree.
